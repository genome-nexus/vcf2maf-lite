#!/usr/bin/env python3
"""maf2vcf - Reformat variants in a given MAF into generic VCFs with GT:AD:DP data if available.

A Python port of maf2vcf.pl (vcf2maf, Cyriac Kandoth, Apache-2.0) with the bugs of that script fixed. 
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, OrderedDict
from datetime import datetime

VERSION = "1.0.0"

DEFAULT_REF_FASTA = os.path.join(
    os.path.expanduser("~"), ".vep", "homo_sapiens", "112_GRCh37",
    "Homo_sapiens.GRCh37.dna.toplevel.fa.gz",
)

# A MAF header line is one whose first column is one of these (case insensitive)
HEADER_FIRST_COLS = frozenset(("hugo_symbol", "chromosome", "tumor_sample_barcode"))
REQUIRED_COLS = ("chromosome", "start_position", "reference_allele", "tumor_sample_barcode")

# MAF uses "-" (and occasionally "?", "0", ".") to mean "no bases here"
NULL_ALLELE_RE = re.compile(r"^[-?0.]+$")
# Accept the full IUPAC alphabet, not just ACGT: N and ambiguity codes are legal in a VCF
IUPAC_RE = re.compile(r"^[ACGTUNRYSWKMBDHV]*$")
INT_RE = re.compile(r"^\d+$")
QUAL_RE = re.compile(r"^\d+(\.\d+)?([eE][-+]?\d+)?$")
# Characters that may not appear in a VCF FILTER/ID string
FILTER_BAD_CHARS_RE = re.compile(r"[\s;,=<>\"]+")
# Barcodes end up in filenames, so refuse anything that could escape the output directory
UNSAFE_BARCODE_RE = re.compile(r"[/\\\x00-\x1f]")

# samtools runs faster when passed many loci at a time, but the argument list is bounded
FAIDX_CHUNK_SIZE = 4000


class Maf2VcfError(Exception):
    """A fatal, user-facing error. Reported without a traceback."""


class Reporter:
    """Collects warnings/skips so a 10M-line MAF cannot produce 10M lines of stderr."""

    MAX_PER_KIND = 3

    def __init__(self):
        self.warn_counts = Counter()
        self.skip_counts = Counter()

    def warn(self, kind, detail=""):
        self.warn_counts[kind] += 1
        if self.warn_counts[kind] <= self.MAX_PER_KIND:
            suffix = f": {detail}" if detail else ""
            sys.stderr.write(f"WARNING: {kind}{suffix}\n")
        elif self.warn_counts[kind] == self.MAX_PER_KIND + 1:
            sys.stderr.write(f"WARNING: {kind} (further occurrences will only be counted)\n")

    def skip(self, kind):
        self.skip_counts[kind] += 1

    def summarize(self):
        for kind, count in sorted(self.warn_counts.items()):
            if count > self.MAX_PER_KIND:
                sys.stderr.write(f"WARNING: {kind} occurred {count} times\n")
        for kind, count in sorted(self.skip_counts.items()):
            sys.stderr.write(f"SKIPPED: {count} variant(s) - {kind}\n")


# --------------------------------------------------------------------------------------
# Reference FASTA access
# --------------------------------------------------------------------------------------

class Reference:
    """Random-access reader for an (optionally bgzipped) reference FASTA.

    Prefers pysam, which needs no external binary and needs no batching. Falls back to
    invoking `samtools faidx` with a real argument list (never a shell string), fetching
    windows in batches because process startup dominates otherwise.
    """

    def __init__(self, path, samtools=None):
        self.path = path
        self.fai_path = path + ".fai"
        self._pysam_fa = None
        self._samtools = samtools
        self._cache = {}
        self._prefetched = None

        try:
            import pysam  # noqa: F401  (optional dependency)
            self._pysam = pysam
        except ImportError:
            self._pysam = None

        if self._pysam is not None:
            if not os.path.exists(self.fai_path):
                try:
                    self._pysam.faidx(path)
                except Exception as err:  # pysam raises a family of subclassed errors
                    raise Maf2VcfError(f"Failed to index reference FASTA {path}: {err}")
            self._pysam_fa = self._pysam.FastaFile(path)
        else:
            if self._samtools is None:
                raise Maf2VcfError(
                    "Neither the pysam module nor a samtools binary is available.\n"
                    "Install one of them: `pip install pysam` or put samtools in your PATH."
                )
            if not os.path.exists(self.fai_path):
                run_samtools([self._samtools, "faidx", path])

        self.contig_lengths = self._load_fai()
        if not self.contig_lengths:
            raise Maf2VcfError(f"Reference FASTA index is empty or unreadable: {self.fai_path}")
        # Contig order in the .fai is the canonical order for sorting the output VCF
        self.contig_rank = {name: i for i, name in enumerate(self.contig_lengths)}

    @property
    def random_access(self):
        return self._pysam_fa is not None

    def _load_fai(self):
        lengths = OrderedDict()
        try:
            with open(self.fai_path) as fh:
                for line in fh:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 2 and INT_RE.match(fields[1]):
                        lengths[fields[0]] = int(fields[1])
        except OSError as err:
            raise Maf2VcfError(f"Couldn't read reference FASTA index {self.fai_path}: {err}")
        return lengths

    def prefetch(self, windows):
        """Batch-fetch (chrom, start, end) windows. A no-op when random access is available."""
        if self.random_access or not windows:
            return
        self._prefetched = {}
        windows = sorted(set(windows))
        regions = [f"{c}:{s}-{e}" for c, s, e in windows]
        for i in range(0, len(regions), FAIDX_CHUNK_SIZE):
            chunk = regions[i:i + FAIDX_CHUNK_SIZE]
            stdout = run_samtools([self._samtools, "faidx", self.path] + chunk)
            self._prefetched.update(parse_fasta(stdout))

    def get(self, chrom, start, end):
        """Return uppercased bases for the 1-based inclusive window, or None if unavailable."""
        if start > end:
            return ""
        if self.random_access:
            key = (chrom, start, end)
            seq = self._cache.get(key)
            if seq is None:
                try:
                    seq = self._pysam_fa.fetch(chrom, start - 1, end).upper()
                except (KeyError, ValueError):
                    return None
                if len(self._cache) > 100000:  # bound the cache on pathological MAFs
                    self._cache.clear()
                self._cache[key] = seq
            return seq or None
        seq = (self._prefetched or {}).get(f"{chrom}:{start}-{end}")
        return seq or None

    def close(self):
        if self._pysam_fa is not None:
            self._pysam_fa.close()


def run_samtools(argv):
    """Run samtools with an argument list (no shell) and fail loudly on a non-zero exit."""
    try:
        proc = subprocess.run(argv, capture_output=True, text=True)
    except OSError as err:
        raise Maf2VcfError(f"Failed to run {argv[0]}: {err}")
    if proc.returncode != 0:
        stderr = proc.stderr.strip().splitlines()
        detail = stderr[-1] if stderr else f"exit status {proc.returncode}"
        raise Maf2VcfError(f"samtools faidx failed: {detail}")
    return proc.stdout


def parse_fasta(text):
    """Parse multi-FASTA text into {header: uppercased sequence}."""
    seqs = {}
    header = None
    chunks = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                seqs[header] = "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        seqs[header] = "".join(chunks).upper()
    return seqs


# --------------------------------------------------------------------------------------
# MAF reading
# --------------------------------------------------------------------------------------

def open_maf(path):
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    try:
        # newline="" disables translation so we can detect bare-CR (Mac OS 9) line breaks
        return opener(path, "rt", newline="", errors="replace")
    except OSError as err:
        raise Maf2VcfError(f"Couldn't open input MAF {path}: {err}")


def iter_maf_lines(path):
    """Yield (line_number, line) with line terminators stripped."""
    with open_maf(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            # With newline="", a CRLF line ends in "\r\n" and a bare-CR line ends in "\r"
            if raw.endswith("\r"):
                raise Maf2VcfError(
                    "Your MAF uses CR line breaks, which we can't support. Please use LF or CRLF."
                )
            yield lineno, raw.rstrip("\r\n")


def split_fields(line):
    """Split on tabs, keeping trailing empty fields, and trim whitespace from each."""
    return [field.strip() for field in line.split("\t")]


def get_col(cols, col_idx, name):
    """Value of a MAF column, or "" when the column is absent or the row is short."""
    idx = col_idx.get(name)
    if idx is None or idx >= len(cols):
        return ""
    return cols[idx]


def parse_header(cols, reporter):
    col_idx = {}
    for idx, name in enumerate(cols):
        key = name.strip().lower()
        if not key:
            continue
        if key in col_idx:
            reporter.warn("Duplicate MAF column name; using the leftmost", name)
            continue
        col_idx[key] = idx

    missing = [c for c in REQUIRED_COLS if c not in col_idx]
    if missing:
        raise Maf2VcfError(
            "Missing required MAF column(s): " + ", ".join(missing)
        )
    if "tumor_seq_allele1" not in col_idx and "tumor_seq_allele2" not in col_idx:
        raise Maf2VcfError("At least one Tumor_Seq_Allele column must be defined!")
    return col_idx


# --------------------------------------------------------------------------------------
# Field normalization
# --------------------------------------------------------------------------------------

def clean_allele(value):
    """Uppercase, drop Oncotator's "|" DNP separators, and blank out MAF null markers."""
    allele = value.replace("|", "").strip().upper()
    if NULL_ALLELE_RE.match(allele):
        return ""
    return allele


def parse_depth(value):
    """Return a non-negative integer as a string, or "." when unusable."""
    text = value.strip()
    if not text or text == ".":
        return "."
    if INT_RE.match(text):
        return str(int(text))
    try:
        number = float(text)
    except ValueError:
        return "."
    if number != number or number < 0:  # NaN or negative
        return "."
    return str(int(round(number)))


def sanitize_filter_tags(value):
    """Split a MAF FILTER value into clean VCF-legal tags."""
    tags = []
    for raw in re.split(r"[,;]", value):
        tag = FILTER_BAD_CHARS_RE.sub("_", raw.strip()).strip("_")
        if tag and tag != "." and tag not in tags:
            tags.append(tag)
    return tags


def format_filter(tags):
    if not tags:
        return "."
    if tags == ["PASS"]:
        return "PASS"
    # PASS is not meaningful alongside other tags
    tags = [t for t in tags if t != "PASS"] or ["PASS"]
    return ";".join(tags)  # VCF delimits FILTER with ";", never ","


def check_barcode(barcode, kind):
    if not barcode or barcode.startswith("#") or UNSAFE_BARCODE_RE.search(barcode):
        raise Maf2VcfError(f'Invalid {kind} in MAF: "{barcode}"')
    return barcode


# --------------------------------------------------------------------------------------
# Pass 1: scan the MAF for its header, sample pairs, FILTER tags, and loci
# --------------------------------------------------------------------------------------

class MafScan:
    def __init__(self):
        self.header_line = None
        self.col_idx = None
        self.pairs = []          # sorted, unique (tumor_id, normal_id)
        self.filter_tags = set()
        self.chroms = set()
        self.windows = set()     # only populated when the reference needs batch prefetching
        self.variant_lines = 0


def scan_maf(path, need_windows, reporter):
    scan = MafScan()
    pairs = set()

    for lineno, line in iter_maf_lines(path):
        if not line or line.startswith("#"):
            continue
        cols = split_fields(line)

        if cols[0].lower() in HEADER_FIRST_COLS:
            if scan.col_idx is None:
                scan.header_line = line
                scan.col_idx = parse_header(cols, reporter)
            elif line != scan.header_line:
                raise Maf2VcfError(
                    f"MAF line {lineno} looks like a second, different header line. "
                    "Concatenated MAFs with mismatched columns are not supported."
                )
            continue

        if scan.col_idx is None:
            raise Maf2VcfError(
                "Couldn't find a header line (must start with Hugo_Symbol, Chromosome, "
                f"or Tumor_Sample_Barcode): {path}"
            )

        scan.variant_lines += 1
        col_idx = scan.col_idx

        chrom = get_col(cols, col_idx, "chromosome")
        pos_text = get_col(cols, col_idx, "start_position")
        ref = clean_allele(get_col(cols, col_idx, "reference_allele"))
        t_id = get_col(cols, col_idx, "tumor_sample_barcode")
        n_id = get_col(cols, col_idx, "matched_norm_sample_barcode")

        # These four are structural: without them the row cannot be placed in a VCF at all
        for name, value in (("Chromosome", chrom), ("Start_Position", pos_text),
                            ("Reference_Allele", get_col(cols, col_idx, "reference_allele")),
                            ("Tumor_Sample_Barcode", t_id)):
            if value == "":
                raise Maf2VcfError(f"{name} is empty in MAF line {lineno}!")
        if not INT_RE.match(pos_text) or int(pos_text) < 1:
            raise Maf2VcfError(
                f"Start_Position is not a positive integer in MAF line {lineno}: {pos_text}"
            )

        check_barcode(t_id, "Tumor_Sample_Barcode")
        # An unmatched tumor gets its own placeholder so that several unmatched tumors
        # cannot collide into a single shared "NORMAL" genotype column
        n_id = check_barcode(n_id, "Matched_Norm_Sample_Barcode") if n_id else f"{t_id}_NORMAL"

        pairs.add((t_id, n_id))
        scan.chroms.add(chrom)
        scan.filter_tags.update(sanitize_filter_tags(get_col(cols, col_idx, "filter")))

        if need_windows:
            pos = int(pos_text)
            scan.windows.add(window_for(chrom, pos, len(ref)))

    if scan.col_idx is None:
        raise Maf2VcfError(f"Couldn't find a header line in {path} - is it really a MAF?")
    if scan.variant_lines == 0:
        raise Maf2VcfError("No variant lines in the input MAF!")

    scan.filter_tags.discard("PASS")
    scan.pairs = sorted(pairs)
    return scan


def window_for(chrom, pos, ref_len, contig_len=None):
    """1-based inclusive window covering REF plus one flanking base on each side."""
    start = max(1, pos - 1)
    end = pos + ref_len
    if contig_len is not None:
        end = min(end, contig_len)
    return (chrom, start, end)


# --------------------------------------------------------------------------------------
# Pass 2: convert each MAF row into VCF-ready data
# --------------------------------------------------------------------------------------

class Converter:
    def __init__(self, args, scan, reference, reporter):
        self.args = args
        self.scan = scan
        self.reference = reference
        self.reporter = reporter
        # key -> [id, qual, filter, {sample_id: "GT:AD:DP"}]
        self.records = {}
        # (tumor, normal) -> ordered set of keys, only tracked for --per-tn-vcfs
        self.pair_keys = {pair: {} for pair in scan.pairs} if args.per_tn_vcfs else None
        self.skipped_fh = None
        self.skipped_path = None
        self.converted = 0

    # -- output file for rows we could not convert -------------------------------------
    def _record_skip(self, reason, line):
        self.reporter.skip(reason)
        if self.skipped_fh is None:
            self.skipped_path = derive_output_path(
                self.args.output_dir, self.args.input_maf, ".skipped.tsv")
            self.skipped_fh = open(self.skipped_path, "w")
            self.skipped_fh.write(self.scan.header_line + "\n")
            sys.stderr.write(
                f"WARNING: Some variants could not be converted. "
                f"Storing them here for debugging: {self.skipped_path}\n")
        self.skipped_fh.write(line + "\n")

    def close(self):
        if self.skipped_fh is not None:
            self.skipped_fh.close()

    # -- main loop ---------------------------------------------------------------------
    def run(self):
        col_idx = self.scan.col_idx
        args = self.args

        for lineno, line in iter_maf_lines(args.input_maf):
            if not line or line.startswith("#"):
                continue
            cols = split_fields(line)
            if cols[0].lower() in HEADER_FIRST_COLS:
                continue

            skip = self._convert_row(lineno, line, cols, col_idx)
            if skip:
                self._record_skip(skip, line)
            else:
                self.converted += 1

    def _convert_row(self, lineno, line, cols, col_idx):
        """Convert one MAF row. Returns None on success, or a skip reason."""
        args = self.args
        chrom = get_col(cols, col_idx, "chromosome")
        pos = int(get_col(cols, col_idx, "start_position"))
        ref = clean_allele(get_col(cols, col_idx, "reference_allele"))

        t_id = get_col(cols, col_idx, "tumor_sample_barcode")
        n_id = get_col(cols, col_idx, "matched_norm_sample_barcode") or f"{t_id}_NORMAL"

        t_al1_raw = get_col(cols, col_idx, "tumor_seq_allele1")
        t_al2_raw = get_col(cols, col_idx, "tumor_seq_allele2")
        n_al1_raw = get_col(cols, col_idx, "match_norm_seq_allele1")
        n_al2_raw = get_col(cols, col_idx, "match_norm_seq_allele2")

        # Validate the actual allele values (the Perl original checked the column *indexes*)
        if t_al1_raw.strip() == "" and t_al2_raw.strip() == "":
            return "no Tumor_Seq_Allele value"

        t_al1 = clean_allele(t_al1_raw)
        t_al2 = clean_allele(t_al2_raw)
        n_al1 = clean_allele(n_al1_raw)
        n_al2 = clean_allele(n_al2_raw)

        # When variant alleles are a SNP and a "-", warn the user of misusing "-" as REF
        if (len(ref) == 1 and t_al1_raw.strip() == "-" and len(t_al2) == 1 and t_al2 != ref
                and t_al1 == ""):
            t_al1 = ref
            self.reporter.warn("Replacing '-' with the reference allele", f"MAF line {lineno}")
        elif (len(ref) == 1 and t_al2_raw.strip() == "-" and len(t_al1) == 1 and t_al1 != ref
                and t_al2 == ""):
            t_al2 = ref
            self.reporter.warn("Replacing '-' with the reference allele", f"MAF line {lineno}")

        # If one tumor allele is unset, assume it matches the reference allele
        if t_al1_raw.strip() == "":
            t_al1 = ref
        if t_al2_raw.strip() == "":
            t_al2 = ref

        # Whether the normal was genotyped at all decides between 0/0 and ./. below
        normal_known = n_al1_raw.strip() != "" or n_al2_raw.strip() != ""
        if not normal_known:
            n_al1 = n_al2 = ref
        else:
            if n_al1_raw.strip() == "":
                n_al1 = ref
            if n_al2_raw.strip() == "":
                n_al2 = ref

        for allele in (ref, t_al1, t_al2, n_al1, n_al2):
            if not IUPAC_RE.match(allele):
                return "invalid characters in an allele"

        # Ensure the second allele is always the non-REF one, so it lands at ALT index 1
        # where the MAF's t_alt_count/n_alt_count read depths belong
        if t_al2 == ref and t_al1 != ref:
            t_al1, t_al2 = t_al2, t_al1
        if n_al2 == ref and n_al1 != ref:
            n_al1, n_al2 = n_al2, n_al1

        # -- reference check ------------------------------------------------------------
        contig_len = self.reference.contig_lengths.get(chrom)
        if contig_len is None:
            return "chromosome not in the reference FASTA"
        if pos > contig_len:
            return "Start_Position is past the end of the contig"

        win_chrom, win_start, win_end = window_for(chrom, pos, len(ref), contig_len)
        seq = self.reference.get(win_chrom, win_start, win_end)
        if seq is None or len(seq) != win_end - win_start + 1:
            return "could not fetch reference bases"

        def base_at(position):
            if position < win_start or position > win_end:
                return None
            return seq[position - win_start]

        if ref:
            if pos + len(ref) - 1 > win_end:
                return "REF allele extends past the end of the contig"
            observed = seq[pos - win_start: pos - win_start + len(ref)]
            if observed != ref:
                return "reference allele mismatch"

        self._check_end_position(cols, col_idx, pos, len(ref), lineno)

        # -- represent indels the way VCF wants, with an anchoring base -------------------
        lengths = (len(ref), len(t_al1), len(t_al2), len(n_al1), len(n_al2))
        needs_anchor = (
            min(lengths) == 0
            or (len(ref) != len(t_al2) and ref[0] != t_al2[0])
        )
        if needs_anchor:
            if not ref:
                # For MAF-format simple insertions, pos is already the preceding base
                anchor = base_at(pos)
                if anchor is None:
                    return "could not fetch an anchoring base for an insertion"
                ref, t_al1, t_al2, n_al1, n_al2 = (
                    anchor + a for a in (ref, t_al1, t_al2, n_al1, n_al2))
            else:
                anchor = base_at(pos - 1)
                if anchor is not None:
                    pos -= 1
                    ref, t_al1, t_al2, n_al1, n_al2 = (
                        anchor + a for a in (ref, t_al1, t_al2, n_al1, n_al2))
                else:
                    # At the very start of a contig there is no preceding base, so VCF
                    # allows anchoring on the base that follows the event instead
                    suffix = base_at(pos + len(ref))
                    if suffix is None:
                        return "no flanking base available to anchor an indel"
                    ref, t_al1, t_al2, n_al1, n_al2 = (
                        a + suffix for a in (ref, t_al1, t_al2, n_al1, n_al2))

        # -- allele indexing, exactly as a VCF numbers them ------------------------------
        alleles = []
        al_idx = {}
        for allele in (ref, t_al2, t_al1, n_al2, n_al1):
            if allele not in al_idx:
                al_idx[allele] = len(alleles)
                alleles.append(allele)
        if len(alleles) == 1:
            return "no variant allele (all alleles match the reference)"
        alt = ",".join(alleles[1:])

        # -- read depths -----------------------------------------------------------------
        t_dp = parse_depth(get_col(cols, col_idx, args.tum_depth_col.lower()))
        t_rad = parse_depth(get_col(cols, col_idx, args.tum_rad_col.lower()))
        t_vad = parse_depth(get_col(cols, col_idx, args.tum_vad_col.lower()))
        n_dp = parse_depth(get_col(cols, col_idx, args.nrm_depth_col.lower()))
        n_rad = parse_depth(get_col(cols, col_idx, args.nrm_rad_col.lower()))
        n_vad = parse_depth(get_col(cols, col_idx, args.nrm_vad_col.lower()))

        # AD is Number=R: REF depth, then one value per ALT. Only the first ALT (the MAF's
        # Tumor_Seq_Allele2) has a known depth; any others are unknown.
        pad = ["."] * (len(alleles) - 2)
        t_ad = "." if t_rad == t_vad == "." else ",".join([t_rad, t_vad] + pad)
        n_ad = "." if n_rad == n_vad == "." else ",".join([n_rad, n_vad] + pad)

        # -- genotypes -------------------------------------------------------------------
        t_gt = "/".join(str(i) for i in sorted((al_idx[t_al1], al_idx[t_al2])))
        if normal_known:
            n_gt = "/".join(str(i) for i in sorted((al_idx[n_al1], al_idx[n_al2])))
        elif n_dp != "." or n_rad != "." or n_vad != ".":
            # The normal was sequenced but not genotyped in the MAF; hom-ref is implied
            n_gt = "0/0"
        else:
            # Nothing is known about the normal, so do not assert hom-ref
            n_gt = "./."

        t_fmt = f"{t_gt}:{t_ad}:{t_dp}"
        n_fmt = f"{n_gt}:{n_ad}:{n_dp}"

        # -- ID / QUAL / FILTER ----------------------------------------------------------
        vid = get_col(cols, col_idx, "variant_id").strip() or "."
        vid = FILTER_BAD_CHARS_RE.sub("_", vid) if vid != "." else "."
        qual_raw = get_col(cols, col_idx, "variant_qual").strip()
        qual = qual_raw if QUAL_RE.match(qual_raw) else "."
        filt = format_filter(sanitize_filter_tags(get_col(cols, col_idx, "filter")))

        # -- store -----------------------------------------------------------------------
        key = (chrom, pos, ref, alt)
        record = self.records.get(key)
        if record is None:
            record = self.records[key] = [vid, qual, filt, {}]
        else:
            # Samples shouldn't disagree on ID/QUAL/FILTER for the same locus+alleles
            for i, (new, field) in enumerate(((vid, "ID"), (qual, "QUAL"), (filt, "FILTER"))):
                old = record[i]
                if new == "." or new == old:
                    continue          # nothing new to say
                elif old == ".":
                    record[i] = new   # fill in what the earlier row left blank
                else:
                    self.reporter.warn(
                        f"Conflicting {field} values for the same variant; keeping the first",
                        f"{chrom}:{pos} {ref}>{alt}")

        fmt_map = record[3]
        for sample, fmt in ((t_id, t_fmt), (n_id, n_fmt)):
            if sample in fmt_map and fmt_map[sample] != fmt:
                self.reporter.warn(
                    "Duplicate MAF rows with conflicting genotypes; keeping the first",
                    f"{sample} at {chrom}:{pos}")
            else:
                fmt_map[sample] = fmt

        if self.pair_keys is not None:
            self.pair_keys.setdefault((t_id, n_id), {})[key] = None
        return None

    def _check_end_position(self, cols, col_idx, pos, ref_len, lineno):
        end_text = get_col(cols, col_idx, "end_position")
        if not INT_RE.match(end_text):
            return
        expected = pos + 1 if ref_len == 0 else pos + ref_len - 1
        if int(end_text) != expected:
            self.reporter.warn(
                "End_Position disagrees with Start_Position and Reference_Allele",
                f"MAF line {lineno} (expected {expected}, found {end_text})")


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------

def derive_output_path(output_dir, input_maf, suffix):
    """<output_dir>/<input basename with one known extension stripped><suffix>."""
    name = os.path.basename(input_maf)
    name = re.sub(r"\.(gz|bgz)$", "", name, flags=re.IGNORECASE)
    # Anchored on a literal dot, so "mymaf" no longer becomes "my"
    name = re.sub(r"\.(maf|tsv|txt)$", "", name, flags=re.IGNORECASE)
    return os.path.join(output_dir, name + suffix)


def build_header(reference, filter_tags, sample_ids):
    lines = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now().strftime('%Y%m%d')}",
        f"##source=maf2vcf.py v{VERSION}",
        f"##reference=file://{os.path.abspath(reference.path)}",
    ]
    lines += [f"##contig=<ID={name},length={length}>"
              for name, length in reference.contig_lengths.items()]
    lines += [
        '##FILTER=<ID=%s,Description="FILTER tag carried over from the input MAF">' % tag
        for tag in sorted(filter_tags)
    ]
    lines += [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,'
        'Description="Allelic depths of REF and ALT(s) in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth across this site">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids),
    ]
    return "\n".join(lines) + "\n"


def write_vcf(path, header, keys, records, sample_ids):
    with open(path, "w") as fh:
        fh.write(header)
        for key in keys:
            chrom, pos, ref, alt = key
            vid, qual, filt, fmt_map = records[key]
            fields = [chrom, str(pos), vid, ref, alt, qual, filt, ".", "GT:AD:DP"]
            fields += [fmt_map.get(sample, "./.:.:.") for sample in sample_ids]
            fh.write("\t".join(fields) + "\n")


def sort_keys(keys, contig_rank):
    # Coordinate-sorted output, using the reference's own contig order
    return sorted(keys, key=lambda k: (contig_rank.get(k[0], len(contig_rank)), k[0], k[1], k[2], k[3]))


def write_pairs_tsv(path, pairs):
    with open(path, "w") as fh:
        fh.write("#Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\n")
        for t_id, n_id in pairs:
            fh.write(f"{t_id}\t{n_id}\n")


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="maf2vcf.py",
        description="Reformat variants in a MAF into a sorted multi-sample VCF with "
                    "GT:AD:DP data if available.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input-maf", required=True,
                        help="Path to input file in MAF format (plain or gzipped)")
    parser.add_argument("--output-dir", required=True,
                        help="Path to output directory where VCFs will be stored")
    parser.add_argument("--output-vcf",
                        help="Path to output multi-sample VCF containing all TN-pairs "
                             "[<output-dir>/<input-maf-name>.vcf]")
    parser.add_argument("--ref-fasta", default=DEFAULT_REF_FASTA,
                        help="Path to reference Fasta file")
    parser.add_argument("--per-tn-vcfs", action="store_true",
                        help="Also generate one VCF per TN pair")
    parser.add_argument("--tum-depth-col", default="t_depth",
                        help="Name of MAF column for read depth in tumor BAM")
    parser.add_argument("--tum-rad-col", default="t_ref_count",
                        help="Name of MAF column for reference allele depth in tumor BAM")
    parser.add_argument("--tum-vad-col", default="t_alt_count",
                        help="Name of MAF column for variant allele depth in tumor BAM")
    parser.add_argument("--nrm-depth-col", default="n_depth",
                        help="Name of MAF column for read depth in normal BAM")
    parser.add_argument("--nrm-rad-col", default="n_ref_count",
                        help="Name of MAF column for reference allele depth in normal BAM")
    parser.add_argument("--nrm-vad-col", default="n_alt_count",
                        help="Name of MAF column for variant allele depth in normal BAM")
    parser.add_argument("--version", action="version", version=f"maf2vcf.py {VERSION}")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    reporter = Reporter()

    if not os.path.isfile(args.input_maf):
        raise Maf2VcfError(f"Input MAF is missing: {args.input_maf}")
    if not os.path.isfile(args.ref_fasta) or os.path.getsize(args.ref_fasta) == 0:
        raise Maf2VcfError(
            f"Provided Reference FASTA is missing or empty! Path: {args.ref_fasta}")

    reference = Reference(args.ref_fasta, samtools=shutil.which("samtools"))
    try:
        scan = scan_maf(args.input_maf, need_windows=not reference.random_access, reporter=reporter)

        # Catch a build/naming mismatch up front rather than skipping every single variant
        unknown = sorted(scan.chroms - set(reference.contig_lengths))
        if unknown:
            hint = ""
            prefixed = {"chr" + c for c in unknown} & set(reference.contig_lengths)
            stripped = {c[3:] for c in unknown if c.startswith("chr")} & set(reference.contig_lengths)
            if prefixed or stripped:
                hint = " The reference uses the opposite 'chr' prefix convention from your MAF."
            raise Maf2VcfError(
                "These MAF chromosomes are not in the reference FASTA: "
                + ", ".join(unknown[:10]) + ("..." if len(unknown) > 10 else "") + "."
                + hint
                + f"\nMake sure --ref-fasta is the same genome build as your MAF: {args.ref_fasta}")

        reference.prefetch(scan.windows)

        os.makedirs(args.output_dir, exist_ok=True)
        output_vcf = args.output_vcf or derive_output_path(
            args.output_dir, args.input_maf, ".vcf")
        os.makedirs(os.path.dirname(os.path.abspath(output_vcf)), exist_ok=True)

        converter = Converter(args, scan, reference, reporter)
        try:
            converter.run()
        finally:
            converter.close()

        if not converter.records:
            raise Maf2VcfError(
                "No variants could be converted."
                + (f" See {converter.skipped_path} for the rows that were skipped."
                   if converter.skipped_path else ""))

        # Sample column order: tumor then normal, following the sorted list of TN pairs
        sample_ids = []
        seen = set()
        for t_id, n_id in scan.pairs:
            for sample in (t_id, n_id):
                if sample not in seen:
                    seen.add(sample)
                    sample_ids.append(sample)

        write_pairs_tsv(derive_output_path(args.output_dir, args.input_maf, ".pairs.tsv"),
                        scan.pairs)

        sorted_keys = sort_keys(converter.records.keys(), reference.contig_rank)
        write_vcf(output_vcf, build_header(reference, scan.filter_tags, sample_ids),
                  sorted_keys, converter.records, sample_ids)

        if args.per_tn_vcfs:
            for (t_id, n_id), keys in converter.pair_keys.items():
                path = os.path.join(args.output_dir, f"{t_id}_vs_{n_id}.vcf")
                write_vcf(path,
                          build_header(reference, scan.filter_tags, [t_id, n_id]),
                          sort_keys(keys.keys(), reference.contig_rank),
                          converter.records, [t_id, n_id])
    finally:
        reference.close()

    reporter.summarize()
    sys.stderr.write(
        f"Converted {converter.converted} of {scan.variant_lines} MAF variant line(s) into "
        f"{len(converter.records)} VCF record(s) across {len(sample_ids)} sample(s): {output_vcf}\n")
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Maf2VcfError as error:
        sys.stderr.write(f"ERROR: {error}\n")
        sys.exit(1)
    except BrokenPipeError:
        sys.exit(0)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted\n")
        sys.exit(130)
