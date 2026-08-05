# maf2vcf.py

Convert variants from a MAF into a multi-sample VCF
with `GT:AD:DP` genotype fields. Optionally, the script can also
generate one VCF for each tumor-normal pair.

This is a Python port of `maf2vcf.pl` (from [vcf2maf](https://github.com/mskcc/vcf2maf)) with bugs of that script fixed. 

---

## Requirements

- Python 3.8 or newer
- A reference FASTA matching the MAF's genome build - see [Getting a reference
  FASTA](#getting-a-reference-fasta)
- `pysam` (`pip install pysam`) or `samtools`

The script picks pysam automatically or falls back to `samtools
faidx` otherwise. Both produce identical output.

You do **not** need to index the FASTA. The script does it on the first run.

---

## Quick start

```bash
cd maf2vcf/

python3 maf2vcf.py --input-maf test_mafs/mut_somatic.maf --output-dir vcfs --ref-fasta ref_genome/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
```

Add `--per-tn-vcfs` to also write one VCF per tumor-normal pair:

```bash
python3 maf2vcf.py --input-maf test_mafs/mut_somatic.maf --output-dir vcfs --ref-fasta ref_genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz --per-tn-vcfs
```
---

## Command-line arguments

Required:

| Argument | What it does |
|---|---|
| `--input-maf` | Input MAF. Can be gzipped. |
| `--output-dir` | Where output goes. Created if it doesn't exist. |

Optional:

| Argument | Default | What it does |
|---|---|---|
| `--output-vcf` | `<output-dir>/<maf name>.vcf` | Where to write the multi-sample VCF. |
| `--ref-fasta` | the reference file path | Reference FASTA. Uncompressed or bgzipped, not plain gzip. |
| `--per-tn-vcfs` | off | Also write one VCF per tumor-normal pair. |
| `--version` | | Print version and exit. |
| `-h`, `--help` | | Print usage and exit. |

Read-depth column names:

These tell the script which MAF columns hold read counts. Change them only if your MAF
uses non-standard names. Any column that's missing or unparseable becomes `.` in the VCF.

| Argument | Default | Feeds |
|---|---|---|
| `--tum-depth-col` | `t_depth` | Tumor `DP` |
| `--tum-rad-col` | `t_ref_count` | Tumor `AD`, REF element |
| `--tum-vad-col` | `t_alt_count` | Tumor `AD`, first ALT element |
| `--nrm-depth-col` | `n_depth` | Normal `DP` |
| `--nrm-rad-col` | `n_ref_count` | Normal `AD`, REF element |
| `--nrm-vad-col` | `n_alt_count` | Normal `AD`, first ALT element |

---

## Input MAF file requirements

The MAF must have a header line whose **first column** is `Hugo_Symbol`, `Chromosome`, or
`Tumor_Sample_Barcode` (case-insensitive). Lines starting with `#` are ignored.

Required columns: `Chromosome`, `Start_Position`, `Reference_Allele`,
`Tumor_Sample_Barcode`, and at least one `Tumor_Seq_Allele` column.

Used if present: `Matched_Norm_Sample_Barcode`, `Match_Norm_Seq_Allele1/2`,
`End_Position`, `FILTER`, `variant_id`, `variant_qual`, and the read count columns.

Column names are not case sensitive.

---

## Output files

For an input named `mut_somatic.maf` and `--output-dir vcfs`:

| File | Contents |
|---|---|
| `vcfs/mut_somatic.vcf` | The multi-sample VCF: one column per sample and coordinate-sorted. |
| `vcfs/mut_somatic.pairs.tsv` | The tumor-normal pairing, which a multi-sample VCF cannot represent. Two columns, one row per pair. |
| `vcfs/mut_somatic.skipped.tsv` | **Only created if rows were skipped.** Header plus the skipped rows as-is, so the file is itself a valid MAF you can fix and re-run. |
| `vcfs/<tumor>_vs_<normal>.vcf` | With `--per-tn-vcfs`: one two-sample VCF per pair. |

Messages go to stderr:
```
SKIPPED: 1 variant(s) - reference allele mismatch
Converted 12 of 14 MAF variant line(s) into 11 VCF record(s) across 6 sample(s): vcfs/x.vcf
```

## How the conversion works

**Records are merged across samples.** One VCF record per unique
`(CHROM, POS, REF, ALT)`; a variant seen in 5 samples becomes one record with 5 populated
genotype columns. The multi-sample VCF file record count can be lower than the MAF row count. Nothing is lost.

**Skipped rows go to `.skipped.tsv`.** If that file isn't there, nothing was skipped.

**Normals are `./.` when the MAF says nothing about them.** They're only `0/0` if the MAF
has normal alleles or normal read counts. The Perl script always said `0/0`.

**Duplicate rows: the first one is selected.** If a sample shows up twice at the same variant
with different counts, the script keeps the first and warns. This happens when one normal
is paired with several tumors and the counts don't agree. The Perl script kept the last.

**Unmatched tumors get their own normal column,** named `<tumor>_NORMAL`. The Perl script
merged them all in one shared `NORMAL` column, where they overwrote each other.

---

## Getting a reference FASTA

The reference FASTA must match both the genome build and chromosome naming used in the MAF. Check the MAF first:

```bash
cut -f<NCBI_Build col> your.maf | sort -u | head    # GRCh37 or GRCh38?
cut -f<Chromosome col> your.maf | sort -u | head -3 # 1 or chr1?
```

#### GRCh37 (no `chr` prefix)

```bash
curl -O https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
```

> **Note:** Ensembl FASTA files are plain gzip, so unzip them before using with `samtools` or `pysam`. Recompressing with `bgzip` is optional.

#### GRCh38 (no `chr` prefix)

These files are already bgzipped and indexed.

```bash
base=https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
curl -O $base && curl -O $base.fai && curl -O $base.gzi
```
