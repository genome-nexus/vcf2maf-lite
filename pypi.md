# vcf2maf_lite

vcf2maf_lite is a lightweight Python adaptation of the [vcf2maf Perl tool](https://github.com/mskcc/vcf2maf), designed to convert the VCF format to MAF format without adding variant annotations. For annotating a MAF file, [Genome Nexus](https://github.com/genome-nexus/genome-nexus-annotation-pipeline) can be utilized.

### Installation using pip

```
pip3 install vcf2maf_lite
```

Usage:

```
vcf2maf_lite --help

Usage: vcf2maf_lite [OPTIONS]

Options:
  -i, --input-data TEXT        A list of .vcf files or input data directories,
                               separated by commas  [required]
  -o, --output-directory TEXT  output data directory [optional]
  -c, --center TEXT            name of center (standard MAF field = 'Center')
                               [optional]
  -s, --sequence-source TEXT   Sequencing source (standard MAF field =
                               'Sequencing_Source'), e.g., WXS or WGS
                               [optional]
  -t, --tumor-id TEXT          The ID of the tumor sample utilized in the
                               genotype columns of the VCF file. [optional]
  -n, --normal-id TEXT         The ID of the normal sample utilized in the
                               genotype columns of the VCF file. [optional]
  -a, --retain-info TEXT       Comma-delimited names of INFO fields to retain
                               as extra columns in MAF [optional]
  -f, --retain-fmt TEXT        Comma-delimited names of FORMAT fields to
                               retain as extra columns in MAF [optional]
  --help                       Show this message and exit.
```

Example Usage:
```
vcf2maf_lite --input-data /data/vcf --output-directory /data/maf/ --center CTR --sequence-source WGS --tumor-id Tumor --normal-id Normal --retain-info Custom_filters,AC,AF,AC_nfe_seu,AC_afr,AF_afr --retain-fmt alt_count_raw,ref_count_raw,depth_raw
```

This command converts the VCF files in /vcf folder to MAF format. 
- The `--input-data` option is used to specify either a single VCF file or a directory containing multiple VCF files (separated by commas). This option supports passing multiple input files or directories at once.
- The `--output-directory` option allows you to specify the directory where the MAF files will be saved. If no output path is provided, the default output directory `vcf2maf_output` will be used in the current working directory. 
- The `--tumor-id` option allows you to specify the ID of the tumor sample used in the genotype columns of the VCF file. If the option is not used, the script will automatically identify the tumor ID from either the `tumor_sample` keyword in the meta data lines or the sample columns from VCF header.
- The `--normal-id` option allows you to specify the ID of the normal sample used in the genotype columns of the VCF file. If the option is not used, the script will automatically identify the normal ID from either the `normal_sample` keyword in the meta data lines or the sample columns from VCF header.
- The `--retain-info` option allows you to specify the INFO fields to be retained as additional columns in the MAF. If the option is not used, standard MAF columns are included by default.
- The `--retain-fmt` option allows you to specify the FORMAT fields to be retained as additional columns in the MAF. If the option is not used, standard MAF columns are included by default.


Importing to Python Scripts:
```
from vcf2maf_lite.vcf2maf_lite import main
from click.testing import CliRunner

runner = CliRunner()
runner.invoke(main, ['--input-data','test_vcf.vcf'])
```

