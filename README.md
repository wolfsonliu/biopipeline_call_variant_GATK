# callvariant_GATK4 Pipeline #

* Take one fastq file as input.
* Generate VCF file as output.

## Requirement ##

### Software and packages ###

* FastQC (>=0.11): used in quality control of sequencing data.
* cutadapt (>=1.16): used in quality control of sequencing data.
* BWA (>=0.7): mapping FASTQ files to reference genome.
* samtools (==0.19): processing SAM/BAM files.
* bcftools (==0.19): processing VCF files.
* GATK (==4.0): variant calling.
* Picard (==2.0): processing SAM/BAM/VCF files.
* ANNOVAR (>=2018Apr16): annotate the VCF files.
* Python3 (>=3.6): processing data.
* numpy (>=1.15): python3 library for data computing.
* pandas (>=0.22): python3 library for data analysis.
* matplotlib (>=3.0): python3 library for graphics.

### Reference Data ###

* Reference genome: hg38 from iGenome website.
* Other data will be download automaticly by ANNOVAR software

## Source File Introduction ##

* callvariant_GATK4.sh: main pipeline script for variant calling.

```{shell}
Usage: ./callvariant_GATK4.sh -l outputfilelabel -r reference -p input1.fastq -q input2.fastq -d outputdirectory -t thread

Parameters:
    h: print help
    l: output label
    r: reference setting file path
    p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing
    q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing
    d: output directory path
    t: max numbers of threads to be used
```

* callvariant.ref: configuration file for pipeline script.

This file should contain following variables:
  * PIPELINE_DIR: the directory of this pipeline
  * PIPELINE_FIGURE: the pdf file used in report for the pipeline analysis schema.
    The pdf file should be in the pipeline directory.
  * BWA_INDEX_DIR: the directory for reference genome index used by BWA software.
  * GENOME_FA: reference genome FASTA file.
  * ANNOVAR_DIR: directory for data used by ANNOVAR software.
  * FastQC: the command or path for fastqc software.
  * Cutadapt: the command or path for cutadapt software.
  * Bwa: the command or path for BWA software.
  * Samtools: the command or path for samtools software.
  * Bcftools: the command or path for bcftools software.
  * Picard: the command or path for Picard software.
  * Gatk: the command or path for GATK 4 software.
  * Annotate_variantion: the command or path for annotation_variantion.pl script in ANNOVAR.
  * Table_annovar: the command or path for table_annovar.pl script in ANNOVAR.
  * Makereport: the command or path for the makereport scripts in this pipeline.
    The default path is in the PIPELINE_DIR.

* makereport: python3 scripts to generate reports, used by main pipeline script.

```{shell}
usage: makereport.zip [-h] [--label [LABEL]] --report-title [REPORT_TITLE]
                      --report-author [REPORT_AUTHOR] --fig-pipeline
                      [FIG_PIPELINE] --qczip1 [QCZIP1] [--qczip2 [QCZIP2]]
                      --samstat [SAMSTAT] --bcfstat [BCFSTAT] --vcf [VCF]
                      [--output-directory [OUTPUT_DIRECTORY]]
                      [--ref-genome [REF_GENOME]]

Generate Report TeX file.

optional arguments:
-h, --help            show this help message and exit
--label [LABEL]
--report-title [REPORT_TITLE]
--report-author [REPORT_AUTHOR]
--fig-pipeline [FIG_PIPELINE]
--qczip1 [QCZIP1]     FastQC zip file path, input the only file path for
                      single end sequencing or the first fastq FastQC zip
                      file path for paired end sequencing
--qczip2 [QCZIP2]     FastQC zip file path, not required for single end
                      sequencing or the second fastq FastQC zip file path
                      for paired end sequencing
--samstat [SAMSTAT]   samtools stats output txt file path
--bcfstat [BCFSTAT]   bcftools stats output txt file path
--vcf [VCF]           vcf file annotated by ANNOVAR
--output-directory [OUTPUT_DIRECTORY]
                      output directory path
--ref-genome [REF_GENOME]
```

Scripts in the makereport:
  * `__init__.py`: used as a package.
  * `__main__.py`: used as a zipped file.
  * `plotfuncs.py`: functions for visualizing the data.
  * `repfuncs.py`: assistant functions.
  * `reporter.py`: main script to generate the report TeX file.
  * `texfuncs.py`: TeX wrapper functions.
