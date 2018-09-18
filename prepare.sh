#! /bin/bash
set -e
# Used to download reference file or install softwares

basedir=${HOME}

softwaredir=${basedir}/Software

refdir=${basedir}/Reference

####################
# Reference
####################

# Reference Genome

wget -P ${refdir} "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz"

tar -xzvf ${refdir}/Homo_sapiens_UCSC_hg38.tar.gz

# java -jar picard.jar CreateSequenceDictionary REFERENCE=reference.fa OUTPUT=reference.dict


# ANNOVAR db

if [ ! -e ${refdir}/humandb ]; then
    mkdir -p ${refdir}/humandb
fi

annotate_variation.pl --buildver hg38 --downdb --webfrom annovar refGene ${refdir}/humandb/
annotate_variation.pl --buildver hg38 --downdb --webfrom annovar exac03 ${refdir}/humandb/
annotate_variation.pl --buildver hg38 --downdb --webfrom annovar avsnp147 ${refdir}/humandb/
annotate_variation.pl --buildver hg38 --downdb --webfrom annovar dbnsfp30a ${refdir}/humandb/
annotate_variation.pl --buildver hg38 --downdb --webfrom annovar clinvar_20180603 ${refdir}/humandb/

# GATK

wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz.tbi"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz"
wget -P ${refdir} "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz.tbi"

################################################################################
