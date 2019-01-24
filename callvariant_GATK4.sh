#! /bin/bash
####################
# callvariant_GATK4.sh
#     Description:
#         Run the GATK4 HaplotypeCaller pipeline to call variants of input fastq files.
#     Usage:
#         callvariant_GATK4.sh -l outputfilelabel -r reference -p input1.fastq -q input2.fastq -d outputdirectory -t thread
#     Parameters:
#         h: print help
#         l: output label
#         r: reference setting file path
#         p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing
#         q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing
#         d: output directory path
#         t: max numbers of threads to be used
####################

set -e

function usage {
    echo "Usage: $0 -l outputfilelabel -r reference -p input1.fastq -q input2.fastq -d outputdirectory -t thread" 1>&2
}

function helpinfo {
    usage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "" 1>&2
    echo "     h: print help" 1>&2
    echo "     l: output label" 1>&2
    echo "     r: reference setting file path" 1>&2
    echo "     p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing" 1>&2
    echo "     q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing" 1>&2
    echo "     d: output directory path" 1>&2
    echo "     t: max numbers of threads to be used" 1>&2
}

while getopts "hl:r:p:q:d:t:" opt; do
    case $opt in
        h)
            helpinfo
            exit 0
            ;;
        l)
            # output name
            OUT_LABEL=$OPTARG
            ;;
        r)
            REFERENCEFILE=$OPTARG
            ;;
        p)
            # input fastq file
            INPUTFQ1=$OPTARG
            ;;
        q)
            # input fastq file
            INPUTFQ2=$OPTARG
            ;;
        d)
            # output directory
            OUT_DIR=$OPTARG
            ;;
        t)
            # thread
            PARR=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if [ -L "$0" ] ; then
    DIR=$(dirname $(readlink -f $0))
else
    DIR=$(dirname $0)
fi

shift $((OPTIND-1))
function start_stage {
    echo "[$(date)] < Start >  " $@
}

function end_stage {
    echo "[$(date)] <  End  >  " $@
}


# read reference file path
source $REFERENCEFILE

if [ -n "${INPUTFQ1}" -a -n "${INPUTFQ2}" ]; then
    SEQ_MODE="Paired-End"
elif [ -n "${INPUTFQ1}" -a -z "${INPUTFQ2}" ]; then
    SEQ_MODE="Single-End"
fi
################################################################################


if [ ! -e "${OUT_DIR}/log" ]; then
    mkdir -p ${OUT_DIR}/log
fi

if [ ! -e "${OUT_DIR}/report" ]; then
    mkdir -p ${OUT_DIR}/report
fi

if [ ! -e "${OUT_DIR}/bam" ]; then
    mkdir -p ${OUT_DIR}/bam
fi


if [ ! -e "${OUT_DIR}/vcf" ]; then
    mkdir -p ${OUT_DIR}/vcf
fi

####################
# Quality Control: fastqc
#     input:  fastq
#     output: html
####################

function getblock {
    awk -v sflag="$1" -v eflag="$2" 'BEGIN {inblock = 0;}
        $0 ~ eflag {inblock = 0;}
        inblock > 0 {print $0;}
        $0 ~ sflag {inblock = 1;}' $3
}

function qc_judge {
    fail=0
    for x in `cut -f 1 $1`; do
        if [ ${x} == 'FAIL' ]; then
            fail=$((fail + 1))
        fi
    done
    echo ${fail}
}

function qc_grep_status {
    cat $1 | grep "$2" | cut -f 1
}



start_stage "FastQC: "${INPUTFQ1}
$Fastqc -o ${OUT_DIR}/log/ -f fastq ${INPUTFQ1}
end_stage "FastQC: result in "${OUT_DIR}/log

if [ "${SEQ_MODE}" = "Paired-End" ]; then
    start_stage "FastQC: "${INPUTFQ2}
    $Fastqc -o ${OUT_DIR}/log/ -f fastq ${INPUTFQ2}
    end_stage "FastQC: result in "${OUT_DIR}/log
fi

#

QC_NAME_BASE1=$(basename ${INPUTFQ1%.*})_fastqc
QC1_DIR1=${OUT_DIR}/log/${QC_NAME_BASE1}
unzip ${QC1_DIR1}.zip -d ${OUT_DIR}/log
failnum1=$(qc_judge ${QC1_DIR1}/summary.txt)
if (( failnum1 > 4 )); then
    echo "FAIL: Too many failures of $INPUTFQ1, QC file in:"${QC1_DIR1}.html
    exit 0
fi

if [ "${SEQ_MODE}" = "Paired-End" ]; then
    QC_NAME_BASE2=$(basename ${INPUTFQ2%.*})_fastqc
    QC1_DIR2=${OUT_DIR}/log/${QC_NAME_BASE2}
    unzip ${QC1_DIR2}.zip -d ${OUT_DIR}/log
    failnum2=$(qc_judge ${QC1_DIR2}/summary.txt)
    if (( failnum2 > 4 )); then
        echo "FAIL: Too many failures of $INPUTFQ2, QC file in:"${QC1_DIR2}.html
        exit 0
    fi
fi

TRIMED=0
if [ $(qc_grep_status ${QC1_DIR1}/summary.txt "Per base sequence quality") = 'FAIL' ]; then
    TRIMED=1
    getblock ">>Per base sequence quality" ">>END_MODULE" ${QC1_DIR1}/fastqc_data.txt > ${OUT_DIR}/log/${QC_NAME_BASE1}_base_quality.txt
    INFO_SEQ_LENGTH1=$(cat ${QC1_DIR1}/fastqc_data.txt | grep "Sequence length" | cut -f 2)
    INFO_SEQ_LENGTH1=${INFO_SEQ_LENGTH1##*-}
    MEAN_LENGTH=${INFO_SEQ_LENGTH1}
    if [ "${SEQ_MODE}" == "Paired-End" ]; then
        getblock ">>Per base sequence quality" ">>END_MODULE" ${QC1_DIR2}/fastqc_data.txt > ${OUT_DIR}/log/${QC_NAME_BASE2}_base_quality.txt
        INFO_SEQ_LENGTH2=$(cat ${QC1_DIR2}/fastqc_data.txt | grep "Sequence length" | cut -f 2)
        INFO_SEQ_LENGTH2=${INFO_SEQ_LENGTH2##*-}
        MEAN_LENGTH=$(((INFO_SEQ_LENGTH1 + INFO_SEQ_LENGTH2)/2))
    fi
    if [ "${SEQ_MODE}" = "Single-End" ]; then
        OUT_FQ1_NAME=$(basename $INPUTFQ1)
        OUT_FQ1_NAME=${OUT_FQ1_NAME%.*}.qc.fastq
        $Cutadapt -q 20,20 -m $((MEAN_LENGTH/2)) \
                  -o ${OUT_DIR}/bam/${OUT_FQ1_NAME} ${INPUTFQ1}
        OLD_INPUTFQ1=${INPUTFQ1}
        MAPFQ1=${OUT_DIR}/bam/${OUT_FQ1_NAME}
        $Fastqc -o ${OUT_DIR}/log/ -f fastq ${MAPFQ1}
        QC2_DIR1=${OUT_DIR}/log/$(basename ${MAPFQ1%.*})_fastqc
    elif [ "${SEQ_MODE}" = "Paired-End" ]; then
        OUT_FQ1_NAME=$(basename $INPUTFQ1)
        OUT_FQ1_NAME=${OUT_FQ1_NAME%.*}.qc.fastq
        OUT_FQ2_NAME=$(basename $INPUTFQ2)
        OUT_FQ2_NAME=${OUT_FQ2_NAME%.*}.qc.fastq
        $Cutadapt -q 20,20 -m $((MEAN_LENGTH/2)) \
                  -o ${OUT_DIR}/bam/${OUT_FQ1_NAME} -p ${OUT_DIR}/bam/${OUT_FQ2_NAME} \
                  ${INPUTFQ1} ${INPUTFQ2}
        OLD_INPUTFQ1=${INPUTFQ1}
        OLD_INPUTFQ2=${INPUTFQ2}
        MAPFQ1=${OUT_DIR}/bam/${OUT_FQ1_NAME}
        MAPFQ2=${OUT_DIR}/bam/${OUT_FQ2_NAME}
        $Fastqc -o ${OUT_DIR}/log/ -f fastq ${MAPFQ1}
        $Fastqc -o ${OUT_DIR}/log/ -f fastq ${MAPFQ2}
        QC2_DIR1=${OUT_DIR}/log/$(basename ${MAPFQ1%.*})_fastqc
        QC2_DIR2=${OUT_DIR}/log/$(basename ${MAPFQ2%.*})_fastqc
    fi
else
    if [ "${SEQ_MODE}" = "Single-End" ]; then
        OUT_FQ1_NAME=$(basename $INPUTFQ1)
        OUT_FQ1_NAME=${OUT_FQ1_NAME}
        ln -s ${INPUTFQ1} ${OUT_DIR}/bam/${OUT_FQ1_NAME}
        MAPFQ1=${OUT_DIR}/bam/${OUT_FQ1_NAME}
    else
        OUT_FQ1_NAME=$(basename $INPUTFQ1)
        OUT_FQ1_NAME=${OUT_FQ1_NAME}
        OUT_FQ2_NAME=$(basename $INPUTFQ2)
        OUT_FQ2_NAME=${OUT_FQ2_NAME}
        ln -s ${INPUTFQ1} ${OUT_DIR}/bam/${OUT_FQ1_NAME}
        ln -s ${INPUTFQ2} ${OUT_DIR}/bam/${OUT_FQ2_NAME}
        MAPFQ1=${OUT_DIR}/bam/${OUT_FQ1_NAME}
        MAPFQ2=${OUT_DIR}/bam/${OUT_FQ2_NAME}
    fi
    echo ""
fi


####################
# Mapping: bwa
#     input:  fastq
#     input:  bowtie2index
#     output: sorted.bam
# Description:
#     mapping fastq data to reference genome
####################

if [ "${SEQ_MODE}" = "Paired-End" ]; then
    start_stage Start  "bwa: ${INPUTFQ1} ${INPUTFQ2}"
    $Bwa mem -M ${BWA_INDEX_DIR}/genome.fa ${MAPFQ1} ${MAPFQ2} 2> ${OUT_DIR}/log/bwa.log | \
        samtools view -bS - | \
        samtools sort -@ ${PARR} - -o ${OUT_DIR}/bam/${OUT_LABEL}.bam
elif [ "${SEQ_MODE}" = "Single-End" ]; then
    start_stage "bwa: ${INPUTFQ1}"
    $Bwa mem -M ${BWA_INDEX_DIR}/genome.fa ${MAPFQ1} 2> ${OUT_DIR}/log/bwa.log | \
        samtools view -bS - | \
        samtools sort -@ ${PARR} - -o ${OUT_DIR}/bam/${OUT_LABEL}.bam
fi
end_stage "bwa: "${OUT_DIR}/bam/${OUT_LABEL}.bam

####################
# Read Group: picard
#    input:  bam
#    output: sorted.bam
# Description:
#    edit the reads information for GATK4 pipeline
####################

start_stage "Picard AddOrReplaceReadGroups: "${OUT_DIR}/bam/${OUT_LABEL}.bam
$Picard AddOrReplaceReadGroups \
     I=${OUT_DIR}/bam/${OUT_LABEL}.bam \
     O=${OUT_DIR}/bam/${OUT_LABEL}.sorted.bam \
     SORT_ORDER=coordinate \
     RGID=${OUT_LABEL} \
     RGLB=bwa \
     RGPL=illumina \
     RGPU=unit1 \
     RGSM=20 2> ${OUT_DIR}/log/picard.addgroup.log
end_stage "Picard AddOrReplaceReadGroups: "${OUT_DIR}/bam/${OUT_LABEL}.sorted.bam

####################
# Mark Duplicates: picard
#    input:  sorted.bam
#    output: markdup.sorted.bam
#    output: markdup.txt
# Description:
#    mark posible PCR related duplicates
####################

start_stage "Picard MarkDuplicates: "${OUT_DIR}/bam/${OUT_LABEL}.sorted.bam
$Picard MarkDuplicates I=${OUT_DIR}/bam/${OUT_LABEL}.sorted.bam \
       O=${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam \
       M=${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.txt 2> ${OUT_DIR}/log/picard.markduplicates.log
end_stage "Picard MarkDuplicates: "${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
# Decription:
#    statstics of the mapped BAM file, used for reporter
####################

$Samtools stats ${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam > ${OUT_DIR}/log/${OUT_LABEL}.sorted.markdup.bam.stats

####################
# Build bai: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bai
# Description:
#    index the BAM file
####################

$Samtools index -b ${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam

####################
# Call variants: GATK HaplotypeCaller
#    input:  markdup.sorted.bam
#    input:  reference fa
#    output: vcf
# Description:
#    call variants with GATK4 HaplotypeCaller
####################

start_stage "GATK HaplotypeCaller: "${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam
$Gatk --java-options "-Xmx8G" HaplotypeCaller \
     -R ${GENOME_FA} \
     --output-mode EMIT_VARIANTS_ONLY \
     -ERC GVCF \
     -I ${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam \
     -O ${OUT_DIR}/vcf/${OUT_LABEL}.g.vcf.gz 2> ${OUT_DIR}/log/gatk.HaplotypeCaller.log
end_stage "GATK HaplotypeCaller: "${OUT_DIR}/vcf/${OUT_LABEL}.g.vcf.gz

####################
# Convert gvcf to vcf: GATK GenotypeGVCFs
#    input:  gatk.haplotype.g.vcf
#    output: vcf
# Description:
#    convert the GVCF file to VCF file
####################

start_stage "GATK GenotypeGVCFs: "${OUT_DIR}/vcf/${OUT_LABEL}.g.vcf.gz
$Gatk --java-options "-Xmx8G" GenotypeGVCFs \
     -R ${GENOME_FA} \
     -V ${OUT_DIR}/vcf/${OUT_LABEL}.g.vcf.gz \
     -O ${OUT_DIR}/vcf/${OUT_LABEL}.vcf.gz 2> ${OUT_DIR}/log/gatk.GenotypeGVCFs.log
end_stage "GATK GenotypeGVCFs: "${OUT_DIR}/vcf/${OUT_LABEL}.vcf.gz

####################
# Filter: GATK VariantFiltration
#    input: vcf
#    ouput: vcf
# Description:
#    filter the VCF according to several scores
####################

start_stage "GATK VariantFiltration: "${OUT_DIR}/bam/${OUT_LABEL}.vcf.gz
$Gatk VariantFiltration \
     -R ${GENOME_FA} \
     -V ${OUT_DIR}/vcf/${OUT_LABEL}.vcf.gz \
     -O ${OUT_DIR}/vcf/${OUT_LABEL}.filter0.vcf.gz \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
     --filter-name "SNPfilter" 2> ${OUT_DIR}/log/gatk.VariantFiltration.snp.log
end_stage "GATK VariantFiltration: "${OUT_DIR}/bam/${OUT_LABEL}.filter0.vcf.gz

zcat ${OUT_DIR}/vcf/${OUT_LABEL}.filter0.vcf.gz | grep -v "SNPfilter" > ${OUT_DIR}/vcf/${OUT_LABEL}.filter.vcf

rm -f ${OUT_DIR}/vcf/${OUT_LABEL}.filter0.vcf.gz
rm -f ${OUT_DIR}/vcf/${OUT_LABEL}.filter0.vcf.gz.tbi

####################
# Variants Annotation: annovar
#    input:  vcf.gz(.csi)
#    input:  clinvar.vcf.gz(.csi)
#    output: output dir vcfs
# Description:
#    annotate the VCF file
####################


start_stage "ANNOVAR table_annovar.pl: "${OUT_DIR}/vcf/${OUT_LABEL}.filter.vcf
$Table_annovar --buildver hg38 \
               --remove --protocol refGene,clinvar_20180603,exac03,avsnp147,dbnsfp30a \
               --operation g,f,f,f,f -nastring . --polish --vcfinput \
               --out ${OUT_DIR}/vcf/${OUT_LABEL} \
               ${OUT_DIR}/vcf/${OUT_LABEL}.filter.vcf ${ANNOVAR_DIR}
end_stage "ANNOVAR table_annovar.pl: " ${OUT_DIR}/vcf/${OUT_LABEL}.hg38_multianno.txt

cp ${OUT_DIR}/vcf/${OUT_LABEL}.hg38_multianno.vcf ${OUT_DIR}/${OUT_LABEL}.anno.vcf
cp ${OUT_DIR}/vcf/${OUT_LABEL}.hg38_multianno.txt ${OUT_DIR}/${OUT_LABEL}.anno.txt

echo "[$(date)] Annotated VCF: " ${OUT_DIR}/${OUT_LABEL}.anno.vcf

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
# Decription:
#    statisctis fo the annotated VCF file
####################

$Bcftools stats ${OUT_DIR}/vcf/${OUT_LABEL}.filter.vcf > ${OUT_DIR}/log/${OUT_LABEL}.filter.vcf.stats


####################
# Generate report with the makereport scripts
####################
if [ "${SEQ_MODE}" = "Paired-End" ]; then
    $Makereport --label ${OUT_LABEL} \
                --report-title "${OUT_LABEL} Variant Analysis Report" \
                --report-author "MS Health Care" \
                --fig-pipeline ${PIPELINE_FIGURE}  \
                --qczip1 ${OUT_DIR}/log/$(basename ${MAPFQ1%.*})_fastqc.zip \
                --qczip2 ${OUT_DIR}/log/$(basename ${MAPFQ2%.*})_fastqc.zip \
                --samstat ${OUT_DIR}/log/${OUT_LABEL}.sorted.markdup.bam.stats \
                --bcfstat ${OUT_DIR}/log/${OUT_LABEL}.filter.vcf.stats \
                --vcf ${OUT_DIR}/${OUT_LABEL}.anno.vcf \
                --output-directory ${OUT_DIR}/report --ref-genome "GRCh38"
elif [ "${SEQ_MODE}" = "Single-End" ]; then
    $Makereport --label ${OUT_LABEL} \
                --report-title "${OUT_LABEL} Variant Analysis Report" \
                --report-author "MS Health Care" \
                --fig-pipeline ${PIPELINE_FIGURE}  \
                --qczip1 ${OUT_DIR}/log/$(basename ${MAPFQ1%.*})_fastqc.zip \
                --samstat ${OUT_DIR}/log/${OUT_LABEL}.sorted.markdup.bam.stats \
                --bcfstat ${OUT_DIR}/log/${OUT_LABEL}.filter.vcf.stats \
                --vcf ${OUT_DIR}/${OUT_LABEL}.anno.vcf \
                --output-directory ${OUT_DIR}/report --ref-genome "GRCh38"
fi

xelatex -8bit -interaction=nonstopmode -output-directory=${OUT_DIR}/report ${OUT_DIR}/report/report.tex
xelatex -8bit -interaction=nonstopmode -output-directory=${OUT_DIR}/report ${OUT_DIR}/report/report.tex

cp ${OUT_DIR}/report/report.pdf ${OUT_DIR}/${OUT_LABEL}_report.pdf

echo "All finished, results in: " ${OUT_DIR}
#####################
