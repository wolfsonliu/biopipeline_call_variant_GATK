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
            OUTLABEL=$OPTARG
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

if [ -L $0 ] ; then
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

echo $PATH
################################################################################


if [ ! -e ${OUT_DIR}/log ]; then
    mkdir -p ${OUT_DIR}/log
fi


if [ ! -e ${OUT_DIR}/bam ]; then
    mkdir -p ${OUT_DIR}/bam
fi


if [ ! -e ${OUT_DIR}/vcf ]; then
    mkdir -p ${OUT_DIR}/vcf
fi

####################
# Quality Control: fastqc
#     input:  fastq
#     output: html
####################

start_stage "FastQC: "${INPUTFQ1}
$Fastqc -o ${OUT_DIR}/log/ -f fastq ${INPUTFQ1}
end_stage "FastQC: result in "${OUT_DIR}/log

if [ -n ${INPUTFQ2} ]; then
    start_stage "FastQC: "${INPUTFQ2}
    $Fastqc -o ${OUT_DIR}/log/ -f fastq ${INPUTFQ2}
    end_stage "FastQC: result in "${OUT_DIR}/log
fi


####################
# Mapping: bwa
#     input:  fastq
#     input:  bowtie2index
#     output: sorted.bam
####################

if [ -n ${INPUTFQ1} -a -n ${INPUTFQ2} ]; then
    start_stage Start  "bwa: ${INPUTFQ1} ${INPUTFQ2}"
    $Bwa mem -M ${BWA_INDEX_DIR}/genome.fa ${INPUTFQ1} ${INPUTFQ2} 2> ${OUT_DIR}/log/bwa.log | \
        samtools view -bS - | \
        samtools sort -@ ${PARR} - -o ${OUT_DIR}/bam/${OUTLABEL}.bam
elif [ -n ${INPUTFQ1} -a -z ${INPUTFQ2} ]; then
    start_stage "bwa: ${INPUTFQ1}"
    $Bwa mem -M ${BWA_INDEX_DIR}/genome.fa ${INPUTFQ1} 2> ${OUT_DIR}/log/bwa.log | \
        samtools view -bS - | \
        samtools sort -@ ${PARR} - -o ${OUT_DIR}/bam/${OUTLABEL}.bam
fi
end_stage "bwa: "${OUT_DIR}/bam/${OUTLABEL}.bam

####################
# Read Group: picard
#    input:  bam
#    output: sorted.bam
####################

start_stage "Picard AddOrReplaceReadGroups: "${OUT_DIR}/bam/${OUT_LABEL}.bam
$Picard AddOrReplaceReadGroups \
     I=${OUT_DIR}/bam/${OUTLABEL}.bam \
     O=${OUT_DIR}/bam/${OUTLABEL}.sorted.bam \
     SORT_ORDER=coordinate \
     RGID=${OUTLABEL} \
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
####################

start_stage "Picard MarkDuplicates: "${OUT_DIR}/bam/${OUT_LABEL}.sorted.bam
$Picard MarkDuplicates I=${OUT_DIR}/bam/${OUTLABEL}.sorted.bam \
       O=${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam \
       M=${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.txt 2> ${OUT_DIR}/log/picard.markduplicates.log
end_stage "Picard MarkDuplicates: "${OUT_DIR}/bam/${OUT_LABEL}.sorted.markdup.bam

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
####################

$Samtools stats ${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam > ${OUT_DIR}/log/${OUTLABEL}.sorted.markdup.bam.stats

####################
# Build bai: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bai
####################

$Samtools index -b ${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam

####################
# Call variants: GATK HaplotypeCaller
#    input:  markdup.sorted.bam
#    input:  reference fa
#    output: vcf
####################

start_stage "GATK HaplotypeCaller: "${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam
$Gatk --java-options "-Xmx8G" HaplotypeCaller \
     -R ${GENOME_FA} \
     --output-mode EMIT_VARIANTS_ONLY \
     -ERC GVCF \
     -I ${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam \
     -O ${OUT_DIR}/vcf/${OUTLABEL}.g.vcf.gz 2> ${OUT_DIR}/log/gatk.HaplotypeCaller.log
end_stage "GATK HaplotypeCaller: "${OUT_DIR}/bam/${OUTLABEL}.g.vcf.gz

####################
# Convert gvcf to vcf: GATK GenotypeGVCFs
#    input:  gatk.haplotype.g.vcf
#    output: vcf
####################

start_stage "GATK GenotypeGVCFs: "${OUT_DIR}/bam/${OUTLABEL}.g.vcf.gz
$Gatk --java-options "-Xmx8G" GenotypeGVCFs \
     -R ${GENOME_FA} \
     -V ${OUT_DIR}/vcf/${OUTLABEL}.g.vcf.gz \
     -O ${OUT_DIR}/vcf/${OUTLABEL}.vcf.gz 2> ${OUT_DIR}/log/gatk.GenotypeGVCFs.log
end_stage "GATK GenotypeGVCFs: "${OUT_DIR}/bam/${OUTLABEL}.vcf.gz

####################
# Filter: GATK VariantFiltration
#    input: vcf
#    ouput: vcf
####################

start_stage "GATK VariantFiltration: "${OUT_DIR}/bam/${OUTLABEL}.vcf.gz
$Gatk VariantFiltration \
     -R ${GENOME_FA} \
     -V ${OUT_DIR}/vcf/${OUTLABEL}.vcf.gz \
     -O ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
     --filter-name "SNPfilter" 2> ${OUT_DIR}/log/gatk.VariantFiltration.snp.log
end_stage "GATK VariantFiltration: "${OUT_DIR}/bam/${OUTLABEL}.filter0.vcf.gz

zcat ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz | grep -v "SNPfilter" > ${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf

rm -f ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz
rm -f ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz.tbi

####################
# Variants Annotation: annovar
#    input:  vcf.gz(.csi)
#    input:  clinvar.vcf.gz(.csi)
#    output: output dir vcfs
####################


start_stage "ANNOVAR table_annovar.pl: "${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf
$Table_annovar --buildver hg38 \
               --remove --protocol refGene,clinvar_20180603,exac03,avsnp147,dbnsfp30a \
               --operation g,f,f,f,f -nastring . --polish --vcfinput \
               --out ${OUT_DIR}/vcf/${OUTLABEL} \
               ${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf ${ANNOVAR_DIR}
end_stage "ANNOVAR table_annovar.pl: " ${OUT_DIR}/vcf/${OUTLABEL}.hg38_multianno.txt


cp ${OUT_DIR}/vcf/${OUTLABEL}.hg38_multianno.vcf ${OUT_DIR}/${OUTLABEL}.anno.vcf
cp ${OUT_DIR}/vcf/${OUTLABEL}.hg38_multianno.txt ${OUT_DIR}/${OUTLABEL}.anno.txt

echo "[$(date)] Annotated VCF: " ${OUT_DIR}/${OUTLABEL}.anno.vcf

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
####################

$Bcftools stats ${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf > ${OUT_DIR}/log/${OUTLABEL}.filter.vcf.stats


####################
# Generate report
####################


# $Makereport -l ${OUTLABEL} -r ${pipelinedir}/gatk4.pdf \
#      -f ${OUT_DIR}/log/$(basename ${INPUTFQ%[.]*})_fastqc.zip \
#      -s ${OUT_DIR}/log/${OUTLABEL}.sorted.markdup.bam.stats \
#      -v ${OUT_DIR}/log/${OUTLABEL}.filter.vcf.stats \
#      -p ${INPUTFQ1} -q ${INPUTFQ2} -b ${OUT_DIR}/${OUTLABEL}.anno.vcf \
#      -d ${OUT_DIR}/report

# cp ${OUT_DIR}/report/report.pdf ${OUT_DIR}/${OUTLABEL}_report_$(date -I"date").pdf

echo "All finished, results in: " ${OUT_DIR}
#####################
