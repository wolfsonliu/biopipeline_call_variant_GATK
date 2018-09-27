#! /bin/bash
set -e
function usage {
    echo "Usage: $0 -r reference -f input.fastq -d outputdirectory -l outputfilelabel -t thread" 1>&2
}


while getopts "hr:f:d:l:t:" opt; do
    case $opt in
        h)
            usage
            ;;
        r)
            REFERENCEFILE=$OPTARG
            ;;
        f)
            # input fastq file
            INPUTFQ=$OPTARG
            ;;
        d)
            # output directory
            OUT_DIR=$OPTARG
            ;;
        l)
            # output name
            OUTLABEL=$OPTARG
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

echo "[$(date)] < Start >  FastQC: "${INPUTFQ}
$Fastqc -o ${OUT_DIR}/log/ -f fastq ${INPUTFQ}
echo "[$(date)] <  End  >  FastQC: result in "${OUT_DIR}/log

####################
# Mapping: bwa
#     input:  fastq
#     input:  bowtie2index
#     output: sorted.bam
####################

echo "[$(date)] < Start >  bwa: "${INPUTFQ}
$Bwa mem -M ${BWA_INDEX_DIR}/genome.fa ${INPUTFQ} 2> ${OUT_DIR}/log/bwa.log | \
    samtools view -bS - | \
    samtools sort -@ ${PARR} - -o ${OUT_DIR}/bam/${OUTLABEL}.bam
echo "[$(date)] <  End  >  bwa: "${OUT_DIR}/bam/${OUTLABEL}.bam

####################
# Read Group: picard
#    input:  bam
#    output: sorted.bam
####################

echo "[$(date)] < Start >  Picard AddOrReplaceReadGroups: "${OUT_DIR}/bam/${OUTLABEL}.bam
$Picard AddOrReplaceReadGroups \
     I=${OUT_DIR}/bam/${OUTLABEL}.bam \
     O=${OUT_DIR}/bam/${OUTLABEL}.sorted.bam \
     SORT_ORDER=coordinate \
     RGID=${OUTLABEL} \
     RGLB=bwa \
     RGPL=illumina \
     RGPU=unit1 \
     RGSM=20 2> ${OUT_DIR}/log/picard.addgroup.log
echo "[$(date)] <  End  >  Picard AddOrReplaceReadGroups: "${OUT_DIR}/bam/${OUTLABEL}.sorted.bam

####################
# Mark Duplicates: picard
#    input:  sorted.bam
#    output: markdup.sorted.bam
#    output: markdup.txt
####################

echo "[$(date)] < Start >  Picard MarkDuplicates: "${OUT_DIR}/bam/${OUTLABEL}.sorted.bam
$Picard MarkDuplicates I=${OUT_DIR}/bam/${OUTLABEL}.sorted.bam \
       O=${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam \
       M=${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.txt 2> ${OUT_DIR}/log/picard.markduplicates.log
echo "[$(date)] <  End  >  Picard MarkDuplicates: "${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam

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

echo "[$(date)] < Start >  GATK HaplotypeCaller: "${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam
$Gatk --java-options "-Xmx8G" HaplotypeCaller \
     -R ${GENOME_FA} \
     --output-mode EMIT_VARIANTS_ONLY \
     -ERC GVCF \
     -I ${OUT_DIR}/bam/${OUTLABEL}.sorted.markdup.bam \
     -O ${OUT_DIR}/vcf/${OUTLABEL}.g.vcf.gz 2> ${OUT_DIR}/log/gatk.HaplotypeCaller.log
echo "[$(date)] <  End  >  GATK HaplotypeCaller: "${OUT_DIR}/bam/${OUTLABEL}.g.vcf.gz

####################
# Convert gvcf to vcf: GATK GenotypeGVCFs
#    input:  gatk.haplotype.g.vcf
#    output: vcf
####################

echo "[$(date)] < Start >  GATK GenotypeGVCFs: "${OUT_DIR}/bam/${OUTLABEL}.g.vcf.gz
$Gatk --java-options "-Xmx8G" GenotypeGVCFs \
     -R ${GENOME_FA} \
     -V ${OUT_DIR}/vcf/${OUTLABEL}.g.vcf.gz \
     -O ${OUT_DIR}/vcf/${OUTLABEL}.vcf.gz 2> ${OUT_DIR}/log/gatk.GenotypeGVCFs.log
echo "[$(date)] <  End  >  GATK GenotypeGVCFs: "${OUT_DIR}/bam/${OUTLABEL}.vcf.gz

####################
# Filter: GATK VariantFiltration
#    input: vcf
#    ouput: vcf
####################

echo "[$(date)] < Start >  GATK VariantFiltration: "${OUT_DIR}/bam/${OUTLABEL}.vcf.gz
$Gatk VariantFiltration \
     -R ${GENOME_FA} \
     -V ${OUT_DIR}/vcf/${OUTLABEL}.vcf.gz \
     -O ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
     --filter-name "SNPfilter" 2> ${OUT_DIR}/log/gatk.VariantFiltration.snp.log
echo "[$(date)] <  End  >  GATK VariantFiltration: "${OUT_DIR}/bam/${OUTLABEL}.filter0.vcf.gz

zcat ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz | grep -v "SNPfilter" > ${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf

rm -f ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz
rm -f ${OUT_DIR}/vcf/${OUTLABEL}.filter0.vcf.gz.tbi

####################
# Variants Annotation: annovar
#    input:  vcf.gz(.csi)
#    input:  clinvar.vcf.gz(.csi)
#    output: output dir vcfs
####################


echo "[$(date)] < Start >  ANNOVAR table_annovar.pl: "${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf
$Table_annovar --buildver hg38 \
               --remove --protocol refGene,clinvar_20180603,exac03,avsnp147,dbnsfp30a \
               --operation g,f,f,f,f -nastring . --polish --vcfinput \
               --out ${OUT_DIR}/vcf/${OUTLABEL} \
               ${OUT_DIR}/vcf/${OUTLABEL}.filter.vcf ${ANNOVAR_DIR}
echo "[$(date)] <  End  >  ANNOVAR table_annovar.pl: " ${OUT_DIR}/vcf/${OUTLABEL}.hg38_multianno.txt


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


$Makereport -l ${OUTLABEL} -p ${pipelinedir}/gatk4.pdf \
     -q ${OUT_DIR}/log/$(basename ${INPUTFQ%[.]*})_fastqc.zip \
     -s ${OUT_DIR}/log/${OUTLABEL}.sorted.markdup.bam.stats \
     -v ${OUT_DIR}/log/${OUTLABEL}.filter.vcf.stats \
     -f ${INPUTFQ} -b ${OUT_DIR}/${OUTLABEL}.anno.vcf \
     -d ${OUT_DIR}/report

cp ${OUT_DIR}/report/report.pdf ${OUT_DIR}/${OUTLABEL}_report_$(date -I"date").pdf

echo "All finished, results in: " ${OUT_DIR}
#####################
