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
            referencefile=$OPTARG
            ;;
        f)
            # input fastq file
            inputfq=$OPTARG
            ;;
        d)
            # output directory
            outputdir=$OPTARG
            ;;
        l)
            # output name
            outputlabel=$OPTARG
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

outdir=${outputdir}
outlabel=${outputlabel}
outprefix=${outdir}/${outlabel}

# read reference file path
source $referencefile

echo $PATH
################################################################################


if [ ! -e ${outdir}/log ]; then
    mkdir -p ${outdir}/log
fi


if [ ! -e ${outdir}/bam ]; then
    mkdir -p ${outdir}/bam
fi


if [ ! -e ${outdir}/vcf ]; then
    mkdir -p ${outdir}/vcf
fi

####################
# Quality Control: fastqc
#     input:  fastq
#     output: html
####################

echo "[$(date)] < Start >  FastQC: "${inputfq}
$Fastqc -o ${outdir}/log/ -f fastq ${inputfq}
echo "[$(date)] <  End  >  FastQC: result in "${outdir}/log

####################
# Mapping: bwa
#     input:  fastq
#     input:  bowtie2index
#     output: sorted.bam
####################

echo "[$(date)] < Start >  bwa: "${inputfq}
$Bwa mem -M ${bwaindexdir}/genome.fa ${inputfq} 2> ${outdir}/log/bwa.log | \
    samtools view -bS - | \
    samtools sort -@ ${PARR} - -o ${outdir}/bam/${outlabel}.bam
echo "[$(date)] <  End  >  bwa: "${outdir}/bam/${outlabel}.bam

####################
# Read Group: picard
#    input:  bam
#    output: sorted.bam
####################

echo "[$(date)] < Start >  Picard AddOrReplaceReadGroups: "${outdir}/bam/${outlabel}.bam
$Picard AddOrReplaceReadGroups \
     I=${outdir}/bam/${outlabel}.bam \
     O=${outdir}/bam/${outlabel}.sorted.bam \
     SORT_ORDER=coordinate \
     RGID=${outlabel} \
     RGLB=bwa \
     RGPL=illumina \
     RGPU=unit1 \
     RGSM=20 2> ${outdir}/log/picard.addgroup.log
echo "[$(date)] <  End  >  Picard AddOrReplaceReadGroups: "${outdir}/bam/${outlabel}.sorted.bam

####################
# Mark Duplicates: picard
#    input:  sorted.bam
#    output: markdup.sorted.bam
#    output: markdup.txt
####################

echo "[$(date)] < Start >  Picard MarkDuplicates: "${outdir}/bam/${outlabel}.sorted.bam
$Picard MarkDuplicates I=${outdir}/bam/${outlabel}.sorted.bam \
       O=${outdir}/bam/${outlabel}.sorted.markdup.bam \
       M=${outdir}/bam/${outlabel}.sorted.markdup.txt 2> ${outdir}/log/picard.markduplicates.log
echo "[$(date)] <  End  >  Picard MarkDuplicates: "${outdir}/bam/${outlabel}.sorted.markdup.bam

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
####################

$Samtools stats ${outdir}/bam/${outlabel}.sorted.markdup.bam > ${outdir}/log/${outlabel}.sorted.markdup.bam.stats

####################
# Build bai: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bai
####################

$Samtools index -b ${outdir}/bam/${outlabel}.sorted.markdup.bam

####################
# Call variants: GATK HaplotypeCaller
#    input:  markdup.sorted.bam
#    input:  reference fa
#    output: vcf
####################

echo "[$(date)] < Start >  GATK HaplotypeCaller: "${outdir}/bam/${outlabel}.sorted.markdup.bam
$Gatk --java-options "-Xmx8G" HaplotypeCaller \
     -R ${genomefa} \
     --output-mode EMIT_VARIANTS_ONLY \
     -ERC GVCF \
     -I ${outdir}/bam/${outlabel}.sorted.markdup.bam \
     -O ${outdir}/vcf/${outlabel}.g.vcf.gz 2> ${outdir}/log/gatk.HaplotypeCaller.log
echo "[$(date)] <  End  >  GATK HaplotypeCaller: "${outdir}/bam/${outlabel}.g.vcf.gz

####################
# Convert gvcf to vcf: GATK GenotypeGVCFs
#    input:  gatk.haplotype.g.vcf
#    output: vcf
####################

echo "[$(date)] < Start >  GATK GenotypeGVCFs: "${outdir}/bam/${outlabel}.g.vcf.gz
$Gatk --java-options "-Xmx8G" GenotypeGVCFs \
     -R ${genomefa} \
     -V ${outdir}/vcf/${outlabel}.g.vcf.gz \
     -O ${outdir}/vcf/${outlabel}.vcf.gz 2> ${outdir}/log/gatk.GenotypeGVCFs.log
echo "[$(date)] <  End  >  GATK GenotypeGVCFs: "${outdir}/bam/${outlabel}.vcf.gz

####################
# Filter: GATK VariantFiltration
#    input: vcf
#    ouput: vcf
####################

echo "[$(date)] < Start >  GATK VariantFiltration: "${outdir}/bam/${outlabel}.vcf.gz
$Gatk VariantFiltration \
     -R ${genomefa} \
     -V ${outdir}/vcf/${outlabel}.vcf.gz \
     -O ${outdir}/vcf/${outlabel}.filter0.vcf.gz \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
     --filter-name "SNPfilter" 2> ${outdir}/log/gatk.VariantFiltration.snp.log
echo "[$(date)] <  End  >  GATK VariantFiltration: "${outdir}/bam/${outlabel}.filter0.vcf.gz

zcat ${outdir}/vcf/${outlabel}.filter0.vcf.gz | grep -v "SNPfilter" > ${outdir}/vcf/${outlabel}.filter.vcf

rm -f ${outdir}/vcf/${outlabel}.filter0.vcf.gz
rm -f ${outdir}/vcf/${outlabel}.filter0.vcf.gz.tbi

####################
# Variants Annotation: annovar
#    input:  vcf.gz(.csi)
#    input:  clinvar.vcf.gz(.csi)
#    output: output dir vcfs
####################


echo "[$(date)] < Start >  ANNOVAR table_annovar.pl: "${outdir}/vcf/${outlabel}.filter.vcf
$Table_annovar --buildver hg38 \
               --remove --protocol refGene,clinvar_20180603,exac03,avsnp147,dbnsfp30a \
               --operation g,f,f,f,f -nastring . --polish --vcfinput \
               --out ${outdir}/vcf/${outlabel} \
               ${outdir}/vcf/${outlabel}.filter.vcf ${ANNOVARdir}
echo "[$(date)] <  End  >  ANNOVAR table_annovar.pl: " ${outdir}/vcf/${outlabel}.hg38_multianno.txt


cp ${outdir}/vcf/${outlabel}.hg38_multianno.vcf ${outdir}/${outlabel}.anno.vcf
cp ${outdir}/vcf/${outlabel}.hg38_multianno.txt ${outdir}/${outlabel}.anno.txt

echo "[$(date)] Annotated VCF: " ${outdir}/${outlabel}.anno.vcf

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
####################

$Bcftools stats ${outdir}/vcf/${outlabel}.filter.vcf > ${outdir}/log/${outlabel}.filter.vcf.stats


####################
# Generate report
####################


$Makereport -l ${outlabel} -p ${pipelinedir}/pipeline.png \
     -q ${outdir}/log/$(basename ${inputfq%[.]*})_fastqc.zip \
     -s ${outdir}/log/${outlabel}.sorted.markdup.bam.stats \
     -v ${outdir}/log/${outlabel}.filter.vcf.stats \
     -f ${inputfq} -b ${outdir}/${outlabel}.anno.vcf \
     -d ${outdir}/report

cp ${outdir}/report/report.pdf ${outdir}/${outlabel}_report_$(date -I"date").pdf

echo "All finished, results in: " ${outdir}
#####################
