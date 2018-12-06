#! /bin/bash
####################
#report_prepare_vcf.sh
#     Description:
#         make variant calling analysis report
#     Usage:
#         report_prepare_vcf.sh -l label -r pipeline.png -p input1.fastq -q input2.fastq -f input1.qc.zip -g input2.qc.zip -s samtools_stats.txt -t bcftools_stats.txt -v anno.vcf -d outdir
#     Parameters:
#         h: print help
#         l: output label
#         r: analysis pipeline figure used in report
#         p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing
#         q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing
#         f: FastQC zip file path, input the only file path for single end sequencing or the first fastq FastQC zip file path for paired end sequencing
#         g: FastQC zip file path, not required for single end sequencing or the second fastq FastQC zip file path for paired end sequencing
#         s: samtools stats output txt file path
#         t: bcftools stats output txt file path
#         v: annotated vcf file path
#         d: output directory path
####################
set -eu

function usage {
    echo "Usage: $0 -l label -t bcftools_stats.txt -v anno.vcf -d outdir" 1>&2
}

function helpinfo {
    usage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "" 1>&2
    echo "    h: print help " 1>&2
    echo "    l: output label " 1>&2
    echo "    t: bcftools stats output txt file path " 1>&2
    echo "    v: annotated vcf file path " 1>&2
    echo "    d: output directory path " 1>&2

}

while getopts "hl:r:p:q:f:g:s:v:b:d:" opt; do
    case ${opt} in
        h)
            helpinfo
            exit 0
            ;;
        l)
            OPT_LABEL=${OPTARG}
            ;;
        t)
            OPT_VSTATFILE=${OPTARG}
            ;;
        v)
            # input fastq file
            OPT_INPUTVCF=${OPTARG}
            ;;
        d)
            # output directory
            OPT_OUTDIR=${OPTARG}
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

shift $((OPTIND-1))

if [ -z ${OPT_LABEL} -o -z ${OPT_INPUTVCF} -o -z ${OPT_VSTATFILE} ]; then
   usage
   exit 1
fi

####################

function info_stage {
    echo "[$(date)] " $@
}

####################

if [ -z ${OPT_OUTDIR} ]; then
    OUT_DIR=$(pwd)
else
    OUT_DIR=${OPT_OUTDIR}
fi

if [ ! -e ${OUT_DIR}/fig ]; then
    mkdir -p ${OUT_DIR}/fig
fi

####################

echo "[$(date) ] Fetch VCF statistics data to output directory"
# get vcfstas
echo -e "id\tkey\tvalue" > ${OUT_DIR}/vcf_summary_number.txt
cat ${OPT_VSTATFILE} | grep ^SN | cut -f 2- >> ${OUT_DIR}/vcf_summary_number.txt
echo -e "id\tts\ttv\tts/tv\tts (1st ALT)\ttv (1st ALT)\tts/tv (1st ALT)" > ${OUT_DIR}/vcf_tstv.txt
cat ${OPT_VSTATFILE} | grep ^TSTV | cut -f 2- >> ${OUT_DIR}/vcf_tstv.txt
echo -e "id\tQuality\tnumber of SNPs\tnumber of transitions (1st ALT)\tnumber of transversions (1st ALT)\tnumber of indels" > ${OUT_DIR}/vcf_qual.txt
cat ${OPT_VSTATFILE} | grep ^QUAL | cut -f 2- >> ${OUT_DIR}/vcf_qual.txt
echo -e "id\tlength (deletions negative)\tcount" > ${OUT_DIR}/vcf_indel_distribution.txt
cat ${OPT_VSTATFILE} | grep ^IDD | cut -f 2- >> ${OUT_DIR}/vcf_indel_distribution.txt
echo -e "id\ttype\tcount" > ${OUT_DIR}/vcf_substitution_types.txt
cat ${OPT_VSTATFILE} | grep ^ST | cut -f 2- >> ${OUT_DIR}/vcf_substitution_types.txt
echo -e "id\tbin\tnumber of genotypes\tfraction of genotypes (%)\tnumber of sites\tfraction of sites (%)" > ${OUT_DIR}/vcf_depth_distribution.txt
cat ${OPT_VSTATFILE} | grep ^DP | cut -f 2- >> ${OUT_DIR}/vcf_depth_distribution.txt

plotvcf_piesnp.py --input ${OPT_INPUTVCF} --output ${OUT_DIR}/fig/vc_snvtype.pdf
plotvcf_annovcf.py --input ${OPT_INPUTVCF} --output ${OUT_DIR}/fig/vc_annovcf.pdf

####################

MKRP_TOTAL_VARIANT=$(grep "number of records:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)
MKRP_SNP=$(grep "number of SNPs:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)
MKRP_INDEL=$(grep "number of indels:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)

echo "total_variant:" ${MKRP_TOTAL_VARIANT} | tee -a ${OUT_DIR}/stats_vcf.txt
echo "snp:" ${MKRP_SNP} | tee -a ${OUT_DIR}/stats_vcf.txt
echo "indel:" ${MKRP_INDEL} | tee -a ${OUT_DIR}/stats_vcf.txt

####################
