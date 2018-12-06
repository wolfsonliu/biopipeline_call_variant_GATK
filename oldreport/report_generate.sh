#! /bin/bash
####################
# makereport.sh
#     Description:
#         make variant calling analysis report
#     Usage:
#         makereport.sh -l label -r pipeline.png -p input1.fastq -q input2.fastq -f input1.qc.zip -g input2.qc.zip -s samtools_stats.txt -t bcftools_stats.txt -v anno.vcf -d outdir
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
    echo "Usage: $0 -l label -r pipeline.png -p input1.fastq -q input2.fastq -f input1.qc.zip -g input2.qc.zip -s samtools_stats.txt -t bcftools_stats.txt -v anno.vcf -d outdir" 1>&2
}

function helpinfo {
    usage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "" 1>&2
    echo "    h: print help " 1>&2
    echo "    l: output label " 1>&2
    echo "    r: analysis pipeline figure used in report " 1>&2
    echo "    p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing " 1>&2
    echo "    q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing " 1>&2
    echo "    f: FastQC zip file path, input the only file path for single end sequencing or the first fastq FastQC zip file path for paired end sequencing " 1>&2
    echo "    g: FastQC zip file path, not required for single end sequencing or the second fastq FastQC zip file path for paired end sequencing " 1>&2
    echo "    s: samtools stats output txt file path " 1>&2
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
        r)
            OPT_FIGPIPELINE=${OPTARG}
            ;;
        p)
            # input fastq file
            OPT_INPUTFQ1=${OPTARG}
            ;;
        q)
            # input fastq file
            OPT_INPUTFQ2=${OPTARG}
            ;;
        f)
            OPT_QCFILE1=${OPTARG}
            ;;
        g)
            OPT_QCFILE2=${OPTARG}
            ;;
        s)
            OPT_SSTATFILE=${OPTARG}
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

if [ -z ${OPT_LABEL} -o -z ${OPT_QCFILE1} -o -z ${OPT_INPUTFQ1} -o -z ${OPT_INPUTVCF} -o -z ${OPT_SSTATFILE} -o -z ${OPT_VSTATFILE} ]; then
   usage
   exit 1
fi

####################
# Make report.tex
reporter.py --output ${OUT_DIR}/report.tex \
       --report-title "${OPT_LABEL} Variant Analysis Report" \
       --report-author "MS Health Care Team" \
       --sum-total-seq ${MKRP_TOTAL_SEQ} \
       --sum-mapping-rate ${MKRP_MAPPING_RATE} \
       --sum-refgenome GRCh38 \
       --sum-variants-number ${MKRP_TOTAL_VARIANT} \
       --sum-snp-number ${MKRP_SNP} \
       --sum-indel-number ${MKRP_INDEL} \
       --qc-seq-mode ${SEQ_MODE} \
       --qc-total-seq ${MKRP_TOTAL_SEQ} \
       --qc-total-pair 0 \
       --qc-seq-length-mean ${MKRP_SEQ_LENGTH_MEAN} \
       --qc-seq-length-min ${MKRP_SEQ_LENGTH_MIN} \
       --qc-seq-length-median ${MKRP_SEQ_LENGTH_MEDIAN} \
       --qc-seq-length-max ${MKRP_SEQ_LENGTH_MAX} \
       --qc-seq-length-q1 ${MKRP_SEQ_LENGTH_Q1} \
       --qc-seq-length-q3 ${MKRP_SEQ_LENGTH_Q3} \
       --qc-seq-gc-mean ${MKRP_SEQ_GC_MEAN} \
       --qc-seq-gc-min ${MKRP_SEQ_GC_MIN} \
       --qc-seq-gc-median ${MKRP_SEQ_GC_MEDIAN} \
       --qc-seq-gc-max ${MKRP_SEQ_GC_MAX} \
       --qc-seq-gc-q1 ${MKRP_SEQ_GC_Q1} \
       --qc-seq-gc-q3 ${MKRP_SEQ_GC_Q3} \
       --map-maptool "BWA MEM" \
       --map-maptool-v "0.7.17" \
       --map-refgenome "GRCh38" \
       --map-refgenome-v "UCSC hg38" \
       --map-samtool "samtools" \
       --map-samtool-v "1.7" \
       --map-mkdup "Picard" \
       --map-mkdup-v "2.18.11" \
       --map-total-seq ${MKRP_TOTAL_SEQ} \
       --map-mapped-seq ${MKRP_MAPPED_SEQ} \
       --map-unmapped-seq ${MKRP_UNMAPPED_SEQ} \
       --vc-total-variant ${MKRP_TOTAL_VARIANT} \
       --vc-snv-number ${MKRP_SNP} \
       --vc-indel-number ${MKRP_INDEL} \
       --figpath "${OUT_DIR}/fig/" \
       --fig-pipeline $(basename ${OPT_FIGPIPELINE}) \
       --fig-qc-seq-length-distribution "qc_seq_length_distribution.pdf" \
       --fig-qc-base-quality-boxplot "qc_base_quality_boxplot.pdf" \
       --fig-qc-seq-quality-distribution "qc_seq_quality_distribution.pdf" \
       --fig-qc-base-content "qc_base_content.pdf" \
       --fig-qc-seq-gc "qc_seq_gc.pdf" \
       --fig-vc-snvtype "vc_snvtype.pdf" \
       --fig-vc-annovcf "vc_annovcf.pdf"


# make pdf
xelatex -8bit -interaction=nonstopmode -output-directory=${OUT_DIR} ${OUT_DIR}/report.tex
xelatex -8bit -interaction=nonstopmode -output-directory=${OUT_DIR} ${OUT_DIR}/report.tex

####################
