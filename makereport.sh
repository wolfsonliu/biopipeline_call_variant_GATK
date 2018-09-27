#! /bin/bash
set -eu

function usage {
    echo "Usage: $0 -l label -p pipeline.png -q qc.zip -s samtools_stats.txt -v bcftools_stats.txt -f input.fastq -b anno.vcf -d outdir" 1>&2
}


while getopts "hl:p:q:s:v:f:b:d:" opt; do
    case ${opt} in
        h)
            usage
            ;;
        l)
            OPT_LABEL=${OPTARG}
            ;;
        p)
            OPT_FIGPIPELINE=${OPTARG}
            ;;
        q)
            OPT_QCFILE=${OPTARG}
            ;;
        s)
            OPT_SSTATFILE=${OPTARG}
            ;;
        v)
            OPT_VSTATFILE=${OPTARG}
            ;;
        f)
            # input fastq file
            OPT_INPUTFQ=${OPTARG}
            ;;
        b)
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

if [ -z ${OPT_LABEL} -o -z ${OPT_QCFILE} -o -z ${OPT_INPUTFQ} -o -z ${OPT_INPUTVCF} -o -z ${OPT_SSTATFILE} -o -z ${OPT_VSTATFILE} ]; then
   usage
   exit 1
fi

####################

function fq_seq_length_mean {
    awk 'BEGIN {
            lengthsum = 0;
            records = 0;
            }
        FNR %4 == 2 {
            records += 1;
            lengthsum += length($0);
            }
        END {
            print lengthsum/records;
            }' $@
}


function fq_seq_length_stat {
    awk 'BEGIN {}
        FNR %4 == 2 {
            readlength[FNR] = length($0);
            }
        END {
            asort(readlength);
            records = length(readlength);
            q1 = readlength[int(records/4)];
            q3 = readlength[int(records/4*3)];
            median = readlength[int(records/2 + 1)];
            if (records % 2 == 0) {
                median = (readlength[int(records/2)] + readlength[int(records/2 + 1)]) / 2;
            }
            print readlength[1] "," q1 "," median "," q3 "," readlength[records];
       }' $@
}

function fq_seq_gc_mean {
    awk 'BEGIN {
        gc = 0;
        records = 0;
        }
        FNR %4 == 2 {
            records += 1;
            seq = $0;
            gc += gsub(/[GCgc]/, "", seq) / length($0);
        }
        END {
            print gc/records;
        }' $@
}

function fq_seq_gc_stat {
    awk 'BEGIN {}
        FNR %4 == 2 {
v            seq = $0;
            gc[FNR] = gsub(/[GCgc]/, "", seq) / length($0);
        }
        END {
            asort(gc);
            records = length(gc);
            q1 = gc[int(records/4)];
            q3 = gc[int(records/4*3)];
            median = gc[int(records/2 + 1)];
            if (records % 2 == 0) {
                median = (gc[int(records/2)] + gc[int(records/2 + 1)]) / 2;
            }
            print gc[1] "," q1 "," median "," q3 "," gc[records];
        }' $@
}

function getblock {
    awk -v sflag="$1" -v eflag="$2" 'BEGIN {inblock = 0;}
        $0 ~ eflag {inblock = 0;}
        inblock > 0 {print $0;}
        $0 ~ sflag {inblock = 1;}' $3
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


echo "[$(date) ] Unzip FastQC zip file"
unzip ${OPT_QCFILE} -d ${OUT_DIR}

QCDIR=${OUT_DIR}/$(basename ${OPT_QCFILE%[.]zip})

####################

echo "[$(date) ] Fetch FastQC data to output directory"
# get data
getblock ">>Per base sequence quality" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_base_quality.txt
getblock ">>Per sequence quality scores" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_seq_quality_distribution.txt
getblock ">>Per base sequence content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_base_content.txt
getblock ">>Per sequence GC content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_seq_gc.txt
getblock ">>Per base N content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_base_n.txt
getblock ">>Sequence Length Distribution" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_seq_length_distribution.txt
getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_seq_duplication.txt
getblock ">>Overrepresented sequences" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_seq_overrepresented.txt
getblock ">>Adapter Content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUT_DIR}/qc_adapter.txt

####################

echo "[$(date) ] Fetch FastQC data to output directory"
# get justification
echo "base_quality:" $(grep ">>Per base sequence quality" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "seq_quality:" $(grep ">>Per sequence quality scores" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "base_content:" $(grep ">>Per base sequence content" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc:" $(grep ">>Per sequence GC content" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "base_n:" $(grep ">>Per base N content" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_distribution:" $(grep ">>Sequence Length Distribution" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "seq_duplication:" $(grep ">>Sequence Duplication Levels" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "seq_overrepresented:" $(grep ">>Overrepresented sequences" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo "adapter:" $(grep ">>Adapter Content" ${QCDIR}/fastqc_data.txt | cut -f2)

####################

echo "[$(date) ] Fetch SAM statistics data to output directory"
# get samstats
cat ${OPT_SSTATFILE} | grep ^SN | cut -f 2- >> ${OUT_DIR}/sam_summary_number.txt
echo -e "gc\tcount" > ${OUT_DIR}/sam_GC_first.txt
cat ${OPT_SSTATFILE} | grep ^GCF | cut -f 2- >> ${OUT_DIR}/sam_GC_first.txt
echo -e "gc\tcount" > ${OUT_DIR}/sam_GC_last.txt
cat ${OPT_SSTATFILE} | grep ^GCL | cut -f 2- >> ${OUT_DIR}/sam_GC_last.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUT_DIR}/sam_ACGT.txt
cat ${OPT_SSTATFILE} | grep ^GCC | cut -f 2- >> ${OUT_DIR}/sam_ACGT.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUT_DIR}/sam_ACGT_first.txt
cat ${OPT_SSTATFILE} | grep ^FBC | cut -f 2- >> ${OUT_DIR}/sam_ACGT_first.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUT_DIR}/sam_ACGT_last.txt
cat ${OPT_SSTATFILE} | grep ^LBC | cut -f 2- >> ${OUT_DIR}/sam_ACGT_last.txt
echo -e "read_length\tcount" > ${OUT_DIR}/sam_read_length.txt
cat ${OPT_SSTATFILE} | grep ^RL | cut -f 2- >> ${OUT_DIR}/sam_read_length.txt
echo -e "read_length\tcount" > ${OUT_DIR}/sam_read_length_first.txt
cat ${OPT_SSTATFILE} | grep ^FRL | cut -f 2- >> ${OUT_DIR}/sam_read_length_first.txt
echo -e "read_length\tcount" > ${OUT_DIR}/sam_read_length_last.txt
cat ${OPT_SSTATFILE} | grep ^LRL | cut -f 2- >> ${OUT_DIR}/sam_read_length_last.txt
echo -e "length\tinsertion\tdeletion" > ${OUT_DIR}/sam_indel.txt
cat ${OPT_SSTATFILE} | grep ^ID | cut -f 2- >> ${OUT_DIR}/sam_indel.txt
echo -e "cycle\tinsertion_fwd\tinsertion_rev\tdeletion_fwd\tdelection_rev" > ${OUT_DIR}/sam_cycle_indel.txt
cat ${OPT_SSTATFILE} | grep ^IC | cut -f 2- >> ${OUT_DIR}/sam_cycle_indel.txt
cat ${OPT_SSTATFILE} | grep ^COV | cut -f 2- >> ${OUT_DIR}/sam_coverage_distribution.txt

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

####################

echo "[$(date) ] Copy FastQC figures to output directory"
# get figure
cp ${OPT_FIGPIPELINE} ${OUT_DIR}/fig/$(basename ${OPT_FIGPIPELINE})

plotqc_seq_length.py --input ${OUT_DIR}/qc_seq_length_distribution.txt \
                     --output ${OUT_DIR}/fig/qc_seq_length_distribution.pdf
plotqc_base_quality.py --input ${OUT_DIR}/qc_base_quality.txt \
                       --output ${OUT_DIR}/fig/qc_base_quality_boxplot.pdf
plotqc_seq_quality.py --input ${OUT_DIR}/qc_seq_quality_distribution.txt \
                      --output ${OUT_DIR}/fig/qc_seq_quality_distribution.pdf
plotqc_base_content.py --input ${OUT_DIR}/qc_base_content.txt \
                       --output ${OUT_DIR}/fig/qc_base_content.pdf
plotqc_seq_gc.py --input ${OUT_DIR}/qc_seq_gc.txt \
   --output ${OUT_DIR}/fig/qc_seq_gc.pdf
plotvcf_piesnp.py --input ${OPT_INPUTVCF} --output ${OUT_DIR}/fig/vc_snvtype.pdf
plotvcf_annovcf.py --input ${OPT_INPUTVCF} --output ${OUT_DIR}/fig/vc_annovcf.pdf

####################

echo "[$(date) ] Store statistics data"
# make variablen
MKRP_TOTAL_SEQ=$(grep "raw total sequences:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
MKRP_SEQ_LENGTH_MEAN=$(fq_seq_length_mean ${OPT_INPUTFQ})
MKRP_SEQ_LENGTH_STAT=$(fq_seq_length_stat ${OPT_INPUTFQ})
MKRP_SEQ_LENGTH_MIN=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 1)
MKRP_SEQ_LENGTH_Q1=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 2)
MKRP_SEQ_LENGTH_MEDIAN=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 3)
MKRP_SEQ_LENGTH_Q3=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 4)
MKRP_SEQ_LENGTH_MAX=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 5)
MKRP_SEQ_GC_MEAN=$(fq_seq_gc_mean ${OPT_INPUTFQ})
MKRP_SEQ_GC_STAT=$(fq_seq_gc_stat ${OPT_INPUTFQ})
MKRP_SEQ_GC_MIN=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 1)
MKRP_SEQ_GC_Q1=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 2)
MKRP_SEQ_GC_MEDIAN=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 3)
MKRP_SEQ_GC_Q3=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 4)
MKRP_SEQ_GC_MAX=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 5)
MKRP_MAPPED_SEQ=$(grep "reads mapped:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
MKRP_UNMAPPED_SEQ=$(grep "reads unmapped:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
MKRP_MAPPING_RATE=$(echo 0 | awk -v num=${MKRP_MAPPED_SEQ} -v deno=$MKRP_TOTAL_SEQ '{print num/deno}')
MKRP_TOTAL_VARIANT=$(grep "number of records:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)
MKRP_SNP=$(grep "number of SNPs:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)
MKRP_INDEL=$(grep "number of indels:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)

# store statistics
echo "total_seq:" ${MKRP_TOTAL_SEQ} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_mean:" ${MKRP_SEQ_LENGTH_MEAN} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_stat:" ${MKRP_SEQ_LENGTH_STAT} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_min:" ${MKRP_SEQ_LENGTH_MIN} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_q1:" ${MKRP_SEQ_LENGTH_Q1} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_median:" ${MKRP_SEQ_LENGTH_MEDIAN} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_q3:" ${MKRP_SEQ_LENGTH_Q3} | tee -a ${OUT_DIR}/stats.txt
echo "seq_length_max:" ${MKRP_SEQ_LENGTH_MAX} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_mean:" ${MKRP_SEQ_GC_MEAN} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_stat:" ${MKRP_SEQ_GC_STAT} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_min:" ${MKRP_SEQ_GC_MIN} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_q1:" ${MKRP_SEQ_GC_Q1} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_median:" ${MKRP_SEQ_GC_MEDIAN} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_q3:" ${MKRP_SEQ_GC_Q3} | tee -a ${OUT_DIR}/stats.txt
echo "seq_gc_max:" ${MKRP_SEQ_GC_MAX} | tee -a ${OUT_DIR}/stats.txt
echo "mapped_seq:" ${MKRP_MAPPED_SEQ} | tee -a ${OUT_DIR}/stats.txt
echo "unmapped_seq:" ${MKRP_UNMAPPED_SEQ} | tee -a ${OUT_DIR}/stats.txt
echo "mapping_rate:" ${MKRP_MAPPING_RATE} | tee -a ${OUT_DIR}/stats.txt
echo "total_variant:" ${MKRP_TOTAL_VARIANT} | tee -a ${OUT_DIR}/stats.txt
echo "snp:" ${MKRP_SNP} | tee -a ${OUT_DIR}/stats.txt
echo "indel:" ${MKRP_INDEL} | tee -a ${OUT_DIR}/stats.txt

####################

echo "[$(date) ] Generate report file in output directory"
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
       --qc-seq-mode "Single End" \
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

xelatex -interaction=nonstopmode -output-directory=${OUT_DIR} ${OUT_DIR}/report.tex
xelatex -interaction=nonstopmode -output-directory=${OUT_DIR} ${OUT_DIR}/report.tex

####################
