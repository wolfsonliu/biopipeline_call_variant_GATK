#! /bin/bash
set -eu

function usage {
    echo "Usage: $0 -l label -p pipeline.png -q qc.zip -s samtools_stats.txt -v bcftools_stats.txt -f input.fastq -d outdir" 1>&2
}


while getopts "hl:p:q:s:v:f:d:" opt; do
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

if [ -z ${OPT_LABEL} -o -z ${OPT_QCFILE} -o -z ${OPT_INPUTFQ} -o -z ${OPT_SSTATFILE} -o -z ${OPT_VSTATFILE} ]; then
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
            seq = $0;
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
    OUTDIR=$(pwd)
else
    OUTDIR=${OPT_OUTDIR}
fi

if [ ! -e ${OUTDIR}/fig ]; then
    mkdir -p ${OUTDIR}/fig
fi


echo "[$(date) ] Unzip FastQC zip file"
unzip ${OPT_QCFILE} -d ${OUTDIR}

QCDIR=${OUTDIR}/$(basename ${OPT_QCFILE%[.]zip})

####################
echo "[$(date) ] Copy FastQC figures to output directory"
# get figure
cp ${OPT_FIGPIPELINE} ${OUTDIR}/fig/$(basename ${OPT_FIGPIPELINE})
cp ${QCDIR}/Images/sequence_length_distribution.png ${OUTDIR}/fig/qc_seq_length_distribution.png
cp ${QCDIR}/Images/per_base_quality.png ${OUTDIR}/fig/qc_base_quality_boxplot.png
cp ${QCDIR}/Images/per_sequence_quality.png ${OUTDIR}/fig/qc_seq_quality_distribution.png
cp ${QCDIR}/Images/per_base_sequence_content.png ${OUTDIR}/fig/qc_base_content.png
cp ${QCDIR}/Images/per_sequence_gc_content.png ${OUTDIR}/fig/qc_seqe_gc.png

####################
echo "[$(date) ] Fetch FastQC data to output directory"
# get data
getblock ">>Per base sequence quality" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_base_quality.txt
getblock ">>Per sequence quality scores" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_seq_quality_distribution.txt
getblock ">>Per base sequence content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_base_content.txt
getblock ">>Per sequence GC content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_seq_gc.txt
getblock ">>Per base N content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_base_n.txt
getblock ">>Sequence Length Distribution" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_seq_length_distribution.txt
getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_seq_duplication.txt
getblock ">>Overrepresented sequences" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_seq_overrepresented.txt
getblock ">>Adapter Content" ">>END_MODULE" ${QCDIR}/fastqc_data.txt > ${OUTDIR}/qc_adapter.txt

####################
echo "[$(date) ] Fetch FastQC data to output directory"
# get justification
echo "base_quality:" $(grep ">>Per base sequence quality" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "seq_quality:" $(grep ">>Per sequence quality scores" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "base_content:" $(grep ">>Per base sequence content" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "seq_gc:" $(grep ">>Per sequence GC content" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "base_n:" $(grep ">>Per base N content" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "seq_length_distribution:" $(grep ">>Sequence Length Distribution" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "seq_duplication:" $(grep ">>Sequence Duplication Levels" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "seq_overrepresented:" $(grep ">>Overrepresented sequences" ${QCDIR}/fastqc_data.txt | cut -f2) | tee -a ${OUTDIR}/stats.txt
echo "adapter:" $(grep ">>Adapter Content" ${QCDIR}/fastqc_data.txt | cut -f2)

####################
echo "[$(date) ] Fetch SAM statistics data to output directory"
# get samstats
cat ${OPT_SSTATFILE} | grep ^SN | cut -f 2- >> ${OUTDIR}/sam_summary_number.txt
echo -e "gc\tcount" > ${OUTDIR}/sam_GC_first.txt
cat ${OPT_SSTATFILE} | grep ^GCF | cut -f 2- >> ${OUTDIR}/sam_GC_first.txt
echo -e "gc\tcount" > ${OUTDIR}/sam_GC_last.txt
cat ${OPT_SSTATFILE} | grep ^GCL | cut -f 2- >> ${OUTDIR}/sam_GC_last.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUTDIR}/sam_ACGT.txt
cat ${OPT_SSTATFILE} | grep ^GCC | cut -f 2- >> ${OUTDIR}/sam_ACGT.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUTDIR}/sam_ACGT_first.txt
cat ${OPT_SSTATFILE} | grep ^FBC | cut -f 2- >> ${OUTDIR}/sam_ACGT_first.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUTDIR}/sam_ACGT_last.txt
cat ${OPT_SSTATFILE} | grep ^LBC | cut -f 2- >> ${OUTDIR}/sam_ACGT_last.txt
echo -e "read_length\tcount" > ${OUTDIR}/sam_read_length.txt
cat ${OPT_SSTATFILE} | grep ^RL | cut -f 2- >> ${OUTDIR}/sam_read_length.txt
echo -e "read_length\tcount" > ${OUTDIR}/sam_read_length_first.txt
cat ${OPT_SSTATFILE} | grep ^FRL | cut -f 2- >> ${OUTDIR}/sam_read_length_first.txt
echo -e "read_length\tcount" > ${OUTDIR}/sam_read_length_last.txt
cat ${OPT_SSTATFILE} | grep ^LRL | cut -f 2- >> ${OUTDIR}/sam_read_length_last.txt
echo -e "length\tinsertion\tdeletion" > ${OUTDIR}/sam_indel.txt
cat ${OPT_SSTATFILE} | grep ^ID | cut -f 2- >> ${OUTDIR}/sam_indel.txt
echo -e "cycle\tinsertion_fwd\tinsertion_rev\tdeletion_fwd\tdelection_rev" > ${OUTDIR}/sam_cycle_indel.txt
cat ${OPT_SSTATFILE} | grep ^IC | cut -f 2- >> ${OUTDIR}/sam_cycle_indel.txt
cat ${OPT_SSTATFILE} | grep ^COV | cut -f 2- >> ${OUTDIR}/sam_coverage_distribution.txt

####################
echo "[$(date) ] Fetch VCF statistics data to output directory"
# get vcfstas
echo -e "id\tkey\tvalue" > ${OUTDIR}/vcf_summary_number.txt
cat ${OPT_VSTATFILE} | grep ^SN | cut -f 2- >> ${OUTDIR}/vcf_summary_number.txt
echo -e "id\tts\ttv\tts/tv\tts (1st ALT)\ttv (1st ALT)\tts/tv (1st ALT)" > ${OUTDIR}/vcf_tstv.txt
cat ${OPT_VSTATFILE} | grep ^TSTV | cut -f 2- >> ${OUTDIR}/vcf_tstv.txt
echo -e "id\tQuality\tnumber of SNPs\tnumber of transitions (1st ALT)\tnumber of transversions (1st ALT)\tnumber of indels" > ${OUTDIR}/vcf_qual.txt
cat ${OPT_VSTATFILE} | grep ^QUAL | cut -f 2- >> ${OUTDIR}/vcf_qual.txt
echo -e "id\tlength (deletions negative)\tcount" > ${OUTDIR}/vcf_indel_distribution.txt
cat ${OPT_VSTATFILE} | grep ^IDD | cut -f 2- >> ${OUTDIR}/vcf_indel_distribution.txt
echo -e "id\ttype\tcount" > ${OUTDIR}/vcf_substitution_types.txt
cat ${OPT_VSTATFILE} | grep ^ST | cut -f 2- >> ${OUTDIR}/vcf_substitution_types.txt
echo -e "id\tbin\tnumber of genotypes\tfraction of genotypes (%)\tnumber of sites\tfraction of sites (%)" > ${OUTDIR}/vcf_depth_distribution.txt
cat ${OPT_VSTATFILE} | grep ^DP | cut -f 2- >> ${OUTDIR}/vcf_depth_distribution.txt

####################
echo "[$(date) ] Store statistics data"
# make variablen
MKRP_TOTAL_SEQ=$(grep "raw total sequences:" ${OUTDIR}/sam_summary_number.txt | cut -f 2)
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
MKRP_MAPPED_SEQ=$(grep "reads mapped:" ${OUTDIR}/sam_summary_number.txt | cut -f 2)
MKRP_UNMAPPED_SEQ=$(grep "reads unmapped:" ${OUTDIR}/sam_summary_number.txt | cut -f 2)
MKRP_MAPPING_RATE=0$(echo ${MKRP_MAPPED_SEQ} / $MKRP_TOTAL_SEQ | bc -l)
MKRP_TOTAL_VARIANT=$(grep "number of records:" ${OUTDIR}/vcf_summary_number.txt | cut -f 3)
MKRP_SNP=$(grep "number of SNPs:" ${OUTDIR}/vcf_summary_number.txt | cut -f 3)
MKRP_INDEL=$(grep "number of indels:" ${OUTDIR}/vcf_summary_number.txt | cut -f 3)

# store statistics
echo "total_seq:" ${MKRP_TOTAL_SEQ} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_mean:" ${MKRP_SEQ_LENGTH_MEAN} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_stat:" ${MKRP_SEQ_LENGTH_STAT} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_min:" ${MKRP_SEQ_LENGTH_MIN} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_q1:" ${MKRP_SEQ_LENGTH_Q1} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_median:" ${MKRP_SEQ_LENGTH_MEDIAN} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_q3:" ${MKRP_SEQ_LENGTH_Q3} | tee -a ${OUTDIR}/stats.txt
echo "seq_length_max:" ${MKRP_SEQ_LENGTH_MAX} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_mean:" ${MKRP_SEQ_GC_MEAN} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_stat:" ${MKRP_SEQ_GC_STAT} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_min:" ${MKRP_SEQ_GC_MIN} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_q1:" ${MKRP_SEQ_GC_Q1} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_median:" ${MKRP_SEQ_GC_MEDIAN} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_q3:" ${MKRP_SEQ_GC_Q3} | tee -a ${OUTDIR}/stats.txt
echo "seq_gc_max:" ${MKRP_SEQ_GC_MAX} | tee -a ${OUTDIR}/stats.txt
echo "mapped_seq:" ${MKRP_MAPPED_SEQ} | tee -a ${OUTDIR}/stats.txt
echo "unmapped_seq:" ${MKRP_UNMAPPED_SEQ} | tee -a ${OUTDIR}/stats.txt
echo "mapping_rate:" ${MKRP_MAPPING_RATE} | tee -a ${OUTDIR}/stats.txt
echo "total_variant:" ${MKRP_TOTAL_VARIANT} | tee -a ${OUTDIR}/stats.txt
echo "snp:" ${MKRP_SNP} | tee -a ${OUTDIR}/stats.txt
echo "indel:" ${MKRP_INDEL} | tee -a ${OUTDIR}/stats.txt


####################
echo "[$(date) ] Generate report file in output directory"
# Make report.tex
reporter.py --output ${OUTDIR}/report.tex \
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
       --figpath "${OUTDIR}/fig/" \
       --fig-pipeline $(basename ${OPT_FIGPIPELINE}) \
       --fig-qc-seq-length-distribution "qc_seq_length_distribution.png" \
       --fig-qc-base-quality-boxplot "qc_base_quality_boxplot.png" \
       --fig-qc-seq-quality-distribution "qc_seq_quality_distribution.png" \
       --fig-qc-base-content "qc_base_content.png" \
       --fig-qc-seq-gc "qc_seqe_gc.png" \
       --fig-vc-snvtype "" \
       --fig-vc-snvanno ""


# make pdf

xelatex -output-directory=${OUTDIR} ${OUTDIR}/report.tex
xelatex -output-directory=${OUTDIR} ${OUTDIR}/report.tex

####################
