#! /bin/bash
####################
# report_prepare_fq.sh
#     Description:
#         make variant calling analysis report
#     Usage:
#         makereport.sh -l label -r pipeline.png -p input1.fastq -q input2.fastq -f input1.qc.zip -g input2.qc.zip -s samtools_stats.txt -t bcftools_stats.txt -v anno.vcf -d outdir
#     Parameters:
#         h: print help
#         l: output label
#         p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing
#         q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing
#         f: FastQC zip file path, input the only file path for single end sequencing or the first fastq FastQC zip file path for paired end sequencing
#         g: FastQC zip file path, not required for single end sequencing or the second fastq FastQC zip file path for paired end sequencing
#         d: output directory path
####################
set -eu

function usage {
    echo "Usage: $0 -l label -p input1.fastq -q input2.fastq -f input1.qc.zip -g input2.qc.zip -d outdir" 1>&2
}

function helpinfo {
    usage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "" 1>&2
    echo "    h: print help " 1>&2
    echo "    l: output label " 1>&2
    echo "    p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing " 1>&2
    echo "    q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing " 1>&2
    echo "    f: FastQC zip file path, input the only file path for single end sequencing or the first fastq FastQC zip file path for paired end sequencing " 1>&2
    echo "    g: FastQC zip file path, not required for single end sequencing or the second fastq FastQC zip file path for paired end sequencing " 1>&2
    echo "    d: output directory path " 1>&2

}

while getopts "hl:p:q:f:g:d:" opt; do
    case ${opt} in
        h)
            helpinfo
            exit 0
            ;;
        l)
            OPT_LABEL=${OPTARG}
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

if [ -z ${OPT_LABEL} -o -z ${OPT_QCFILE1} -o -z ${OPT_INPUTFQ1} ]; then
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

# test the mode: paired end or single end
if [ -n ${OPT_INPUTFQ1} -a -n ${OPT_INPUTFQ2} ]; then
    SEQ_MODE="Paired-End"
elif [ -n ${OPT_INPUTFQ1} -a -z ${OPT_INPUTFQ2} ]; then
    SEQ_MODE="Single-End"
fi


info_stage "Unzip FastQC zip file"
unzip ${OPT_QCFILE1} -d ${OUT_DIR}
QCDIR1=${OUT_DIR}/$(basename ${OPT_QCFILE1%[.]zip})
if [ -n ${OPT_QCFILE2} ]; then
    info_stage "Unzip Paired FastQC zip file"
    unzip ${OPT_QCFILE2} -d ${OUT_DIR}
    QCDIR2=${OUT_DIR}/$(basename ${OPT_QCFILE2%[.]zip})
fi


####################

info_stage "Fetch FastQC data to output directory"
# get data
getblock ">>Per base sequence quality" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_base_quality.txt
getblock ">>Per sequence quality scores" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_seq_quality_distribution.txt
getblock ">>Per base sequence content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_base_content.txt
getblock ">>Per sequence GC content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_seq_gc.txt
getblock ">>Per base N content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_base_n.txt
getblock ">>Sequence Length Distribution" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_seq_length_distribution.txt
getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_seq_duplication.txt
getblock ">>Overrepresented sequences" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_seq_overrepresented.txt
getblock ">>Adapter Content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc1_adapter.txt

info_stage "Generate QC figures to output directory"
# generate figures
plotqc_seq_length.py --input ${OUT_DIR}/${OPT_LABEL}.qc1_seq_length_distribution.txt \
                     --output ${OUT_DIR}/fig/${OPT_LABEL}.qc1_seq_length_distribution.pdf
plotqc_base_quality.py --input ${OUT_DIR}/${OPT_LABEL}.qc1_base_quality.txt \
                       --output ${OUT_DIR}/fig/${OPT_LABEL}.qc1_base_quality_boxplot.pdf
plotqc_seq_quality.py --input ${OUT_DIR}/${OPT_LABEL}.qc1_seq_quality_distribution.txt \
                      --output ${OUT_DIR}/fig/${OPT_LABEL}.qc1_seq_quality_distribution.pdf
plotqc_base_content.py --input ${OUT_DIR}/${OPT_LABEL}.qc1_base_content.txt \
                       --output ${OUT_DIR}/fig/${OPT_LABEL}.qc1_base_content.pdf
plotqc_seq_gc.py --input ${OUT_DIR}/${OPT_LABEL}.qc1_seq_gc.txt \
                 --output ${OUT_DIR}/fig/${OPT_LABEL}.qc1_seq_gc.pdf


if [ "${SEQ_MODE}" = "Paired-End" ]; then
    info_stage "Fetch Paired FastQC data to output directory"
    getblock ">>Per base sequence quality" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_base_quality.txt
    getblock ">>Per sequence quality scores" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_seq_quality_distribution.txt
    getblock ">>Per base sequence content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_base_content.txt
    getblock ">>Per sequence GC content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_seq_gc.txt
    getblock ">>Per base N content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_base_n.txt
    getblock ">>Sequence Length Distribution" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_seq_length_distribution.txt
    getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_seq_duplication.txt
    getblock ">>Overrepresented sequences" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_seq_overrepresented.txt
    getblock ">>Adapter Content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/${OPT_LABEL}.qc2_adapter.txt
    info_stage "Generate Paired QC figures to output directory"
    # generate figures
    plotqc_seq_length.py --input ${OUT_DIR}/${OPT_LABEL}.qc2_seq_length_distribution.txt \
                         --output ${OUT_DIR}/fig/${OPT_LABEL}.qc2_seq_length_distribution.pdf
    plotqc_base_quality.py --input ${OUT_DIR}/${OPT_LABEL}.qc2_base_quality.txt \
                           --output ${OUT_DIR}/fig/${OPT_LABEL}.qc2_base_quality_boxplot.pdf
    plotqc_seq_quality.py --input ${OUT_DIR}/${OPT_LABEL}.qc2_seq_quality_distribution.txt \
                          --output ${OUT_DIR}/fig/${OPT_LABEL}.qc2_seq_quality_distribution.pdf
    plotqc_base_content.py --input ${OUT_DIR}/${OPT_LABEL}.qc2_base_content.txt \
                           --output ${OUT_DIR}/fig/${OPT_LABEL}.qc2_base_content.pdf
    plotqc_seq_gc.py --input ${OUT_DIR}/${OPT_LABEL}.qc2_seq_gc.txt \
                     --output ${OUT_DIR}/fig/${OPT_LABEL}.qc2_seq_gc.pdf
fi


####################

info_stage "Fetch FastQC data to output directory"
# get justification
echo $(basename ${OPT_INPUTFQ1}) "base_quality1:" $(grep ">>Per base sequence quality" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_quality1:" $(grep ">>Per sequence quality scores" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "base_content1:" $(grep ">>Per base sequence content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_gc1:" $(grep ">>Per sequence GC content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "base_n1:" $(grep ">>Per base N content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_length_distribution1:" $(grep ">>Sequence Length Distribution" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_duplication1:" $(grep ">>Sequence Duplication Levels" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_overrepresented1:" $(grep ">>Overrepresented sequences" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo $(basename ${OPT_INPUTFQ1}) "adapter1:" $(grep ">>Adapter Content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt

if [ "${SEQ_MODE}" = "Paired-End" ]; then
    echo $(basename ${OPT_INPUTFQ2}) "base_quality2:" $(grep ">>Per base sequence quality" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_quality2:" $(grep ">>Per sequence quality scores" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "base_content2:" $(grep ">>Per base sequence content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_gc2:" $(grep ">>Per sequence GC content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "base_n2:" $(grep ">>Per base N content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_length_distribution2:" $(grep ">>Sequence Length Distribution" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_duplication2:" $(grep ">>Sequence Duplication Levels" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_overrepresented2:" $(grep ">>Overrepresented sequences" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
    echo $(basename ${OPT_INPUTFQ2}) "adapter2:" $(grep ">>Adapter Content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
fi

####################

info_stage "Store statistics data"
# make variablen
MKRP_SEQ_LENGTH_MEAN=$(fq_seq_length_mean ${OPT_INPUTFQ1})
MKRP_SEQ_LENGTH_STAT=$(fq_seq_length_stat ${OPT_INPUTFQ1})
MKRP_SEQ_LENGTH_MIN=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 1)
MKRP_SEQ_LENGTH_Q1=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 2)
MKRP_SEQ_LENGTH_MEDIAN=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 3)
MKRP_SEQ_LENGTH_Q3=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 4)
MKRP_SEQ_LENGTH_MAX=$(echo ${MKRP_SEQ_LENGTH_STAT} | cut -d"," -f 5)
MKRP_SEQ_GC_MEAN=$(fq_seq_gc_mean ${OPT_INPUTFQ1})
MKRP_SEQ_GC_STAT=$(fq_seq_gc_stat ${OPT_INPUTFQ1})
MKRP_SEQ_GC_MIN=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 1)
MKRP_SEQ_GC_Q1=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 2)
MKRP_SEQ_GC_MEDIAN=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 3)
MKRP_SEQ_GC_Q3=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 4)
MKRP_SEQ_GC_MAX=$(echo ${MKRP_SEQ_GC_STAT} | cut -d"," -f 5)

# store statistics
echo "seq_length_mean:"${MKRP_SEQ_LENGTH_MEAN} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_length_stat:"${MKRP_SEQ_LENGTH_STAT} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_length_min:"${MKRP_SEQ_LENGTH_MIN} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_length_q1:"${MKRP_SEQ_LENGTH_Q1} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_length_median:"${MKRP_SEQ_LENGTH_MEDIAN} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_length_q3:"${MKRP_SEQ_LENGTH_Q3} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_length_max:"${MKRP_SEQ_LENGTH_MAX} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_mean:"${MKRP_SEQ_GC_MEAN} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_stat:"${MKRP_SEQ_GC_STAT} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_min:"${MKRP_SEQ_GC_MIN} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_q1:"${MKRP_SEQ_GC_Q1} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_median:"${MKRP_SEQ_GC_MEDIAN} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_q3:"${MKRP_SEQ_GC_Q3} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt
echo "seq_gc_max:"${MKRP_SEQ_GC_MAX} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_fq.txt

####################
