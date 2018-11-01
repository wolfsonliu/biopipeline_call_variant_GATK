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
set -e

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

while getopts "hl:r:p:q:f:g:s:t:v:b:d:" opt; do
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
SEQ_MODE="Single-End"
if [ -n "${OPT_INPUTFQ2}" ]; then
    SEQ_MODE="Paired-End"
fi


info_stage "Unzip FastQC zip file"
unzip ${OPT_QCFILE1} -d ${OUT_DIR}
QCDIR1=${OUT_DIR}/$(basename ${OPT_QCFILE1%[.]zip})
echo ${SEQ_MODE}
if [ "${SEQ_MODE}" = "Paired-End" ]; then
    info_stage "Unzip Paired FastQC zip file"
    echo unzip ${OPT_QCFILE2} -d ${OUT_DIR}
    QCDIR2=${OUT_DIR}/$(basename ${OPT_QCFILE2%[.]zip})
fi

#  copy pipeline figure
cp ${OPT_FIGPIPELINE} ${OUT_DIR}/fig/$(basename ${OPT_FIGPIPELINE})

####################

info_stage "Fetch FastQC data to output directory"
# get data
getblock ">>Per base sequence quality" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_base_quality.txt
getblock ">>Per sequence quality scores" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_seq_quality_distribution.txt
getblock ">>Per base sequence content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_base_content.txt
getblock ">>Per sequence GC content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_seq_gc.txt
getblock ">>Per base N content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_base_n.txt
getblock ">>Sequence Length Distribution" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_seq_length_distribution.txt
getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_seq_duplication.txt
getblock ">>Overrepresented sequences" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_seq_overrepresented.txt
getblock ">>Adapter Content" ">>END_MODULE" ${QCDIR1}/fastqc_data.txt > ${OUT_DIR}/qc_adapter.txt

info_stage "Generate QC figures to output directory"
# generate figures
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


if [ ${SEQ_MODE} = "Paired End" ]; then
    info_stage "Fetch Paired FastQC data to output directory"
    getblock ">>Per base sequence quality" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_base_quality.txt
    getblock ">>Per sequence quality scores" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_seq_quality_distribution.txt
    getblock ">>Per base sequence content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_base_content.txt
    getblock ">>Per sequence GC content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_seq_gc.txt
    getblock ">>Per base N content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_base_n.txt
    getblock ">>Sequence Length Distribution" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_seq_length_distribution.txt
    getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_seq_duplication.txt
    getblock ">>Overrepresented sequences" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_seq_overrepresented.txt
    getblock ">>Adapter Content" ">>END_MODULE" ${QCDIR2}/fastqc_data.txt > ${OUT_DIR}/qc2_adapter.txt
    info_stage "Generate Paired QC figures to output directory"
    # generate figures
    plotqc_seq_length.py --input ${OUT_DIR}/qc2_seq_length_distribution.txt \
                         --output ${OUT_DIR}/fig/qc2_seq_length_distribution.pdf
    plotqc_base_quality.py --input ${OUT_DIR}/qc2_base_quality.txt \
                           --output ${OUT_DIR}/fig/qc2_base_quality_boxplot.pdf
    plotqc_seq_quality.py --input ${OUT_DIR}/qc2_seq_quality_distribution.txt \
                          --output ${OUT_DIR}/fig/qc2_seq_quality_distribution.pdf
    plotqc_base_content.py --input ${OUT_DIR}/qc2_base_content.txt \
                           --output ${OUT_DIR}/fig/qc2_base_content.pdf
    plotqc_seq_gc.py --input ${OUT_DIR}/qc2_seq_gc.txt \
                     --output ${OUT_DIR}/fig/qc2_seq_gc.pdf
fi


####################

info_stage "Fetch FastQC data to output directory"
# get justification
echo $(basename ${OPT_INPUTFQ1}) "base_quality:" $(grep ">>Per base sequence quality" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_quality:" $(grep ">>Per sequence quality scores" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo $(basename ${OPT_INPUTFQ1}) "base_content:" $(grep ">>Per base sequence content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_gc:" $(grep ">>Per sequence GC content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
#echo $(basename ${OPT_INPUTFQ1}) "base_n:" $(grep ">>Per base N content" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
echo $(basename ${OPT_INPUTFQ1}) "seq_length_distribution:" $(grep ">>Sequence Length Distribution" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
#echo $(basename ${OPT_INPUTFQ1}) "seq_duplication:" $(grep ">>Sequence Duplication Levels" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
#echo $(basename ${OPT_INPUTFQ1}) "seq_overrepresented:" $(grep ">>Overrepresented sequences" ${QCDIR1}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
#echo $(basename ${OPT_INPUTFQ1}) "adapter:" $(grep ">>Adapter Content" ${QCDIR1}/fastqc_data.txt | cut -f2)

if [ ${SEQ_MODE} = "Paired End" ]; then
    echo $(basename ${OPT_INPUTFQ2}) "base_quality:" $(grep ">>Per base sequence quality" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_quality:" $(grep ">>Per sequence quality scores" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    echo $(basename ${OPT_INPUTFQ2}) "base_content:" $(grep ">>Per base sequence content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_gc:" $(grep ">>Per sequence GC content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    #echo $(basename ${OPT_INPUTFQ2}) "base_n:" $(grep ">>Per base N content" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    echo $(basename ${OPT_INPUTFQ2}) "seq_length_distribution:" $(grep ">>Sequence Length Distribution" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    #echo $(basename ${OPT_INPUTFQ2}) "seq_duplication:" $(grep ">>Sequence Duplication Levels" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    #echo $(basename ${OPT_INPUTFQ2}) "seq_overrepresented:" $(grep ">>Overrepresented sequences" ${QCDIR2}/fastqc_data.txt | cut -f2) | tee -a ${OUT_DIR}/stats.txt
    #echo $(basename ${OPT_INPUTFQ2}) "adapter:" $(grep ">>Adapter Content" ${QCDIR2}/fastqc_data.txt | cut -f2)
fi

####################

info_stage "Fetch SAM statistics data to output directory"
# get samstats
cat ${OPT_SSTATFILE} | grep ^SN | cut -f 2- > ${OUT_DIR}/sam_summary_number.txt
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

echo "[$(date)] Fetch VCF statistics data to output directory"
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

info_stage "Store statistics data"
# make variablen
MKRP_TOTAL_SEQ=$(grep "raw total sequences:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
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
MKRP_MAPPED_SEQ=$(grep "reads mapped:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
MKRP_UNMAPPED_SEQ=$(grep "reads unmapped:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
MKRP_MAPPING_RATE=$(echo 0 | awk -v num=${MKRP_MAPPED_SEQ} -v deno=$MKRP_TOTAL_SEQ '{print num/deno}')
MKRP_TOTAL_VARIANT=$(grep "number of records:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)
MKRP_SNP=$(grep "number of SNPs:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)
MKRP_INDEL=$(grep "number of indels:" ${OUT_DIR}/vcf_summary_number.txt | cut -f 3)

if [ SEQ_MODE = "Single End" ]; then
    MKRP_MAPPED_PAIR=0
else
    MKRP_MAPPED_PAIR=$(grep "reads mapped and paired:" ${OUT_DIR}/sam_summary_number.txt | cut -f 2)
fi
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

info_stage "Generate report file in output directory"
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
       --seq-mode "${SEQ_MODE}" \
       --total-seq ${MKRP_TOTAL_SEQ} \
       --total-pair 0 \
       --qc-warn $(grep warn ${OUT_DIR}/stats.txt | cut -d " " -f 2 | tr -d ":") \
       --qc-fail $(grep fail ${OUT_DIR}/stats.txt | cut -d " " -f 2 | tr -d ":") \
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
info_stage "Generate report pdf"
xelatex -8bit -interaction=nonstopmode -output-directory=${OUT_DIR} ${OUT_DIR}/report.tex
xelatex -8bit -interaction=nonstopmode -output-directory=${OUT_DIR} ${OUT_DIR}/report.tex

####################
