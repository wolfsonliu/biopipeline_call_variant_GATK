#! /bin/bash
####################
# report_prepare_sam.sh
#     Description:
#         make variant calling analysis report
#     Usage:
#         makereport.sh -l label -r pipeline.png -p input1.fastq -q input2.fastq -f input1.qc.zip -g input2.qc.zip -s samtools_stats.txt -t bcftools_stats.txt -v anno.vcf -d outdir
#     Parameters:
#         h: print help
#         l: output label
#         s: samtools stats output txt file path
#         d: output directory path
####################
set -eu

function usage {
    echo "Usage: $0 -l label -s samtools_stats.txt -d outdir" 1>&2
}

function helpinfo {
    usage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "" 1>&2
    echo "    h: print help " 1>&2
    echo "    l: output label " 1>&2
    echo "    s: samtools stats output txt file path " 1>&2
    echo "    d: output directory path " 1>&2

}

while getopts "hl:s:d:" opt; do
    case ${opt} in
        h)
            helpinfo
            exit 0
            ;;
        l)
            OPT_LABEL=${OPTARG}
            ;;
        s)
            OPT_SSTATFILE=${OPTARG}
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

if [ -z ${OPT_LABEL} -o -z ${OPT_SSTATFILE} ]; then
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

####################

info_stage "Fetch SAM statistics data to output directory"
# get samstats
cat ${OPT_SSTATFILE} | grep ^SN | cut -f 2- > ${OUT_DIR}/${OPT_LABEL}.sam_summary_number.txt
echo -e "gc\tcount" > ${OUT_DIR}/${OPT_LABEL}.sam_GC_first.txt
cat ${OPT_SSTATFILE} | grep ^GCF | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_GC_first.txt
echo -e "gc\tcount" > ${OUT_DIR}/${OPT_LABEL}.sam_GC_last.txt
cat ${OPT_SSTATFILE} | grep ^GCL | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_GC_last.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUT_DIR}/${OPT_LABEL}.sam_ACGT.txt
cat ${OPT_SSTATFILE} | grep ^GCC | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_ACGT.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUT_DIR}/${OPT_LABEL}.sam_ACGT_first.txt
cat ${OPT_SSTATFILE} | grep ^FBC | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_ACGT_first.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${OUT_DIR}/${OPT_LABEL}.sam_ACGT_last.txt
cat ${OPT_SSTATFILE} | grep ^LBC | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_ACGT_last.txt
echo -e "read_length\tcount" > ${OUT_DIR}/${OPT_LABEL}.sam_read_length.txt
cat ${OPT_SSTATFILE} | grep ^RL | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_read_length.txt
echo -e "read_length\tcount" > ${OUT_DIR}/${OPT_LABEL}.sam_read_length_first.txt
cat ${OPT_SSTATFILE} | grep ^FRL | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_read_length_first.txt
echo -e "read_length\tcount" > ${OUT_DIR}/${OPT_LABEL}.sam_read_length_last.txt
cat ${OPT_SSTATFILE} | grep ^LRL | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_read_length_last.txt
echo -e "length\tinsertion\tdeletion" > ${OUT_DIR}/${OPT_LABEL}.sam_indel.txt
cat ${OPT_SSTATFILE} | grep ^ID | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_indel.txt
echo -e "cycle\tinsertion_fwd\tinsertion_rev\tdeletion_fwd\tdelection_rev" > ${OUT_DIR}/${OPT_LABEL}.sam_cycle_indel.txt
cat ${OPT_SSTATFILE} | grep ^IC | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_cycle_indel.txt
cat ${OPT_SSTATFILE} | grep ^COV | cut -f 2- >> ${OUT_DIR}/${OPT_LABEL}.sam_coverage_distribution.txt

####################

info_stage "Store statistics data"
# make variablen
MKRP_TOTAL_SEQ=$(grep "raw total sequences:" ${OUT_DIR}/${OPT_LABEL}.sam_summary_number.txt | cut -f 2)
MKRP_MAPPED_SEQ=$(grep "reads mapped:" ${OUT_DIR}/${OPT_LABEL}.sam_summary_number.txt | cut -f 2)
MKRP_UNMAPPED_SEQ=$(grep "reads unmapped:" ${OUT_DIR}/${OPT_LABEL}.sam_summary_number.txt | cut -f 2)
MKRP_MAPPING_RATE=$(echo 0 | awk -v num=${MKRP_MAPPED_SEQ} -v deno=$MKRP_TOTAL_SEQ '{print num/deno}')
MKRP_MAPPED_PAIR=$(grep "reads mapped and paired:" ${OUT_DIR}/${OPT_LABEL}.sam_summary_number.txt | cut -f 2)

# store statistics
echo "total_seq:"${MKRP_TOTAL_SEQ} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_sam.txt
echo "mapped_seq:"${MKRP_MAPPED_SEQ} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_sam.txt
echo "unmapped_seq:"${MKRP_UNMAPPED_SEQ} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_sam.txt
echo "mapping_rate:"${MKRP_MAPPING_RATE} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_sam.txt
echo "mapping_pair:"${MKRP_MAPPED_PAIR} | tee -a ${OUT_DIR}/${OPT_LABEL}.stats_sam.txt

####################
