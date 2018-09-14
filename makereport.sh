#! /bin/bash
function usage {
    echo "Usage: $0 -l label -p pipeline.png -q qc.zip -s samtools_stats.txt -v bcftools_stats.txt -f input.fastq -d outdir" 1>&2
}


while getopts "hl:p:q:s:v:f:d:" opt; do
    case $opt in
        h)
            usage
            ;;
        l)
            label=$OPTARG
            ;;
        p)
            figpipeline=$OPTARG
            ;;
        q)
            qcfile=$OPTARG
            ;;
        s)
            sstatfile=$OPTARG
            ;;
        v)
            vstatfile=$OPTARG
            ;;
        f)
            # input fastq file
            inputfq=$OPTARG
            ;;
        d)
            # output directory
            outputdir=$OPTARG
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

if [ -z $label -o -z $qcfile -o -z $inputfq -o -z $sstatfile -o -z $vstatfile ]; then
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
if [ -z $outputdir ]; then
    outdir=$(pwd)
else
    outdir=${outputdir}
fi

if [ ! -e ${outputdir}/fig ]; then
    mkdir -p ${outputdir}/fig
fi

echo "[$(date) ] Unzip FastQC zip file"
unzip ${qcfile} -d ${outputdir}

qcdir=${outputdir}/$(basename ${qcfile%[.]zip})

####################
echo "[$(date) ] Copy FastQC figures to output directory"
# get figure
cp ${figpipeline} ${outputdir}/fig/$(basename ${figpipeline})
cp ${qcdir}/Images/sequence_length_distribution.png ${outputdir}/fig/qc_seq_length_distribution.png
cp ${qcdir}/Images/per_base_quality.png ${outputdir}/fig/qc_base_quality_boxplot.png
cp ${qcdir}/Images/per_sequence_quality.png ${outputdir}/fig/qc_seq_quality_distribution.png
cp ${qcdir}/Images/per_base_sequence_content.png ${outputdir}/fig/qc_base_content.png
cp ${qcdir}/Images/per_sequence_gc_content.png ${outputdir}/fig/qc_seqe_gc.png

####################
echo "[$(date) ] Fetch FastQC data to output directory"
# get data
getblock ">>Per base sequence quality" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_base_quality.txt
getblock ">>Per sequence quality scores" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_seq_quality_distribution.txt
getblock ">>Per base sequence content" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_base_content.txt
getblock ">>Per sequence GC content" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_seq_gc.txt
getblock ">>Per base N content" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_base_n.txt
getblock ">>Sequence Length Distribution" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_seq_length_distribution.txt
getblock ">>Sequence Duplication Levels" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_seq_duplication.txt
getblock ">>Overrepresented sequences" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_seq_overrepresented.txt
getblock ">>Adapter Content" ">>END_MODULE" ${qcdir}/fastqc_data.txt > ${outputdir}/qc_adapter.txt

####################
echo "[$(date) ] Fetch FastQC data to output directory"
# get justification
echo "base_quality:" $(grep ">>Per base sequence quality" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "seq_quality:" $(grep ">>Per sequence quality scores" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "base_content:" $(grep ">>Per base sequence content" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "seq_gc:" $(grep ">>Per sequence GC content" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "base_n:" $(grep ">>Per base N content" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "seq_length_distribution:" $(grep ">>Sequence Length Distribution" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "seq_duplication:" $(grep ">>Sequence Duplication Levels" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "seq_overrepresented:" $(grep ">>Overrepresented sequences" ${qcdir}/fastqc_data.txt | cut -f2) | tee -a ${outputdir}/stats.txt
echo "adapter:" $(grep ">>Adapter Content" ${qcdir}/fastqc_data.txt | cut -f2)

####################
echo "[$(date) ] Fetch SAM statistics data to output directory"
# get samstats
cat ${sstatfile} | grep ^SN | cut -f 2- >> ${outputdir}/sam_summary_number.txt
echo -e "gc\tcount" > ${outputdir}/sam_GC_first.txt
cat ${sstatfile} | grep ^GCF | cut -f 2- >> ${outputdir}/sam_GC_first.txt
echo -e "gc\tcount" > ${outputdir}/sam_GC_last.txt
cat ${sstatfile} | grep ^GCL | cut -f 2- >> ${outputdir}/sam_GC_last.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${outputdir}/sam_ACGT.txt
cat ${sstatfile} | grep ^GCC | cut -f 2- >> ${outputdir}/sam_ACGT.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${outputdir}/sam_ACGT_first.txt
cat ${sstatfile} | grep ^FBC | cut -f 2- >> ${outputdir}/sam_ACGT_first.txt
echo -e "cycle\tA\tC\tG\tT\tN\t0" > ${outputdir}/sam_ACGT_last.txt
cat ${sstatfile} | grep ^LBC | cut -f 2- >> ${outputdir}/sam_ACGT_last.txt
echo -e "read_length\tcount" > ${outputdir}/sam_read_length.txt
cat ${sstatfile} | grep ^RL | cut -f 2- >> ${outputdir}/sam_read_length.txt
echo -e "read_length\tcount" > ${outputdir}/sam_read_length_first.txt
cat ${sstatfile} | grep ^FRL | cut -f 2- >> ${outputdir}/sam_read_length_first.txt
echo -e "read_length\tcount" > ${outputdir}/sam_read_length_last.txt
cat ${sstatfile} | grep ^LRL | cut -f 2- >> ${outputdir}/sam_read_length_last.txt
echo -e "length\tinsertion\tdeletion" > ${outputdir}/sam_indel.txt
cat ${sstatfile} | grep ^ID | cut -f 2- >> ${outputdir}/sam_indel.txt
echo -e "cycle\tinsertion_fwd\tinsertion_rev\tdeletion_fwd\tdelection_rev" > ${outputdir}/sam_cycle_indel.txt
cat ${sstatfile} | grep ^IC | cut -f 2- >> ${outputdir}/sam_cycle_indel.txt
cat ${sstatfile} | grep ^COV | cut -f 2- >> ${outputdir}/sam_coverage_distribution.txt

####################
echo "[$(date) ] Fetch VCF statistics data to output directory"
# get vcfstas
echo -e "id\tkey\tvalue" > ${outputdir}/vcf_summary_number.txt
cat ${vstatfile} | grep ^SN | cut -f 2- >> ${outputdir}/vcf_summary_number.txt
echo -e "id\tts\ttv\tts/tv\tts (1st ALT)\ttv (1st ALT)\tts/tv (1st ALT)" > ${outputdir}/vcf_tstv.txt
cat ${vstatfile} | grep ^TSTV | cut -f 2- >> ${outputdir}/vcf_tstv.txt
echo -e "id\tQuality\tnumber of SNPs\tnumber of transitions (1st ALT)\tnumber of transversions (1st ALT)\tnumber of indels" > ${outputdir}/vcf_qual.txt
cat ${vstatfile} | grep ^QUAL | cut -f 2- >> ${outputdir}/vcf_qual.txt
echo -e "id\tlength (deletions negative)\tcount" > ${outputdir}/vcf_indel_distribution.txt
cat ${vstatfile} | grep ^IDD | cut -f 2- >> ${outputdir}/vcf_indel_distribution.txt
echo -e "id\ttype\tcount" > ${outputdir}/vcf_substitution_types.txt
cat ${vstatfile} | grep ^ST | cut -f 2- >> ${outputdir}/vcf_substitution_types.txt
echo -e "id\tbin\tnumber of genotypes\tfraction of genotypes (%)\tnumber of sites\tfraction of sites (%)" > ${outputdir}/vcf_depth_distribution.txt
cat ${vstatfile} | grep ^DP | cut -f 2- >> ${outputdir}/vcf_depth_distribution.txt

####################
echo "[$(date) ] Store statistics data"
# make variable
mkrp_total_seq=$(grep "raw total sequences:" ${outputdir}/sam_summary_number.txt | cut -f 2)
mkrp_seq_length_mean=$(fq_seq_length_mean ${inputfq})
mkrp_seq_length_stat=$(fq_seq_length_stat ${inputfq})
mkrp_seq_length_min=$(echo $mkrp_seq_length_stat | cut -d"," -f 1)
mkrp_seq_length_q1=$(echo $mkrp_seq_length_stat | cut -d"," -f 2)
mkrp_seq_length_median=$(echo $mkrp_seq_length_stat | cut -d"," -f 3)
mkrp_seq_length_q3=$(echo $mkrp_seq_length_stat | cut -d"," -f 4)
mkrp_seq_length_max=$(echo $mkrp_seq_length_stat | cut -d"," -f 5)
mkrp_seq_gc_mean=$(fq_seq_gc_mean ${inputfq})
mkrp_seq_gc_stat=$(fq_seq_gc_stat ${inputfq})
mkrp_seq_gc_min=$(echo $mkrp_seq_gc_stat | cut -d"," -f 1)
mkrp_seq_gc_q1=$(echo $mkrp_seq_gc_stat | cut -d"," -f 2)
mkrp_seq_gc_median=$(echo $mkrp_seq_gc_stat | cut -d"," -f 3)
mkrp_seq_gc_q3=$(echo $mkrp_seq_gc_stat | cut -d"," -f 4)
mkrp_seq_gc_max=$(echo $mkrp_seq_gc_stat | cut -d"," -f 5)
mkrp_mapped_seq=$(grep "reads mapped:" ${outputdir}/sam_summary_number.txt | cut -f 2)
mkrp_unmapped_seq=$(grep "reads unmapped:" ${outputdir}/sam_summary_number.txt | cut -f 2)
mkrp_mapping_rate=0$(echo $mkrp_mapped_seq / $mkrp_total_seq | bc -l)
mkrp_total_variant=$(grep "number of records:" ${outputdir}/vcf_summary_number.txt | cut -f 3)
mkrp_snp=$(grep "number of SNPs:" ${outputdir}/vcf_summary_number.txt | cut -f 3)
mkrp_indel=$(grep "number of indels:" ${outputdir}/vcf_summary_number.txt | cut -f 3)

# store statistics
echo "total_seq:" $mkrp_total_seq | tee -a ${outputdir}/stats.txt
echo "seq_length_mean:" $mkrp_seq_length_mean | tee -a ${outputdir}/stats.txt
echo "seq_length_stat:" $mkrp_seq_length_stat | tee -a ${outputdir}/stats.txt
echo "seq_length_min:" $mkrp_seq_length_min | tee -a ${outputdir}/stats.txt
echo "seq_length_q1:" $mkrp_seq_length_q1 | tee -a ${outputdir}/stats.txt
echo "seq_length_median:" $mkrp_seq_length_median | tee -a ${outputdir}/stats.txt
echo "seq_length_q3:" $mkrp_seq_length_q3 | tee -a ${outputdir}/stats.txt
echo "seq_length_max:" $mkrp_seq_length_max | tee -a ${outputdir}/stats.txt
echo "seq_gc_mean:" $mkrp_seq_gc_mean | tee -a ${outputdir}/stats.txt
echo "seq_gc_stat:" $mkrp_seq_gc_stat | tee -a ${outputdir}/stats.txt
echo "seq_gc_min:" $mkrp_seq_gc_min | tee -a ${outputdir}/stats.txt
echo "seq_gc_q1:" $mkrp_seq_gc_q1 | tee -a ${outputdir}/stats.txt
echo "seq_gc_median:" $mkrp_seq_gc_median | tee -a ${outputdir}/stats.txt
echo "seq_gc_q3:" $mkrp_seq_gc_q3 | tee -a ${outputdir}/stats.txt
echo "seq_gc_max:" $mkrp_seq_gc_max | tee -a ${outputdir}/stats.txt
echo "mapped_seq:" $mkrp_mapped_seq | tee -a ${outputdir}/stats.txt
echo "unmapped_seq:" $mkrp_unmapped_seq | tee -a ${outputdir}/stats.txt
echo "mapping_rate:" $mkrp_mapping_rate | tee -a ${outputdir}/stats.txt
echo "total_variant:" $mkrp_total_variant | tee -a ${outputdir}/stats.txt
echo "snp:" $mkrp_snp | tee -a ${outputdir}/stats.txt
echo "indel:" $mkrp_indel | tee -a ${outputdir}/stats.txt

####################
echo "[$(date) ] Generate report file in output directory"
# Make report.tex
python3 reporter.py --output ${outputdir}/report.tex \
       --report-title "${label} Variant Analysis Report" \
       --report-author "MS Health Care Team" \
       --sum-total-seq $mkrp_total_seq \
       --sum-mapping-rate $mkrp_mapping_rate \
       --sum-refgenome GRCh38 \
       --sum-variants-number $mkrp_total_variant \
       --sum-snp-number $mkrp_snp \
       --sum-indel-number $mkrp_indel \
       --qc-seq-mode "Single End" \
       --qc-total-seq $mkrp_total_seq \
       --qc-total-pair 0 \
       --qc-seq-length-mean $mkrp_seq_length_mean \
       --qc-seq-length-min $mkrp_seq_length_min \
       --qc-seq-length-median $mkrp_seq_length_median \
       --qc-seq-length-max $mkrp_seq_length_max \
       --qc-seq-length-q1 $mkrp_seq_length_q1 \
       --qc-seq-length-q3 $mkrp_seq_length_q3 \
       --qc-seq-gc-mean $mkrp_seq_gc_mean \
       --qc-seq-gc-min $mkrp_seq_gc_min \
       --qc-seq-gc-median $mkrp_seq_gc_median \
       --qc-seq-gc-max $mkrp_seq_gc_max \
       --qc-seq-gc-q1 $mkrp_seq_gc_q1 \
       --qc-seq-gc-q3 $mkrp_seq_gc_q3 \
       --map-maptool "BWA MEM" \
       --map-maptool-v "0.7.17" \
       --map-refgenome "GRCh38" \
       --map-refgenome-v "UCSC hg38" \
       --map-samtool "samtools" \
       --map-samtool-v "1.7" \
       --map-mkdup "Picard" \
       --map-mkdup-v "2.18.11" \
       --map-total-seq $mkrp_total_seq \
       --map-mapped-seq $mkrp_mapped_seq \
       --map-unmapped-seq $mkrp_unmapped_seq \
       --vc-total-variant $mkrp_total_variant \
       --vc-snv-number $mkrp_snp \
       --vc-indel-number $mkrp_indel \
       --figpath "${outputdir}/fig/" \
       --fig-pipeline $(basename $figpipeline) \
       --fig-qc-seq-length-distribution "qc_seq_length_distribution.png" \
       --fig-qc-base-quality-boxplot "qc_base_quality_boxplot.png" \
       --fig-qc-seq-quality-distribution "qc_seq_quality_distribution.png" \
       --fig-qc-base-content "qc_base_content.png" \
       --fig-qc-seq-gc "qc_seqe_gc.png" \
       --fig-vc-snvtype "" \
       --fig-vc-snvanno ""


# make pdf

xelatex -output-directory=${outputdir} ${outputdir}/report.tex
xelatex -output-directory=${outputdir} ${outputdir}/report.tex

####################
