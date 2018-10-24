#! /bin/bash
set -e

function usage {
    echo "Usage: $0 -p input1.fastq -q input2.fastq -z fastqc1.zip -s fastqc2.zip -d outputdirectory" 1>&2
}

function helpinfo {
    usage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "" 1>&2
    echo "     h: print help" 1>&2
    echo "     p: fastq file path, input the fastq file path for single end sequencing or the first fastq file path for paired end sequencing" 1>&2
    echo "     q: fastq file path, not required for single end sequencing or the second fastq file path for paired end sequencing" 1>&2
    echo "     z: zip file path for FastQC result of fastq file 1" 1>&2
    echo "     s: zip file path for FastQC result of fastq file 2" 1>&2
    echo "     d: output directory path" 1>&2
    echo "     t: max numbers of threads to be used" 1>&2
}

while getopts "hp:q:z:s:d:" opt; do
    case $opt in
        h)
            helpinfo
            exit 0
            ;;
        p)
            OPT_INPUTFQ1=$OPTARG
            ;;
        q)
            OPT_INPUTFQ2=$OPTARG
            ;;
        z)
            OPT_INPUTZIP1=$OPTARG
            ;;
        s)
            OPT_INPUTZIP2=$OPTARG
            ;;
        d)
            OPT_OUT_DIR=$OPTARG
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

####################

function getblock {
    awk -v sflag="$1" -v eflag="$2" 'BEGIN {inblock = 0;}
        $0 ~ eflag {inblock = 0;}
        inblock > 0 {print $0;}
        $0 ~ sflag {inblock = 1;}' $3
}

function qc_judge {
    fail=0
    for x in `cut -f 1 $1`; do
        if [ ${x} == 'FAIL' ]; then
            fail=$((fail + 1))
        fi
    done
    echo ${fail}
}

function qc_grep_status {
    cat $1 | grep "$2" | cut -f 1
}

####################

if [ -n ${OPT_INPUTFQ1} -a -n ${OPT_INPUTFQ2} ]; then
    SEQ_MODE="Paired End"
elif [ -n ${OPT_INPUTFQ1} -a -z ${OPT_INPUTFQ2} ]; then
    SEQ_MODE="Single End"
fi

Cutadapt=cutadapt
NAME_BASE1=$(basename ${OPT_INPUTZIP1%.zip})
QC_DIR1=${OPT_OUT_DIR}/${NAME_BASE1}
unzip ${OPT_INPUTZIP1} -d ${OPT_OUT_DIR}
failnum1=$(qc_judge ${QC_DIR1}/summary.txt)
if (( failnum1 > 4 )); then
    echo "FAIL: Too many failures of $OPT_INPUTZIP1"
    exit 0
fi

if [ "${SEQ_MODE}" = "Paired End" ]; then
    NAME_BASE2=$(basename ${OPT_INPUTZIP1%.zip})
    QC_DIR2=${OPT_OUT_DIR}/${NAME_BASE2}
    unzip ${OPT_INPUTZIP2} -d ${OPT_OUT_DIR}
    failnum2=$(qc_judge ${QC_DIR2}/summary.txt)
    if (( failnum2 > 4 )); then
        echo "FAIL: Too many failures of $OPT_INPUTZIP2"
        exit 0
    fi
fi


if [ $(qc_grep_status ${QC_DIR1}/summary.txt "Per base sequence quality") != 'FAIL' ]; then
    getblock ">>Per base sequence quality" ">>END_MODULE" ${QC_DIR1}/fastqc_data.txt > ${OPT_OUT_DIR}/${NAME_BASE1}_base_quality.txt
    INFO_SEQ_LENGTH1=$(cat ${QC_DIR1}/fastqc_data.txt | grep "Sequence length" | cut -f 2)
    INFO_SEQ_LENGTH1=${INFO_SEQ_LENGTH1##*-}
    MEAN_LENGTH=${INFO_SEQ_LENGTH1}
    if [ "${SEQ_MODE}" == "Paired End" ]; then
        getblock ">>Per base sequence quality" ">>END_MODULE" ${QC_DIR2}/fastqc_data.txt > ${OPT_OUT_DIR}/${NAME_BASE2}_base_quality.txt
        INFO_SEQ_LENGTH2=$(cat ${QC_DIR2}/fastqc_data.txt | grep "Sequence length" | cut -f 2)
        INFO_SEQ_LENGTH2=${INFO_SEQ_LENGTH2##*-}
        MEAN_LENGTH=$(((INFO_SEQ_LENGTH1 + INFO_SEQ_LENGTH2)/2))
    fi
    if [ "${SEQ_MODE}" = "Single End" ]; then
        OUT_FQ1_NAME=$(basename $OPT_INPUTFQ1)
        OUT_FQ1_NAME=${OUT_FQ1_NAME/fastq/qc.fastq}
        $Cutadapt -q 20,20 -m $((MEAN_LENGTH/2)) \
                  -o ${OPT_OUT_DIR}/${OUT_FQ1_NAME} ${OPT_INPUTFQ1}
    elif [ "${SEQ_MODE}" = "Paired End" ]; then
        OUT_FQ1_NAME=$(basename $OPT_INPUTFQ1)
        OUT_FQ1_NAME=${OUT_FQ1_NAME/fastq/qc.fastq}
        OUT_FQ2_NAME=$(basename $OPT_INPUTFQ2)
        OUT_FQ2_NAME=${OUT_FQ2_NAME/fastq/qc.fastq}
        $Cutadapt -q 20,20 -m $((MEAN_LENGTH/2)) \
                  -o ${OPT_OUT_DIR}/${OUT_FQ1_NAME} -p ${OPT_OUT_DIR}/${OUT_FQ2_NAME} \
                  ${OPT_INPUTFQ1} ${OPT_INPUTFQ2}
    fi
else
    ln -s ${OPT_INPUTFQ1} ${OPT_OUT_DIR}/${OUT_FQ1_NAME}
    if [ "${SEQ_MODE}" = "Paired End" ]; then
        ln -s ${OPT_INPUTFQ2} ${OPT_OUT_DIR}/${OUT_FQ2_NAME}
    fi
fi
