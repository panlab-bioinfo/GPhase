#!/bin/bash

# Function for displaying usage
usage() {
    echo "-------------------------------------------------------------------------------------------------------"
    echo "|Usage: $0 -f <fa_file>  -p <output_prefix> -t <threads>"
    echo "|"
    echo "|Required Parameters:"
    echo "|  -f     <fa_file>                : The FASTA file containing the genome sequences."
    echo "|  -p     <output_prefix>          : The prefix for the output files."
    echo "|  -t     <threads>                : The number of threads (default: 32)."
    echo "|"
    echo "|Example:"
    echo "|  bash $0 -f genome.fa -p output_prefix -t 32"
    echo "--------------------------------------------------------------------------------------------------------"
    exit 1
}

# Initialize variables
fa_file=""
output_prefix=""
threads=""
default_threads=32

TEMP=$(getopt -o f:p:t:: -- "$@")

if [ $? != 0 ]; then
    echo "Error: Invalid arguments."
    usage
fi

eval set -- "$TEMP"

while true; do
    case "$1" in
        -f) fa_file="$2"; shift 2 ;;
        -p) output_prefix="$2"; shift 2 ;;
        -t) threads="$2"; shift 2 ;;
        --) shift; break ;;
        *) usage ;;
    esac
done

# If threads is not set, use the default value
if [ -z "$threads" ]; then
    threads=$default_threads
fi

# Validate that all required arguments are provided
if [ -z "$fa_file" ] || [ -z "$output_prefix" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if input files exist
if [ ! -f "$fa_file" ]; then
    echo "Error: File '$file' does not exist."
    exit 1
fi

log_path=$(pwd)
log_file="${log_path}/popCNV_pipeline.log"

LOG_INFO() {
    time=$(date "+%Y-%m-%d %H:%M:%S")
    log_file=$1
    flag=$2
    msg=$3
    echo "${time} <popCNV_pipeline> [${flag}] ${msg}" >> ${log_file}

}

# Directory setup
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
LOG_INFO ${log_file} "path" "Script dir : ${SCRIPT_DIR}"

LOG_INFO ${log_file} "run" "samtools faidx ${fa_file}"
samtools faidx ${fa_file}

LOG_INFO ${log_file} "run" "Mapping HiFi data to ${fa_file} using minimap2"
minimap2 -t ${threads} -ax map-hifi ${fa_file} *.fastq.gz \
	| samtools view -bhS -t ${fa_file}.fai \
    | samtools sort -@ ${threads} > ${output_prefix}.bam

LOG_INFO ${log_file} "run" "samtools index -@ ${threads} ${output_prefix}.bam"
samtools index -@ ${threads} ${output_prefix}.bam

LOG_INFO ${log_file} "run" "Creat directory: bam_files"
mkdir bam_files && cd bam_files
ln -s ../${output_prefix}.bam
ln -s ../${output_prefix}.bam.bai
cd ../


# gene.list
LOG_INFO ${log_file} "run" "Creat gene.list"
awk '{print $1"\t"1"\t"$2"\t"$1}' ${fa_file}.fai > gene.list

# group.list
LOG_INFO ${log_file} "run" "Creat group.list"
echo -e "${output_prefix}\t${output_prefix}" > group.list

LOG_INFO ${log_file} "run" "Run the popcnv pipeline: python ${SCRIPT_DIR}/../src/popCNV/popCNV.py -g ${fa_file} -s 1000  -b bam_files/ -l gene.list -w ./popcnv --group group.list  --wild 0 -t ${threads}"
python ${SCRIPT_DIR}/../src/popCNV/popCNV.py -g ${fa_file} -s 1000  -b bam_files/ -l gene.list -w ./popcnv --group group.list  --wild 0 -t ${threads}
