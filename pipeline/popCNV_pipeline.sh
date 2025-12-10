#!/bin/bash

# Function for displaying usage
usage() {
    echo "|"
    echo "|The popcnv software was run to estimate the collapse number of unitigs based on unitigs fasta file and HiFi reads"
    echo "|    Usage: $0 -f <fa_file>  -r <reads1.fastq.gz [reads2.fastq.gz ...]> -p <output_prefix> -t <threads>"
    echo "|"
    echo "|Required Parameters:"
    echo "|  -f     <fa_file>                : The FASTA file containing the utg sequences."
    echo "|  -p     <output_prefix>          : The prefix for the output files."
    echo "|  -t     <threads>                : The number of threads (default: 32)."
    echo "|  -r     <reads_files>            : One or more FASTQ(.gz) files, if there are multiple hifi reads files, the -r parameter needs to be placed last."
    echo "|"
    echo "|Example:"
    echo "|  bash $0 -f asm.fa -p output_prefix -t 32 -r hifi_read1.fq.gz"
    exit 1
}

# Initialize variables
fa_file=""
reads_files=()
output_prefix=""
threads=""
default_threads=32

r_flag=0

while getopts ":f:p:t:r" opt; do
    case $opt in
        f) fa_file="$OPTARG" ;;
        p) output_prefix="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        r)
            r_flag=1 ;;
        ?) usage; exit 1 ;;
    esac
done

shift $((OPTIND - 1))
if [[ $r_flag -eq 1 ]]; then
    while [[ $# -gt 0 && "$1" != -* ]]; do
        reads_files+=("$1")
        shift
    done
fi

# Validate required arguments
if [ -z "$fa_file" ] || [ -z "$output_prefix" ] || [ ${#reads_files[@]} -eq 0 ]; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi

# Check fasta exists
if [ ! -f "$fa_file" ]; then
    echo "Error: File '$fa_file' does not exist."
    exit 1
fi

set -e
log_path=$(pwd)
log_file="${log_path}/popCNV_pipeline.log"

LOG_INFO() {
    time=$(date "+%Y-%m-%d %H:%M:%S")
    log_file=$1
    flag=$2
    shift 2
    echo "${time} <popCNV_pipeline> [${flag}] $* " >> "$log_file"
}


# Directory setup
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
LOG_INFO ${log_file} "path" "Script dir : ${SCRIPT_DIR}"

LOG_INFO ${log_file} "run" "samtools faidx ${fa_file}"
samtools faidx ${fa_file}

LOG_INFO ${log_file} "info" "input FASTA files: ${fa_file}"
LOG_INFO ${log_file} "info" "input HiFi files: ${reads_files[@]}"

# Mapping
LOG_INFO ${log_file} "run" "Mapping HiFi data to ${fa_file} using minimap2"
minimap2 -t ${threads} -ax map-hifi ${fa_file} "${reads_files[@]}" \
    | samtools view -bhS -t ${fa_file}.fai \
    | samtools sort -@ ${threads} > ${output_prefix}.bam

LOG_INFO ${log_file} "run" "samtools index -@ ${threads} ${output_prefix}.bam"
samtools index -@ ${threads} ${output_prefix}.bam

# Organize bam
LOG_INFO ${log_file} "run" "Creat directory: bam_files"
mkdir -p bam_files && cd bam_files
ln -sf ../${output_prefix}.bam
ln -sf ../${output_prefix}.bam.bai
cd ..

# gene.list
LOG_INFO ${log_file} "run" "Creat gene.list"
awk '{print $1"\t"1"\t"$2"\t"$1}' ${fa_file}.fai > gene.list

# group.list
LOG_INFO ${log_file} "run" "Creat group.list"
echo -e "${output_prefix}\t${output_prefix}" > group.list

LOG_INFO ${log_file} "run" "Run the popcnv pipeline"
python ${SCRIPT_DIR}/../src/popCNV/popCNV.py \
    -g ${fa_file} -s 1000 \
    -b bam_files/ -l gene.list \
    -w ./popcnv --group group.list \
    --wild 0 -t ${threads}
