#!/usr/bin/env bash
set -euo pipefail

############################################
# Logging function
############################################
log_path=$(pwd)
log_file="${log_path}/PPL_pipeline.log"

LOG_INFO() {
    local time=$(date "+%Y-%m-%d %H:%M:%S")
    local log_file="${1:-}" 
    local flag="${2:-INFO}"  
    local msg="${3:-}"  

    if [[ -z "$log_file" ]]; then
        echo "${time} <PPL_pipeline> [${flag}] ${msg}" >&2
    else
        echo "${time} <PPL_pipeline> [${flag}] ${msg}" >> "$log_file"
    fi
}

############################################
# Print usage
############################################
usage() {
    echo "|"
    echo "|Run the PPL software using the fasta file and the fastq.gz file"
    echo "|Usage: $0 -j ./PPL-0.1.1.jar -g <genome.fa> -f <reads.fq> [options]"
    echo "|"
    echo "| Required:"
    echo "|  -j <jar_path>       PPL jar file (required)"
    echo "|  -g <genome.fa>      Genome FASTA file (required)"
    echo "|  -f <reads.fq>       Fastq file (required)"
    echo "|"
    echo "| Optional:"
    echo "|  -s <site>           Restriction enzyme site (default: GATC)"
    echo "|  -o <prefix>         Output prefix (default: PPL)"
    echo "|  -t <threads>        Threads (default: 12)"
    echo "|  -q <mapq>           MAPQ cutoff (default: 1)"
    echo "|  -h                  Show help and exit"
    echo "|"
    echo "|Example:"
    echo "|  bash $0 -j ./PPL-0.1.1.jar -f asm.fa -r reads.fq.gz -p PPL -t 32 -q 1"
    exit 1
}

############################################
# Default optional parameter values
############################################
site="^GATC"
output_prefix="PPL"
threads=12
cutoffMapq=1

############################################
# Parse options
############################################
while getopts "j:g:f:s:o:t:q:h" opt; do
    case $opt in
        j) jar="$OPTARG" ;;
        g) genome="$OPTARG" ;;
        f) fq_file="$OPTARG" ;;
        s) site="$OPTARG" ;;
        o) output_prefix="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        q) cutoffMapq="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

############################################
# Check required parameters
############################################
if [[ -z "${jar:-}" ]] || [[ -z "${genome:-}" ]] || [[ -z "${fq_file:-}" ]]; then
    LOG_INFO ${log_file} "error" "ERROR: -j <jar_path> -g <genome.fa> and -f <reads.fq> are required."
    usage
fi

############################################
# Print final parameter settings
############################################
LOG_INFO ${log_file} "args" "========== PARAMETERS =========="
LOG_INFO ${log_file} "args" "PPL jar path      : $jar"
LOG_INFO ${log_file} "args" "Genome file       : $genome"
LOG_INFO ${log_file} "args" "Fastq file        : $fq_file"
LOG_INFO ${log_file} "args" "Restriction site  : $site"
LOG_INFO ${log_file} "args" "Output prefix     : $output_prefix"
LOG_INFO ${log_file} "args" "Threads           : $threads"
LOG_INFO ${log_file} "args" "MAPQ cutoff       : $cutoffMapq"
LOG_INFO ${log_file} "================================"
LOG_INFO ${log_file} ""

############################################
# Main pipeline code with logging
############################################
LOG_INFO ${log_file} "info" "Running utils.VirDigestTool ..."
clean_site=${site//[\^\$]/}
java -Xmx64g -cp ${jar} utils.VirDigestTool \
    "$genome" \
    "$site" \
    ${genome%%.*}.${clean_site}.res.bed 2>&1 | tee -a "$log_file"

LOG_INFO ${log_file} "info" "Running PPL.jar ..."
java -Xmx64g -jar ${jar} --ligation_type res \
    --genomefile "${genome}" \
    --fastq "${fq_file}" \
    --splitReads N --resRemove N --disRemove N \
    --output ./ \
    --prefix "${output_prefix}" \
    --skipmap N \
    --start_step 2 \
    --restrictionsiteFile "${genome%%.*}.${clean_site}.res.bed" \
    --thread "${threads}" \
    --cutoffMapq "${cutoffMapq}" \
    --filter res 2>&1 | tee -a "$log_file"

LOG_INFO ${log_file} "info" "Generating chromsizes ..."
samtools faidx "${genome}" 2>&1 | tee -a "$log_file"
cut -f1,2 "${genome}.fai" > "${genome}.chromsizes"
LOG_INFO ${log_file} "info" "Chromsizes file: ${genome}.chromsizes"
cd ${output_prefix}
ln -s ../${genome}.chromsizes
LOG_INFO ${log_file} "info" "Running utils.FilterHyper ..."
java -Xmx64g -cp ${jar} utils.FilterHyper ${output_prefix}.final.contacts ${output_prefix}.final.filtered.contacts ${genome}.chromsizes 1000000 0.85
LOG_INFO ${log_file} "info" "Running utils.Contact2Pairs ..."
java -Xmx64g -cp ${jar} utils.Contact2Pairs ${output_prefix}.final.filtered.contacts map.PPL.pairs
LOG_INFO ${log_file} "info" "All steps completed successfully!"