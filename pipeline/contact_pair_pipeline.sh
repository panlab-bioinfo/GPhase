#!/usr/bin/env bash
set -euo pipefail

############################################
# Logging function
############################################
log_path=$(pwd)
log_file="${log_path}/contact_pair_pipeline.log"

: > "${log_file}"

LOG_INFO() {
    local time
    time=$(date "+%Y-%m-%d %H:%M:%S")
    local target="${1:-}"
    local flag="${2:-INFO}"
    local msg="${3:-}"

    if [[ -z "${target}" ]]; then
        echo "${time} <contact_pair_pipeline> [${flag}] ${msg}" >&2
    else
        echo "${time} <contact_pair_pipeline> [${flag}] ${msg}" >> "${target}"
    fi
}

############################################
# Locate concatemer2pe
############################################
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
concatemer2pe_py="${script_dir}/../src/HapHiC/utils/concatemer2pe.py"

############################################
# Print usage
############################################
usage() {
    echo "|"
    echo "|Map contact-pair long reads (Pore-C/CiFi) with minimap2 and convert them to a GPhase-ready BAM file"
    echo "|Usage: $0 <genome.fa> <reads1.fq.gz> [reads2.fq.gz ...] [options]"
    echo "|"
    echo "| Required:"
    echo "|  <genome.fa>             Genome FASTA file"
    echo "|  <reads*.fq.gz>          Raw contact-pair reads in FASTQ format"
    echo "|"
    echo "| Optional:"
    echo "|  -o <prefix>             Output directory/prefix (default: contact_pair)"
    echo "|                          Final BAM path: <prefix>/map.concatemer2pe.bam"
    echo "|  -t <threads>            Threads (default: 16)"
    echo "|  -x <preset>             minimap2 preset [map-ont/map-hifi] (default: map-ont)"
    echo "|  -q <mapq>               MAPQ cutoff for concatemer2pe (default: 0)"
    echo "|  -i <percent_identity>   Percent identity cutoff (default: 0)"
    echo "|  -l <alignment_length>   Alignment length cutoff (default: 0)"
    echo "|  -h                      Show help and exit"
    echo "|"
    echo "| Output:"
    echo "|  <prefix>/map.concatemer2pe.bam can be passed directly to:"
    echo "|  gphase pipeline -m <prefix>/map.concatemer2pe.bam"
    echo "|  gphase will convert BAM to preprocessing/map.pairs internally."
    echo "|"
    echo "| Example:"
    echo "|  bash $0 asm.fa porec.1.fq.gz porec.2.fq.gz -x map-ont -o contact_pair -t 32 -q 1"
    echo "|  bash $0 asm.fa cifi.1.fq.gz cifi.2.fq.gz -x map-hifi -o contact_pair -t 32 -q 1"
    exit 1
}

############################################
# Default optional parameter values
############################################
output_prefix="contact_pair"
threads=16
preset="map-ont"
mapq=0
percent_identity=0
alignment_length=0
genome=""
reads_files=()

is_non_negative_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

############################################
# Parse options
############################################
if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
    usage
fi

if [[ $# -lt 2 ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: <genome.fa> and at least one reads FASTQ file are required."
    usage
fi

genome="$1"
shift

while [[ $# -gt 0 && "$1" != -* ]]; do
    reads_files+=("$1")
    shift
done

while getopts "o:t:x:q:i:l:h" opt; do
    case $opt in
        o) output_prefix="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        x) preset="$OPTARG" ;;
        q) mapq="$OPTARG" ;;
        i) percent_identity="$OPTARG" ;;
        l) alignment_length="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

shift $((OPTIND - 1))

if [[ $# -ne 0 ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: unexpected extra arguments: $*"
    usage
fi

############################################
# Validate parameters
############################################
if [[ -z "${genome}" || ${#reads_files[@]} -eq 0 ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: <genome.fa> and at least one reads FASTQ file are required."
    usage
fi

if ! is_non_negative_integer "${threads}" || [[ "${threads}" -eq 0 ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: -t <threads> must be a positive integer."
    exit 1
fi

if ! is_non_negative_integer "${mapq}"; then
    LOG_INFO "${log_file}" "error" "ERROR: -q <mapq> must be a non-negative integer."
    exit 1
fi

if ! is_non_negative_integer "${alignment_length}"; then
    LOG_INFO "${log_file}" "error" "ERROR: -l <alignment_length> must be a non-negative integer."
    exit 1
fi

if ! awk -v x="${percent_identity}" 'BEGIN { exit !(x >= 0 && x <= 100) }'; then
    LOG_INFO "${log_file}" "error" "ERROR: -i <percent_identity> must be between 0 and 100."
    exit 1
fi

case "${preset}" in
    map-ont|map-hifi) ;;
    *)
        LOG_INFO "${log_file}" "error" "ERROR: unsupported minimap2 preset '${preset}'. Use map-ont or map-hifi."
        exit 1
        ;;
esac

for cmd in python3 samtools awk minimap2; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        LOG_INFO "${log_file}" "error" "ERROR: required command not found: ${cmd}"
        exit 1
    fi
done

if [[ ! -f "${concatemer2pe_py}" ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: concatemer2pe.py not found: ${concatemer2pe_py}"
    exit 1
fi

if [[ ! -f "${genome}" ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: genome FASTA not found: ${genome}"
    exit 1
fi

for fq_file in "${reads_files[@]}"; do
    if [[ ! -f "${fq_file}" ]]; then
        LOG_INFO "${log_file}" "error" "ERROR: FASTQ file not found: ${fq_file}"
        exit 1
    fi
done

############################################
# Print final parameter settings
############################################
LOG_INFO "${log_file}" "args" "========== PARAMETERS =========="
LOG_INFO "${log_file}" "args" "concatemer2pe path : ${concatemer2pe_py}"
LOG_INFO "${log_file}" "args" "Genome file        : ${genome}"
LOG_INFO "${log_file}" "args" "FASTQ files        : ${reads_files[*]}"
LOG_INFO "${log_file}" "args" "minimap2 preset    : ${preset}"
LOG_INFO "${log_file}" "args" "Output prefix      : ${output_prefix}"
LOG_INFO "${log_file}" "args" "Threads            : ${threads}"
LOG_INFO "${log_file}" "args" "MAPQ cutoff        : ${mapq}"
LOG_INFO "${log_file}" "args" "Percent identity   : ${percent_identity}"
LOG_INFO "${log_file}" "args" "Alignment length   : ${alignment_length}"
LOG_INFO "${log_file}" "args" "================================"
LOG_INFO "${log_file}" "info" ""

############################################
# Prepare outputs
############################################
mkdir -p "${output_prefix}"
output_name="$(basename "${output_prefix}")"
local_genome="${output_prefix}/$(basename "${genome}")"
unsorted_bam="${output_prefix}/${output_name}.concatemer.unsorted.bam"
paired_bam="${output_prefix}/map.concatemer2pe.bam"

ln -sfn "$(realpath "${genome}")" "${local_genome}"

############################################
# Main pipeline code with logging
############################################
LOG_INFO "${log_file}" "info" "Preparing local FASTA index in output directory ..."
samtools faidx "${local_genome}" 2>&1 | tee -a "${log_file}"

reference_bases=$(awk '{s+=$2} END {print s+0}' "${local_genome}.fai")
reference_size_gb=$(awk -v s="${reference_bases}" 'BEGIN {print int((s + 1000000000 - 1) / 1000000000)}')
minimap2_I=$(( reference_size_gb > 8 ? reference_size_gb : 8 ))

LOG_INFO "${log_file}" "info" "Detected reference bases: ${reference_bases}"
LOG_INFO "${log_file}" "info" "Rounded reference size (GB): ${reference_size_gb}"
LOG_INFO "${log_file}" "info" "Using minimap2 -I value: ${minimap2_I}G"
LOG_INFO "${log_file}" "info" "Running minimap2 and writing unsorted BAM ..."
minimap2 -t "${threads}" -I "${minimap2_I}G" -ax "${preset}" "${local_genome}" "${reads_files[@]}" \
    2>> "${log_file}" | \
    samtools view -@ "${threads}" -b -o "${unsorted_bam}" - \
    2>> "${log_file}"

LOG_INFO "${log_file}" "info" "Running concatemer2pe.py ..."
python3 "${concatemer2pe_py}" "${unsorted_bam}" \
    --output "${paired_bam}" \
    --mapq "${mapq}" \
    --percent-identity "${percent_identity}" \
    --alignment-length "${alignment_length}" \
    --threads "${threads}" 2>&1 | tee -a "${log_file}"

LOG_INFO "${log_file}" "info" "All steps completed successfully!"
LOG_INFO "${log_file}" "info" "GPhase-ready BAM file: ${paired_bam}"
LOG_INFO "${log_file}" "info" "Run downstream with: gphase pipeline -m ${paired_bam}"
