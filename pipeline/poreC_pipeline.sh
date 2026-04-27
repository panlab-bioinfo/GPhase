#!/usr/bin/env bash
set -euo pipefail

############################################
# Logging function
############################################
log_path=$(pwd)
log_file="${log_path}/poreC_pipeline.log"

: > "${log_file}"

LOG_INFO() {
    local time
    time=$(date "+%Y-%m-%d %H:%M:%S")
    local target="${1:-}"
    local flag="${2:-INFO}"
    local msg="${3:-}"

    if [[ -z "${target}" ]]; then
        echo "${time} <poreC_pipeline> [${flag}] ${msg}" >&2
    else
        echo "${time} <poreC_pipeline> [${flag}] ${msg}" >> "${target}"
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
    echo "|Map ONT concatemer long reads with minimap2 and convert them to a GPhase-ready BAM file"
    echo "|Usage: $0 <genome.fa> <reads.fq.gz> [options]"
    echo "|"
    echo "| Required:"
    echo "|  <genome.fa>             Genome FASTA file"
    echo "|  <reads.fq.gz>           Raw concatemer reads in FASTQ format"
    echo "|"
    echo "| Optional:"
    echo "|  -o <prefix>             Output directory/prefix (default: porec)"
    echo "|                          Final BAM path: <prefix>/map.concatemer2pe.bam"
    echo "|  -t <threads>            Threads (default: 16)"
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
    echo "|  bash $0 asm.fa porec.fq.gz -o concatemer2pe -t 32 -q 1"
    exit 1
}

############################################
# Default optional parameter values
############################################
output_prefix="porec"
threads=16
mapq=0
percent_identity=0
alignment_length=0
genome=""
fq_file=""

############################################
# Parse options
############################################
if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
    usage
fi

if [[ $# -lt 2 ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: <genome.fa> and <reads.fq.gz> are required."
    usage
fi

genome="$1"
fq_file="$2"
shift 2

while getopts "o:t:q:i:l:h" opt; do
    case $opt in
        o) output_prefix="$OPTARG" ;;
        t) threads="$OPTARG" ;;
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
if [[ -z "${genome}" || -z "${fq_file}" ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: <genome.fa> and <reads.fq.gz> are required."
    usage
fi

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

if [[ -n "${fq_file}" && ! -f "${fq_file}" ]]; then
    LOG_INFO "${log_file}" "error" "ERROR: FASTQ file not found: ${fq_file}"
    exit 1
fi

############################################
# Print final parameter settings
############################################
LOG_INFO "${log_file}" "args" "========== PARAMETERS =========="
LOG_INFO "${log_file}" "args" "concatemer2pe path : ${concatemer2pe_py}"
LOG_INFO "${log_file}" "args" "Genome file        : ${genome}"
LOG_INFO "${log_file}" "args" "FASTQ file         : ${fq_file}"
LOG_INFO "${log_file}" "args" "minimap2 preset    : map-ont"
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
unsorted_bam="${output_prefix}/${output_name}.concatemer.unsorted.bam"
paired_bam="${output_prefix}/map.concatemer2pe.bam"

############################################
# Main pipeline code with logging
############################################
LOG_INFO "${log_file}" "info" "Indexing FASTA for minimap2 -I estimation ..."
samtools faidx "${genome}" 2>&1 | tee -a "${log_file}"

reference_bases=$(awk '{s+=$2} END {print s+0}' "${genome}.fai")
reference_size_gb=$(awk -v s="${reference_bases}" 'BEGIN {print int((s + 1000000000 - 1) / 1000000000)}')
minimap2_I=$(( reference_size_gb > 8 ? reference_size_gb : 8 ))

LOG_INFO "${log_file}" "info" "Detected reference bases: ${reference_bases}"
LOG_INFO "${log_file}" "info" "Rounded reference size (GB): ${reference_size_gb}"
LOG_INFO "${log_file}" "info" "Using minimap2 -I value: ${minimap2_I}G"
LOG_INFO "${log_file}" "info" "Running minimap2 and writing unsorted BAM ..."
minimap2 -t "${threads}" -I "${minimap2_I}G" -ax map-ont "${genome}" "${fq_file}" \
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
