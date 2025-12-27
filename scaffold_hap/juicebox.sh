#!/bin/bash

set -euo pipefail

fa=""
agp=""
pairs=""
output_prefix=""
GPhase_path=""

# Usage instructions
usage() {
    cat << EOF
|
|Generate a Hi-C heatmap
|Usage: $(basename "$0") -f <asm.rename.fa> -a <gphase.rename.agp> -p <map.rename.pairs> -o <prefix> -g <GPhase_path> [options]
|
| Required:
|  -f, --fa FILE         FASTA file after copying the collapsed unitigs. (required)
|  -a, --agp FILE        AGP file after copying the collapsed unitigs. (required)
|  -p, --pairs FILE      Pairs file after copying the collapsed unitigs. (required)
|  -o, --output PREFIX   Output prefix for .hic and intermediate files (required)
|  -g, --gphase PATH     Path to GPhase installation directory (required)
|
|Example:
|  $(basename "$0") -f rename.fa -a rename.agp -p rename.pairs -o sample_hic -g /Path/to/GPhase
|
EOF
    exit 0
}

# If no arguments are provided, display help directly
[[ $# -eq 0 ]] && usage
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--fa)
            fa="$2"
            shift 2
            ;;
        -a|--agp)
            agp="$2"
            shift 2
            ;;
        -p|--pairs)
            pairs="$2"
            shift 2
            ;;
        -o|--output)
            output_prefix="$2"
            shift 2
            ;;
        -g|--gphase)
            GPhase_path="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option: $1" >&2
            echo "Use -h or --help for usage information."
            exit 1
            ;;
    esac
done

# Check if all required parameters are provided
for var in fa agp pairs output_prefix GPhase_path; do
    if [[ -z "${!var}" ]]; then
        echo "Error: Missing required parameter --$(echo $var | tr '[:upper:]' '[:lower:]' | tr _ -)" >&2
        echo "Use -h or --help for usage information."
        exit 1
    fi
done

# Check if input files exist
for file in "$fa" "$agp" "$pairs"; do
    [[ -f "$file" ]] || { echo "Error: File not found: $file"; exit 1; }
done

# Check if GPhase directory exists
[[ -d "$GPhase_path" ]] || { echo "Error: GPhase_path directory not found: $GPhase_path"; exit 1; }
[[ -d "$GPhase_path/src/HapHiC/utils" ]] || { echo "Error: GPhase_path/src/HapHiC/utils directory not found"; exit 1; }

# Locate juicer executable and juicer_tools.jar
juicer_path=$(find "$GPhase_path/src/HapHiC/utils" -maxdepth 1 -name 'juicer*' -type f -executable 2>/dev/null | head -n1)
juicer_tools_path=$(find "$GPhase_path/src/HapHiC/utils" -maxdepth 1 -name 'juicer_tools*' -type f 2>/dev/null | head -n1)

[[ -x "$juicer_path" ]] || { echo "Error: juicer executable not found in $GPhase_path/src/HapHiC/utils/"; exit 1; }
[[ -f "$juicer_tools_path" ]] || { echo "Error: juicer_tools.jar not found in $GPhase_path/src/HapHiC/utils/"; exit 1; }

log_file="${output_prefix}.log"
hic_part="${output_prefix}.hic.part"
hic_final="${output_prefix}.hic"
txt_file="${output_prefix}.txt"

echo "Starting conversion to Juicer .hic format..."
echo "Input pairs: $pairs"
echo "Output prefix: $output_prefix"
echo "Log file: $log_file"
echo "Final .hic: $hic_final"

# Step 1: Run juicer pre
samtools faidx ${fa}
fai=${fa}.fai
"$juicer_path" pre \
    -a -q 0 \
    -o "$output_prefix" \
    "$pairs" "$agp" "$fai" \
    --file-type pa5 \
    > "$log_file" 2>&1

# Check if PRE_C_SIZE exists in the log
if ! grep -q '^PRE_C_SIZE' "$log_file"; then
    echo "Error: PRE_C_SIZE not found in log. juicer pre may have failed."
    echo "Check log: $log_file"
    exit 1
fi

# Extract chromosome sizes
chrom_sizes=$(awk '/^PRE_C_SIZE/ {print $2 " " $3}' "$log_file")

# Step 2: Run juicer_tools pre to generate .hic file
echo "Generating .hic file with juicer_tools..."
java -jar -Xmx64g "$juicer_tools_path" pre \
    "$txt_file" \
    "$hic_part" \
    <(echo "$chrom_sizes")

# Rename the final file
if [[ -f "$hic_part" ]]; then
    mv "$hic_part" "$hic_final"
    echo "Successfully generated: $hic_final"
else
    echo "Error: Failed to generate .hic.part file."
    exit 1
fi
echo "HiC heatmap generated successfully!"