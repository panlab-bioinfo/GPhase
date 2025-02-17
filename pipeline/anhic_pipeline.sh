#!/bin/bash

# Function for displaying usage
usage() {
    echo "-------------------------------------------------------------------------------------------------------"
    echo "|Usage: $0 -f <fa_file> -g <gfa> -c <collapse_num_file> -m <map_file> -n_chr <n_chr> -n_hap <n_hap> -p <output_prefix>"
    echo "|"
    echo "|Required Parameters:"
    echo "|  -f     <fa_file>                : The FASTA file containing the genome sequences."
    echo "|  -g     <gfa>                    : The GFA file representing the assembly graph."
    echo "|  -c     <collapse_num_file>      : The file that number information for collapse unitigs."
    echo "|  -m     <map_file>               : The mapping file used to map the Hi-C reads."
    echo "|  -n_chr <n_chr>                  : The number of chromosomes (integer)."
    echo "|  -n_hap <n_hap>                  : The number of haplotypes (integer)."
    echo "|  -p     <output_prefix>          : The prefix for the output files."
    echo "|"
    echo "|Example:"
    echo "|  bash $0 -f genome.fa -g genome.bp.p_utg.gfa -c collapse_num.txt -m map_file.pairs -n 12 -h 4 -p output_prefix"
    echo "--------------------------------------------------------------------------------------------------------"
    exit 1
}

# Initialize variables
fa_file=""
gfa=""
collapse_num_file=""
map_file=""
n_chr=""
n_hap=""
output_prefix=""


TEMP=$(getopt -o f:g:c:m:p: --long n_chr:,n_hap:,f:,g:,c:,m:,p: -- "$@")

if [ $? != 0 ]; then
    echo "Error: Invalid arguments."
    usage
fi

eval set -- "$TEMP"


while true; do
    case "$1" in
        -f) fa_file="$2"; shift 2 ;;
        -g) gfa="$2"; shift 2 ;;
        -c) collapse_num_file="$2"; shift 2 ;;
        -m) map_file="$2"; shift 2 ;;
        --n_chr) n_chr="$2"; shift 2 ;;
        --n_hap) n_hap="$2"; shift 2 ;;
        -p) output_prefix="$2"; shift 2 ;;
        --) shift; break ;;
        *) usage ;;
    esac
done

# Validate that all required arguments are provided
if [ -z "$fa_file" ] || [ -z "$gfa" ] || [ -z "$collapse_num_file" ] || [ -z "$map_file" ] || [ -z "$n_chr" ] || [ -z "$n_hap" ] || [ -z "$output_prefix" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if input files exist
for file in "$fa_file" "$gfa" "$collapse_num_file" "$map_file"; do
    if [ ! -f "$file" ]; then
        echo "Error: File '$file' does not exist."
        exit 1
    fi
done

log_path=$(pwd)
log_file="${log_path}/pipeline.log"

LOG_INFO() {
    time=$(date "+%Y-%m-%d %H:%M:%S")
    log_file=$1
    flag=$2
    msg=$3
    echo "${time} <anhic_pipeline> [${flag}] ${msg}" >> ${log_file}

}
# Directory setup
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
LOG_INFO ${log_file} "path" "Script dir : ${SCRIPT_DIR}"

# Set up logging
LOG_INFO ${log_file} "start" "Pipeline start"

mkdir -p anhic_output && cd anhic_output
LOG_INFO ${log_file} "run" "Created output directory: anhic_output"

mkdir -p preprocessing && cd preprocessing
LOG_INFO ${log_file} "run" "Created output directory: preprocessing"

# Symlink files
ln -s "../../${fa_file}"
ln -s "../../${map_file}"

# Step 1: Run get_RE.py
LOG_INFO ${log_file} "run" "Running get_RE.py..."
python ${SCRIPT_DIR}/../cluster_chr/get_RE.py -f ${fa_file} -e GATC -op ${output_prefix}

if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: get_RE.py failed."
    exit 1
fi

# Generate the links file
grep '^#' -v ${map_file} | awk '{if($2 != $4){print $2,$4}}' | \
        uniq -c | awk '{print $2","$3","$1}' > ${output_prefix}.chromap.links.csv

python  ${SCRIPT_DIR}/../cluster_chr/nor_hic.py -f ${output_prefix}.chromap.links.csv -r ${output_prefix}.RE_counts.txt -o ${output_prefix}.chromap.links.nor.csv

if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: nor_hic.py failed."
    exit 1
fi

RE_file=${output_prefix}.RE_counts.txt
hic_links=${output_prefix}.chromap.links.nor.csv

# Step 2: Cluster chromosomes
LOG_INFO ${log_file} "run" "Created output directory: cluster_chr"
cd ../ && mkdir -p cluster_chr && cd cluster_chr
ln -s "../../${fa_file}"
ln -s "../../${gfa}"
ln -s "../preprocessing/${output_prefix}.RE_counts.txt"
ln -s "../preprocessing/${output_prefix}.chromap.links.nor.csv"

LOG_INFO ${log_file} "run" "Running cluster_chr.py..."
python ${SCRIPT_DIR}/../cluster_chr/cluster_chr.py -f ${fa_file} -r ${RE_file} -l ${hic_links} -op ${output_prefix} -n_chr ${n_chr} -g ${gfa}

if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: cluster_chr.py failed."
    exit 1
fi

# Step 3: Cluster haplotypes
LOG_INFO ${log_file} "run" "Created output directory: cluster_hap"
cd ../ && mkdir -p cluster_hap && cd cluster_hap

ln -s "../cluster_chr/group_ctgs_All.txt"
ln -s "../cluster_chr/rescue.cluster.ctg.txt"
ln -s "../../${collapse_num_file}"
ln -s "../cluster_chr/${fa_file}"
ln -s "../cluster_chr/${RE_file}"
ln -s "../cluster_chr/${hic_links}"
ln -s "../cluster_chr/${output_prefix}.digraph.csv"
ln -s "../cluster_chr/${output_prefix}.chr.cluster.ctg.txt"


LOG_INFO ${log_file} "run" "Running cluster_hap.py..."
python ${SCRIPT_DIR}/../cluster_hap/cluster_hap.py -f ${fa_file} -r ${RE_file} -l ${hic_links} -op ${output_prefix} -n_chr ${n_chr} -n_hap ${n_hap} --collapse_num_file ${collapse_num_file} -d ${output_prefix}.digraph.csv -s group_ctgs_All.txt -c ${output_prefix}.chr.cluster.ctg.txt -cr rescue.cluster.ctg.txt
if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: cluster_hap.py failed."
    exit 1
fi

# Step 4: Scaffold haplotypes
LOG_INFO ${log_file} "run" "Created output directory: scaffold_hap"
cd ../ && mkdir -p scaffold_hap && cd scaffold_hap
ln -s "../cluster_chr/${fa_file}"
ln -s "../cluster_chr/${RE_file}"
ln -s "../cluster_chr/${hic_links}"
ln -s "../cluster_chr/${gfa}"
ln -s "../cluster_chr/group_ctgs_All.txt"
ln -s "../cluster_chr/${output_prefix}.digraph.csv"
ln -s "../../${map_file}"


LOG_INFO ${log_file} "run" "Running scaffold_hap.py..."
python ${SCRIPT_DIR}/../scaffold_hap/scaffold_hap.py  -f ${fa_file} -r ${RE_file} -l ${hic_links} -op ${output_prefix} -n_chr ${n_chr} -n_hap ${n_hap} -CHP ../cluster_hap -s group_ctgs_All.txt -g ${gfa} -d ${output_prefix}.digraph.csv -m ${map_file}
if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: scaffold_hap.py failed."
    exit 1
fi

# Step 5: Final Output
LOG_INFO ${log_file} "run" "Created output directory: anhic_final"
cd ../ && mkdir -p anhic_final && cd anhic_final
ln -s "../preprocessing/${RE_file}"
ln -s "../preprocessing/${hic_links}"
ln -s "../cluster_chr/${output_prefix}.rmTip.split.gfa"
ln -s "../cluster_chr/group_ctgs_All.txt" ${output_prefix}.subgraphs.txt
ln -s "../cluster_chr/rescue.cluster.ctg.txt" ${output_prefix}.chr.cluster.txt
ln -s ../scaffold_hap/HapHiC_sort/scaffolds.sort.fa ${output_prefix}.genome.fasta

LOG_INFO ${log_file} "done" "Pipeline completed successfully!"
