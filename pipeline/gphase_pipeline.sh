#!/bin/bash

# Function for displaying usage
usage() {
    echo "| "
    echo "|Gphase: A phasing assembly tool using assembly graphs and Hi-C data"
    echo "|    Usage: /GPhase/to/path/gphase pipeline -f <fa_file> -g <gfa> -c <collapse_num_file> -m <map_file> --n_chr <n_chr> --n_hap <n_hap> -p <output_prefix>"
    echo "|"
    echo "|>>> Required Parameters:"
    echo "|  -f                 <fa_file>                : The FASTA file containing the genome sequences."
    echo "|  -g                 <gfa>                    : The GFA file representing the assembly graph."
    echo "|  -c                 <collapse_num_file>      : The file that number information for collapse unitigs."
    echo "|  -m                 <map_file>               : The mapping file used to map the Hi-C reads."
    echo "|  -p                 <output_prefix>          : The prefix for the output files."
    echo "|  --n_chr            <n_chr>                  : The number of chromosomes (integer)."
    echo "|  --n_hap            <n_hap>                  : The number of haplotypes (integer)."
    echo "|"
    echo "|>>> Optional parameters:"
    echo "|  -e                 <enzyme_site>            : The restriction enzyme cutting site, default: GATC."
    echo "|"
    echo "|>>> preprocessing Parameters:"
    echo "|  --cluster_q        <cluster_q>              : Filtered mapQ value when using HiC in the clustering step, default: 1"
    echo "|  --scaffold_q       <scaffold_q>             : Filter mapQ value when using HiC in the scaffolding step, default: 0"
    echo "|"
    echo "|>>> clustering chromosomes Parameters:"
    echo "|  --split_gfa_n      <split_gfa_n>            : Number of common neighbors when splitting GFA, default: 5"
    echo "|  --chr_pm           <partig_chr_pm>          : Similarity of partig when clustering chr, default: 0.9"
    echo "|"
    echo "|>>> clustering haplotypes Parameters:"
    echo "|  --hap_pm           <partig_hap_pm>          : Similarity of partig when clustering hap, default: 0.60"
    echo "|  --expand           <resexpandcue>           : Whether to expand the allele, default: False."
    echo "|  --rescue           <rescue>                 : Whether to rescue the subgraph, default: False."
    echo "|  --reassign_number  <reassign_number>        : Number of reassign step, default: 1. [1-3]"
    echo "|"
    echo "|>>> scaffolding haplotypes Parameters:"
    echo "|  --thread           <thread>                 : Number of parallel processes, default: 12."
    echo "|  --no_contig_ec     <no_contig_ec >          : do not do contig error correction in YaHS, default: False."
    echo "|  --no_scaffold_ec   <no_scaffold_ec >        : do not do scaffold error correction in YaHS, default: False."
    echo "|  --min_len          <min_len>                : minimum scaffold length(kb) in haphic sort, default: 0. [0-1000]"
    echo "|  --mutprob          <mutprob>                : mutation probability in the genetic algorithm in haphic sort, default: 0.6. [0.1-0.9]"
    echo "|  --ngen             <ngen>                   : number of generations for convergence in haphic sort, default: 20000."
    echo "|  --npop             <npop>                   : mopulation size in haphic sort, default: 200."
    echo "|  --processes        <processes>              : processes for fast sorting and ALLHiC optimization, default: 32."
    echo "|  -h, --help         Show this help message"
    echo "|"
    echo "|Example:"
    echo "|  /GPhase/to/path/gphase pipeline -f genome.fa -g genome.bp.p_utg.gfa -c collapse_num.txt -m map_file.bam --n_chr 12 --n_hap 4 -p output_prefix"
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
enzyme_site="GATC"
cluster_q=1
scaffold_q=0
split_gfa_n="5"
chr_pm="0.95"
hap_pm="0.60"
thread="12"
no_contig_ec=""
no_scaffold_ec=""
min_len="0" 
mutprob="0.6"   
ngen="20000"
npop="200"
processes="32"




TEMP=$(getopt -o f:g:c:m:p:e:h --long n_chr:,n_hap:,f:,g:,c:,m:,p:e:,cluster_q:,scaffold_q:,chr_pm:,hap_pm:,split_gfa_n:,rescue,expand,reassign_number:,thread:,no_contig_ec,no_scaffold_ec,min_len:,mutprob:,ngen:,processes:,help -- "$@")

if [ $? != 0 ]; then
    echo "Error: Invalid arguments."
    usage
fi

eval set -- "$TEMP"


while true; do
    case "$1" in
        -f|--f) fa_file="$2"; shift 2 ;;
        -g|--g) gfa="$2"; shift 2 ;;
        -c|--c) collapse_num_file="$2"; shift 2 ;;
        -m|--m) map_file="$2"; shift 2 ;;
        --n_chr) n_chr="$2"; shift 2 ;;
        --n_hap) n_hap="$2"; shift 2 ;;
        -p) output_prefix="$2"; shift 2 ;;
        -e) enzyme_site="$2"; shift 2 ;;
        --cluster_q) cluster_q="$2"; shift 2 ;;    
        --scaffold_q) scaffold_q="$2"; shift 2 ;;
        --split_gfa_n)
            if [[ "$2" =~ ^[2-9]+$ ]]; then 
                split_gfa_n="$2"
            else
                echo "Error: --split_gfa_n must be an integer between 2 and 9."
                usage
            fi
            shift 2 ;;
        --reassign_number) 
            if [[ "$2" =~ ^[1-3]+$ ]]; then
                reassign_number="$2"
            else
                echo "Error: --reassign_number must be an integer between 1 and 3."
                usage
            fi
            shift 2 ;;
        --chr_pm)
            if [ "$(echo "$2 >= 0.5 && $2 < 1" | bc -l)" -eq 1 ]; then
                chr_pm="$2"
            else
                echo "Error: --chr_pm must be a float between 0.5 and 1."
                usage
            fi
            shift 2 ;;
        --hap_pm)
            if [ "$(echo "$2 >= 0.5 && $2 < 1" | bc -l)" -eq 1 ]; then
                hap_pm="$2"
            else
                echo "Error: --hap_pm must be a float between 0.5 and 1."
                usage
            fi
            shift 2 ;;
        --rescue) rescue="--rescue"; shift ;;
        --expand) expand="--expand"; shift ;;
        --thread) thread="$2"; shift 2 ;;
        --no_contig_ec) no_contig_ec="--no_contig_ec";shift ;;
        --no_scaffold_ec) no_scaffold_ec="--no_scaffold_ec";shift ;;
        --min_len)
            if [[ "$2" =~ ^([0-9]{1,3}|1000)$ ]]; then 
                min_len="$2"
            else
                echo "Error: --min_len must be an int between 0 and 1000."
                usage
            fi
            shift 2 ;;
        --mutprob) 
            pattern='^0\.[1-9]*[1-9]+$'  # 允许0.1-0.9
            if [[ "$2" =~ $pattern ]]; then
                mutprob="$2"
                else
                    echo "Error: --mutprob must be a float between 0.1 and 0.9."
                    usage
                fi
                shift 2 ;;
        --ngen) ngen="$2"; shift 2 ;;
        --npop) npop="$2"; shift 2 ;;
        --processes) processes="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) usage ;;
    esac
done

set -e 

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
    if [ ! -f "$(basename "$file")" ]; then
        ln -s "$file" "$(basename "$file")"
    fi
done

if (( $(echo "$cluster_q < $scaffold_q" | bc -l) )); then
    echo "Error: --cluster_q ($cluster_q) must be greater than or equal to --scaffold_q ($scaffold_q)."
    exit 1
fi

if [[ -z "$reassign_number" ]]; then
    reassign_number=1
fi


fa_file="$(basename "$fa_file")"
gfa="$(basename "$gfa")"
collapse_num_file="$(basename "$collapse_num_file")"
map_file="$(basename "$map_file")"

current_dir=$(pwd)

log_path=$(pwd)
log_file="${log_path}/GPhase_pipeline.log"

LOG_INFO() {
    time=$(date "+%Y-%m-%d %H:%M:%S")
    log_file=$1
    flag=$2
    msg=$3
    echo "${time} <GPhase_pipeline> [${flag}] ${msg}" >> ${log_file}

}

run_step() {
    local cmd="$1"
    local step_name="$2"
    LOG_INFO "${log_file}" "run" "Running ${step_name}..."
    eval "${cmd}"
    local status=$?
    if [ $status -ne 0 ]; then
        LOG_INFO "${log_file}" "err" "Error: ${step_name} failed with exit code $status."
        exit $status
    fi
}

# Directory setup
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
LOG_INFO ${log_file} "path" "Script dir : ${SCRIPT_DIR}"

# Set up logging
LOG_INFO ${log_file} "start" "Pipeline start"

mkdir -p gphase_output && cd gphase_output
LOG_INFO ${log_file} "run" "Created output directory: gphase_output"

mkdir -p preprocessing && cd preprocessing
LOG_INFO ${log_file} "run" "Created output directory: preprocessing"

# Symlink files
ln -s "../../${fa_file}"
ln -s "../../${map_file}"

# Step 1: Run get_RE.py
run_step "python ${SCRIPT_DIR}/../cluster_chr/get_RE.py -f ${fa_file} -e ${enzyme_site} -op ${output_prefix}"

if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: get_RE.py failed."
    exit 1
fi

# Generate the pairs from bam
samtools faidx ${fa_file}
awk 'BEGIN{print "## pairs format v1.0.0\n##columns: readID chrom1 pos1 chrom2 pos2 mapQ"}{print "##chromsize: "$1" "$2}' ${fa_file}.fai > map.pairs
samtools view -@ 8 "${map_file}" | awk -v scaffold_q="${scaffold_q}" '{if($7=="=")$7=$3;if($5>=scaffold_q)print $1"\t"$3"\t"$4"\t"$7"\t"$8"\t"$5}' >> map.pairs



# get hic links from pairs
run_step "python  ${SCRIPT_DIR}/../cluster_chr/get_links.py -i map.pairs -o ${output_prefix} -q ${cluster_q}"

run_step "python  ${SCRIPT_DIR}/../cluster_chr/nor_hic.py -f ${output_prefix}.map.links.csv -r ${output_prefix}.RE_counts.txt -o ${output_prefix}.map.links.nor.csv"

if [ $? -ne 0 ]; then
    LOG_INFO ${log_file} "err" "Error: nor_hic.py failed."
    exit 1
fi

RE_file=${output_prefix}.RE_counts.txt
hic_links=${output_prefix}.map.links.nor.csv

# Step 2: Cluster chromosomes
LOG_INFO ${log_file} "run" "Created output directory: cluster_chr"
cd ../ && mkdir -p cluster_chr && cd cluster_chr
ln -s "../../${fa_file}"
ln -s "../../${gfa}"
ln -s "../preprocessing/${output_prefix}.RE_counts.txt"
ln -s "../preprocessing/${output_prefix}.map.links.nor.csv"


run_step "python ${SCRIPT_DIR}/../cluster_chr/cluster_chr.py -f ${fa_file} -r ${RE_file} -l ${hic_links} -op ${output_prefix} -n_chr ${n_chr} -g ${gfa}  -pm ${chr_pm} --split_gfa_n ${split_gfa_n}"


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


if [[ -z "$rescue" ]]; then
    rescue=""
    cr=${output_prefix}.chr.cluster.ctg.txt
else
    cr="rescue.cluster.ctg.txt"
fi

run_step "python ${SCRIPT_DIR}/../cluster_hap/cluster_hap.py -f ${fa_file} -r ${RE_file} -l ${hic_links} -op ${output_prefix} -n_chr ${n_chr} -n_hap ${n_hap} --collapse_num_file ${collapse_num_file} -d ${output_prefix}.digraph.csv -s group_ctgs_All.txt -c ${output_prefix}.chr.cluster.ctg.txt -cr ${cr} -pm ${hap_pm} --reassign_number ${reassign_number} ${rescue} ${expand}"


# Step 4: Scaffold haplotypes
LOG_INFO ${log_file} "run" "Created output directory: scaffold_hap"
cd ../ && mkdir -p scaffold_hap && cd scaffold_hap
ln -s "../cluster_chr/${fa_file}"
ln -s "../cluster_chr/${RE_file}"
ln -s "../cluster_chr/${hic_links}"
ln -s "../cluster_chr/${gfa}"
ln -s "../cluster_chr/group_ctgs_All.txt"
ln -s "../cluster_chr/${output_prefix}.digraph.csv"
ln -s "../preprocessing/map.map.pairs"

run_step "python ${SCRIPT_DIR}/../scaffold_hap/scaffold_hap_v2.py  -f ${fa_file} -r ${RE_file} -l ${hic_links} -op ${output_prefix} -n_chr ${n_chr} -n_hap ${n_hap} -CHP ../cluster_hap -s group_ctgs_All.txt -g ${gfa} -d ${output_prefix}.digraph.csv -m map.pairs -t ${thread} ${no_contig_ec} ${no_scaffold_ec} --min_len ${min_len} --mutprob ${mutprob} --ngen ${ngen} --npop ${npop} --processes ${processes}"

# Step 5: Final Output
# LOG_INFO ${log_file} "run" "Created output directory: gphase_final"
# cd ../ && mkdir -p gphase_final && cd gphase_final
# ln -s "../preprocessing/${RE_file}"
# ln -s "../preprocessing/${hic_links}"
# ln -s "../cluster_chr/${output_prefix}.rmTip.split.gfa"
# ln -s "../cluster_chr/group_ctgs_All.txt" ${output_prefix}.subgraphs.txt
# ln -s "../cluster_chr/rescue.cluster.ctg.txt" ${output_prefix}.chr.cluster.txt
# ln -s ../scaffold_hap/HapHiC_sort/scaffolds.sort.fa ${output_prefix}.scaffolds.fasta
# ln -s ../scaffold_hap/HapHiC_sort/final_agp/${output_prefix}.final.agp
