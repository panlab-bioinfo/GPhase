#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# ====== basic helpers & trap ======
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
current_dir=$(pwd)
log_path="${current_dir}"
log_file="${log_path}/gphase_pipeline.log"

: > "${log_file}"


timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

die() {
    echo "$(timestamp) <gphase_pipeline> [FATAL] $*" | tee -a "${log_file}" >&2
    exit 1
}
info() {
    echo "$(timestamp) <gphase_pipeline> [INFO] $*" | tee -a "${log_file}"
}
warn() {
    echo "$(timestamp) <gphase_pipeline> [WARN] $*" | tee -a "${log_file}"
}

_on_error() {
    local rc=$?
    local line="$1"
    echo "$(timestamp) <gphase_pipeline> [ERROR] Exit code ${rc} at line ${line}. Last command: ${BASH_COMMAND}" >> "${log_file}"
    echo "---- Stack trace ----" >> "${log_file}"
    local i=1
    while caller $i >> "${log_file}" 2>&1; do ((i++)); done
    exit "${rc}"
}
trap '_on_error ${LINENO}' ERR

check_file_exists_and_nonempty() {
    local path="$1"
    local desc="${2:-file}"

    if [[ -z "$path" ]]; then
        die "check_file_exists_and_nonempty called with empty path for $desc"
    fi

    if [[ -L "$path" ]]; then
        if [[ ! -e "$path" ]]; then
            die "Broken symlink: $path ($desc)"
        fi
        path="$(readlink -f "$path")"
    fi

    if [[ ! -f "$path" ]]; then
        die "Required $desc does not exist: $path"
    fi
    if [[ ! -s "$path" ]]; then
        die "Required $desc is empty: $path"
    fi
    info "Checked OK: $path ($desc)"
}

# safe symlink
safe_ln() {
    local src="$1"
    local dst="${2:-$(basename "$src")}"
    if [[ ! -e "$src" ]]; then
        warn "Source file for link not found: $src"
        return 1
    fi

    if [ -e "$dst" ] || [ -L "$dst" ]; then
        info "Skip link: ${dst} already exists"
    else
        ln -s "$src" "$dst"
        info "Linked: ${src} -> ${dst}"
    fi
}

run_step() {
    local cmd="$1"
    local step_name="${2:-command}"
    info "==== Start: ${step_name} ===="
    info "CMD: ${cmd}"
    eval "${cmd}" >> "${log_file}" 2>&1
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        warn "==== Failed: ${step_name} (Exit code ${rc}) ===="
    else
        info "==== Done: ${step_name} ===="
    fi
}

# usage
usage() {
    cat <<'USAGE' | sed 's/^/| /'
GPhase: A phasing assembly tool using assembly graph and Hi-C data
Usage: $(basename "$0") pipeline -f <fa_file> -g <gfa> -c <collapse_num_file> -m <map_file> --n_chr <n_chr> --n_hap <n_hap> -p <output_prefix>

>>> Required Parameters:
  -f                  <fa_file>                  : The FASTA file containing the genome sequences.
  -g                  <gfa>                      : The GFA file representing the assembly graph.
  -c                  <collapse_num_file>        : The file that number information for collapse unitigs.
  -m                  <map_file>                 : The mapping file used to map the Hi-C reads (bam or pairs).
  -p                  <output_prefix>            : The prefix for the output files. Only [a-zA-Z0-9.] allowed.
  --n_chr             <n_chr>                    : The number of chromosomes (integer, > 0).
  --n_hap             <n_hap>                    : The number of haplotypes (integer, > 0).

>>> Optional parameters:
  -e                  <enzyme_site>              : The restriction enzyme cutting site, default: GATC.

>>> preprocessing Parameters:
  --cluster_q         <cluster_q>                : Filtered mapQ value for clustering, default: 1 (Enable only when the input mapping file is in BAM format).
  --scaffold_q        <scaffold_q>               : Filter mapQ value for scaffolding, default: 0 (Enable only when the input mapping file is in BAM format).

>>> clustering chromosomes Parameters:
  --split_gfa_n       <split_gfa_n>              : Number of common neighbors when splitting GFA [2-5], default: 5.
  --chr_pm            <partig_chr_pm>            : Similarity of partig when clustering chr [0.8 <= x < 1], default: 0.95.

>>> clustering haplotypes Parameters:
  --hap_pm            <partig_hap_pm>            : Similarity of partig when clustering hap [0.6 <= x < 1], default: 0.7.
  --expand                                       : Whether to expand the allele, default: False.
  --rescue                                       : Whether to rescue the subgraph, default: False.
  --reassign_number   <reassign_number>          : Number of reassign step [1-3], default: 1.

>>> scaffolding haplotypes Parameters:
  --thread            <thread>                   : Number of parallel processes, default: 12.
  --no_contig_ec                                 : do not do contig error correction in YaHS, default: False.
  --no_scaffold_ec                               : do not do scaffold error correction in YaHS, default: False.
  --min_len           <min_len>                  : minimum scaffold length(kb) in haphic sort [0-1000], default: 200.
  --mutprob           <mutprob>                  : mutation probability [0.1-0.9] default: 0.6.
  --ngen              <ngen>                     : generations for GA, default: 20000.
  --npop              <npop>                     : population size, default: 200.
  --processes         <processes>                : processes for fast sorting, default: 32.

-h, --help            Show this help message
USAGE
    exit 1
}

# ====== defaults and flag initialization ======
fa_file=""
gfa=""
collapse_num_file=""
map_file=""
output_prefix="gphase"
n_chr=""
n_hap=""
output_prefix=""
enzyme_site="GATC"
cluster_q=1
scaffold_q=0
split_gfa_n=5
chr_pm="0.95"
hap_pm="0.7"
thread=12
min_len=200
mutprob="0.6"
ngen=20000
npop=200
processes=32
reassign_number=1
rescue_flag=""
expand_flag=""
no_contig_ec=""
no_scaffold_ec=""


# ===== parse args =====
TEMP=$(getopt -o f:g:c:m:p:e:h --long n_chr:,n_hap:,cluster_q:,scaffold_q:,chr_pm:,hap_pm:,split_gfa_n:,rescue,expand,reassign_number:,thread:,no_contig_ec,no_scaffold_ec,min_len:,mutprob:,ngen:,npop:,processes:,help -- "$@")
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
        -e) enzyme_site="$2"; shift 2 ;;
        --cluster_q) cluster_q="$2"; shift 2 ;;
        --scaffold_q) scaffold_q="$2"; shift 2 ;;
        --split_gfa_n) 
            if [[ "$2" =~ ^[0-9]+$ ]] && [ "$2" -ge 2 ] && [ "$2" -le 5 ]; then split_gfa_n="$2"; else die "--split_gfa_n must be integer 2-5"; fi
            shift 2 ;;
        --chr_pm) 
            if (($(echo "$2 >= 0.8 && $2 < 1" | bc -l) )); then chr_pm="$2"; else die "--chr_pm must be 0.8 <= x < 1"; fi
            shift 2 ;;
        --hap_pm) 
            if (($(echo "$2 >= 0.6 && $2 < 1" | bc -l) )); then hap_pm="$2"; else die "--hap_pm must be 0.6 <= x < 1"; fi
            shift 2 ;;
        --reassign_number) 
            if [[ "$2" =~ ^[1-3]$ ]]; then reassign_number="$2"; else die "--reassign_number must be 1, 2 or 3"; fi
            shift 2 ;;
        --rescue) rescue_flag="--rescue"; shift ;; 
        --expand) expand_flag="--expand"; shift ;;
        --thread) thread="$2"; shift 2 ;;
        --no_contig_ec) no_contig_ec="--no_contig_ec"; shift ;;
        --no_scaffold_ec) no_scaffold_ec="--no_scaffold_ec"; shift ;;
        --min_len) 
            if [[ "$2" =~ ^[0-9]+$ ]] && [ "$2" -le 1000 ]; then min_len="$2"; else die "--min_len must be 0-1000"; fi
            shift 2 ;;
        --mutprob) 
            if (($(echo "$2 >= 0.1 && $2 <= 0.9" | bc -l) )); then mutprob="$2"; else die "--mutprob must be 0.1-0.9"; fi
            shift 2 ;;
        --ngen) ngen="$2"; shift 2 ;;
        --npop) npop="$2"; shift 2 ;;
        --processes) processes="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) usage ;;
    esac
done

[[ -z "$fa_file" || -z "$gfa" || -z "$collapse_num_file" || -z "$map_file" || -z "$n_chr" || -z "$n_hap" || -z "$output_prefix" ]] && die "Missing required arguments. See -h/--help."

[[ "$output_prefix" =~ ^[a-zA-Z0-9.]+$ ]] || die "--output_prefix may only contain characters [a-zA-Z0-9.]"
[[ ! "$n_chr" =~ ^[0-9]+$ || "$n_chr" -le 0 ]] && die "--n_chr must be positive integer"
[[ ! "$n_hap" =~ ^[0-9]+$ || "$n_hap" -le 0 ]] && die "--n_hap must be positive integer"
[[ ! "$cluster_q" =~ ^[0-9]+$ ]] && die "--cluster_q must be non-negative integer"
[[ ! "$scaffold_q" =~ ^[0-9]+$ ]] && die "--scaffold_q must be non-negative integer"
[[ "$cluster_q" -lt "$scaffold_q" ]] && die "--cluster_q must be >= --scaffold_q"


for cmd in python samtools bc awk realpath; do 
    command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
done

case "${map_file,,}" in
    *.pairs|*.bam) ;;
    *) die "Error: Hi-C/Pore-C mapping file must be *.pairs (PA5 format) or *.bam (BAM format): $map_file" ;;
esac

[[ ! -f "$fa_file" ]] && die "Input FASTA file not found: $fa_file"
[[ ! -f "$gfa" ]] && die "Input GFA file not found: $gfa"
[[ ! -f "$collapse_num_file" ]] && die "Collapse number file not found: $collapse_num_file"
[[ ! -e "$map_file" ]] && die "Hi-C/Pore-C mapping file (bam/pairs) not found: $map_file"

fa_file="$(realpath "$fa_file")" && check_file_exists_and_nonempty "$fa_file" "input FASTA"
gfa="$(realpath "$gfa")" && check_file_exists_and_nonempty "$gfa" "input GFA"
collapse_num_file="$(realpath "$collapse_num_file")" && check_file_exists_and_nonempty "$collapse_num_file" "collapse number file"
map_file="$(realpath "$map_file")"
[[ -e "$map_file" ]] || die "Hi-C/Pore-C Mapping file not accessible after realpath: $map_file"


workdir="${current_dir}/gphase_output"
mkdir -p "${workdir}"
log_file="${current_dir}/gphase_pipeline.log"
info "Working directory: ${workdir}"
info "Log file: ${log_file}"
cd "${workdir}"
info "Changed directory to: $(pwd)"

# ============================== Preprocessing ==============================
info "=== Starting Preprocessing ==="
workdir_preprocessing="${workdir}/preprocessing"
mkdir -p preprocessing && cd preprocessing

safe_ln "$fa_file"

# 1. get_RE.py
check_file_exists_and_nonempty "$(basename "$fa_file")" "FASTA for get_RE"
run_step "python ${SCRIPT_DIR}/../cluster_chr/get_RE.py -f $(basename "$fa_file") -e ${enzyme_site} -op ${output_prefix}" "get_RE.py"
expected_RE="${output_prefix}.RE_counts.txt"
check_file_exists_and_nonempty "${expected_RE}" "RE_counts output"

# 2. Generate or link map.pairs
map_basename=$(basename "$map_file")
if [[ "$map_basename" =~ \.pairs$ ]]; then 
    if [ ! -e "map.pairs" ]; then
        ln -s "$map_file" map.pairs
        info "Linked: $map_file -> map.pairs"
    else
        info "skip link: map.pairs already exists"
    fi
elif [[ "$map_basename" =~ \.bam$ ]]; then
    safe_ln "$map_file" "$map_file"
    check_file_exists_and_nonempty "$(basename "$fa_file")" "FASTA for samtools faidx"
    
    info "Indexing FASTA for samtools"
    run_step "samtools faidx $(basename "$fa_file")" "samtools faidx"

    info "Converting BAM to map.pairs format (mapQ >= $scaffold_q)"
    tmp_pairs="map.pairs.$$"
    {
        echo "## pairs format v1.0.0"
        echo "##columns: readID chrom1 pos1 chrom2 pos2 mapQ"
        awk '{print "##chromsize: "$1" "$2}' "$(basename "$fa_file").fai"
    } > "$tmp_pairs"
    
    run_step "samtools view -@ ${thread} \"$map_file\" | awk -v q=\"${scaffold_q}\" '(\$5+0 >= q && \$7!=\"=\"){print \$1\"\t\"\$3\"\t\"\$4\"\t\"\$7\"\t\"\$8\"\t\"\$5}' >> \"$tmp_pairs\"" "BAM to pairs conversion"
    mv "$tmp_pairs" map.pairs
else
    die "Unsupported map file extension for: $map_basename. Must be .pairs or .bam"
fi

check_file_exists_and_nonempty "map.pairs" "generated map.pairs"

# 3. get_links.py
run_step "python ${SCRIPT_DIR}/../cluster_chr/get_links.py -i map.pairs -o ${output_prefix} -q ${cluster_q}" "get_links.py"
links_csv="${output_prefix}.map.links.csv"

lines=$(wc -l "${links_csv}" | awk '{print $1}')
if [[ "${lines}" -eq 1 ]]; then
    die "Error: The link file ${links_csv} must contain exactly 1 data line."
fi

# 4. nor_hic.py
run_step "python ${SCRIPT_DIR}/../cluster_chr/nor_hic.py -f ${links_csv} -r ${expected_RE} -o ${output_prefix}.map.links.nor.csv" "nor_hic.py"
normalized_links="${output_prefix}.map.links.nor.csv"
check_file_exists_and_nonempty "$normalized_links" "normalized links"

# ============================== Cluster chromosomes ==============================
info "=== Starting Chromosome Clustering ==="
cd "${workdir}"
mkdir -p cluster_chr && cd cluster_chr

safe_ln "$fa_file"
safe_ln "$gfa"
safe_ln "${workdir_preprocessing}/${expected_RE}"
safe_ln "${workdir_preprocessing}/${normalized_links}"

check_file_exists_and_nonempty "$(basename "$fa_file")" "FASTA in cluster_chr"
check_file_exists_and_nonempty "$(basename "$gfa")" "GFA in cluster_chr"
check_file_exists_and_nonempty "${output_prefix}.RE_counts.txt" "RE_counts"
check_file_exists_and_nonempty "${output_prefix}.map.links.nor.csv" "normalized links"

run_step "python ${SCRIPT_DIR}/../cluster_chr/cluster_chr.py -f $(basename "$fa_file") -r ${output_prefix}.RE_counts.txt -l ${output_prefix}.map.links.nor.csv -op ${output_prefix} -n_chr ${n_chr} -g $(basename "$gfa") -pm ${chr_pm} --split_gfa_n ${split_gfa_n}" "cluster_chr.py"

check_file_exists_and_nonempty "${output_prefix}.digraph.csv" "digraph from cluster_chr"
check_file_exists_and_nonempty "group_ctgs_All.txt" "group_ctgs_All.txt"
check_file_exists_and_nonempty "${output_prefix}.chr.cluster.ctg.txt" "clusters file from cluster_chr"
check_file_exists_and_nonempty "rescue.cluster.ctg.txt" "rescue clusters file from cluster_chr"

# ============================== Cluster haplotypes ==============================
info "=== Starting Haplotype Clustering ==="
cd "${workdir}"
mkdir -p cluster_hap && cd cluster_hap

cluster_chr_dir="${workdir}/cluster_chr"
hap_collapse_dir="${workdir}/cluster_hap"


safe_ln "$collapse_num_file"
safe_ln "$fa_file"
safe_ln "${cluster_chr_dir}/group_ctgs_All.txt"
safe_ln "${cluster_chr_dir}/rescue.cluster.ctg.txt"
safe_ln "${cluster_chr_dir}/${output_prefix}.RE_counts.txt"
safe_ln "${cluster_chr_dir}/${output_prefix}.map.links.nor.csv"
safe_ln "${cluster_chr_dir}/${output_prefix}.digraph.csv"
safe_ln "${cluster_chr_dir}/${output_prefix}.chr.cluster.ctg.txt"

cr_file="${cluster_chr_dir}/rescue.cluster.ctg.txt"

run_step "python ${SCRIPT_DIR}/../cluster_hap/cluster_hap.py -f $(basename "$fa_file") -r ${output_prefix}.RE_counts.txt -l ${output_prefix}.map.links.nor.csv -op ${output_prefix} -n_chr ${n_chr} -n_hap ${n_hap} --collapse_num_file $(basename "$collapse_num_file") -d ${output_prefix}.digraph.csv -s group_ctgs_All.txt -c ${output_prefix}.chr.cluster.ctg.txt -cr ${cr_file} -pm ${hap_pm} --reassign_number ${reassign_number} ${rescue_flag} ${expand_flag}" "cluster_hap.py"

for chr_id in $(seq 1 "$n_chr"); do
    reassign_file="${output_prefix}.reassign.cluster.txt"
    check_file_exists_and_nonempty "chr${chr_id}/$reassign_file" "reassign cluster file for chromosome ${chr_id}"
done

# ============================== Scaffold haplotypes ==============================
info "=== Starting Haplotype Scaffolding ==="
cd "${workdir}"
mkdir -p scaffold_hap && cd scaffold_hap


safe_ln "$fa_file"
safe_ln "${cluster_chr_dir}/${output_prefix}.RE_counts.txt"
safe_ln "${cluster_chr_dir}/${output_prefix}.map.links.nor.csv"
safe_ln "$gfa"
safe_ln "${cluster_chr_dir}/group_ctgs_All.txt"
safe_ln "${cluster_chr_dir}/${output_prefix}.digraph.csv"
safe_ln "${workdir_preprocessing}/map.pairs" "map.pairs"

check_file_exists_and_nonempty "map.pairs" "map.pairs for scaffolding"

run_step "python ${SCRIPT_DIR}/../scaffold_hap/scaffold_hap_v2.py -f $(basename "$fa_file") -r ${output_prefix}.RE_counts.txt -l ${output_prefix}.map.links.nor.csv -op ${output_prefix} -n_chr ${n_chr} -n_hap ${n_hap} -CHP ../cluster_hap -s group_ctgs_All.txt -g $(basename "$gfa") -d ${output_prefix}.digraph.csv -m map.pairs -t ${thread} ${no_contig_ec} ${no_scaffold_ec} --min_len ${min_len} --mutprob ${mutprob} --ngen ${ngen} --npop ${npop} --processes ${processes}" "scaffold_hap_v2.py"

# ====== Final summary ======
cd "${current_dir}"
info "==== Pipeline completed successfully! ===="
info "Key output files:"
for f in \
    "gphase_output/preprocessing/${expected_RE}" \
    "gphase_output/preprocessing/${normalized_links}" \
    "gphase_output/cluster_chr/${output_prefix}.digraph.csv" \
    "gphase_output/cluster_chr/group_ctgs_All.txt" \
    "gphase_output/scaffold_hap/${output_prefix}.scaffolds.fasta" \
    "gphase_output/scaffold_hap/${output_prefix}.scaffolds.agp"; do
    [[ -s "${f}" ]] && info "OK   ${f}" || warn "MISSING ${f}"
done

info "Full log: ${log_file}"
info "All done."

exit 0