#!/usr/bin/env python3

import os
import argparse
import subprocess
import logging
import functools
import sys
import pandas as pd
import argcomplete
import networkx as nx
from argcomplete.completers import FilesCompleter
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from find_knees import find_best_knee
from multilevel_cluster_v2 import multilevel_cluster
from louvain_reassign_allele import louvain_reassign_allele
from louvain_nei import louvain_nei



def setup_logging(log_file: str = "cluster_hap.log") -> logging.Logger:
    """Configure logging to both file and console."""
    logger = logging.getLogger('cluster_hap')
    logger.setLevel(logging.INFO)
    
    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s <%(module)s.py> [%(funcName)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def parse_arguments() -> argparse.Namespace:

    parser = argparse.ArgumentParser(prog='cluster_hap')
    # Required arguments
    base_group  = parser.add_argument_group('>>> Parameters for basic data')
    base_group.add_argument("-f", "--fa_file",metavar='\b', required=True, help="Path to the assembly FASTA file.")
    base_group.add_argument("--collapse_num_file", metavar='\b',required=True, help="Collapse number file")
    base_group.add_argument("-r", "--RE_file", metavar='\b',required=True, help="Restriction enzyme file")

    hic_group  = parser.add_argument_group('>>> Parameters for HiC data alignment')
    hic_group.add_argument("-l", "--HiC_file", metavar='\b',required=True, help="Hi-C file")


    graph_group  = parser.add_argument_group('>>> Parameters related to subgraphs')
    graph_group.add_argument("-d", "--digraph_file", metavar='\b',required=True, help="Directed graph file")
    graph_group.add_argument("-s", "--subgraph_file", metavar='\b',required=True, help="Subgraph file")

    cluster_group  = parser.add_argument_group('>>> Parameters of the file for the chr clustering results')
    cluster_group.add_argument("-c", "--chr_cluster_file", metavar='\b',required=True, help="Chromosome cluster file")
    cluster_group.add_argument("-cr", "--chr_cluster_rescue_file",metavar='\b', required=True, help="Chromosome cluster rescue file")

    partig_group  = parser.add_argument_group('>>> Parameter for partig')
    partig_group.add_argument("-pk", "--partig_k",metavar='\b', type=int, default=17, help="K-mer size for Partig. Default: 17.")
    partig_group.add_argument("-pw", "--partig_w", metavar='\b',type=int, default=17, help="Minimizer window size for Partig. Default: 17.")
    partig_group.add_argument("-pc", "--partig_c", metavar='\b',type=int, default=60, help="Max occurrance for Partig. Default: 60.")
    partig_group.add_argument("-pm", "--partig_m", metavar='\b',type=float, default=0.6, help="Mini k-mer similarity for Partig. Default: 0.6.")

    output_group  = parser.add_argument_group('>>> Parameter for the prefix of the result file')
    output_group.add_argument("-op", "--output_prefix", metavar='\b',required=True, help="Output prefix")

    genome_group  = parser.add_argument_group('>>> Parameters of chromosome and haplotype numbers')
    genome_group.add_argument("-n_chr", "--chr_number",metavar='\b', type=int, required=True, help="Number of chromosomes")
    genome_group.add_argument("-n_hap", "--hap_number", metavar='\b',type=int, required=True, help="Number of haplotypes")

    performance_group  = parser.add_argument_group('>>> Parameters for performance')
    performance_group.add_argument("--process_num", metavar='\b',type=int, default=8,
                               help=f"Number of parallel processes (default: 8)")
    performance_group.add_argument("--thread_num", metavar='\b',type=int,default=8,
                               help=f"Number of threads per process (default: 8)")
    
    reassign_group = parser.add_argument_group('>>> Parameters for reassign step')
    reassign_group.add_argument("--reassign_number", metavar='\b',type=int, default=1,
                               help=f"Number of reassign step (default: 1)")

    Optional_group  = parser.add_argument_group('>>> Optional parameters')
    Optional_group.add_argument('--expand', action='store_true', help='expand the allele')
    Optional_group.add_argument('--rescue', action='store_true', help='rescue filtered subgraph ')


    return parser.parse_args()

def execute_command(command, error_message,logger):
    try:
        # logger.info(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True, executable='/bin/bash')
        # logger.info(f"Running command successfully: {command}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(command)}")
        logger.error(f"Error: {e}")
    except Exception as e:
        logger.error(f"{error_message}")
        logger.error(f"Error executing command: {e}")
        return False

# --- Check if file exists and is not empty ---
def check_file_exists_and_not_empty(filepath: str, logger: logging.Logger, action_name: str) -> None:
    """Check if a file exists and is not empty. Raise an error if it fails."""
    # logger.info(action_name)
    if not os.path.exists(filepath):
        logger.error(f"File check failed for: {filepath}")
        raise FileNotFoundError(f"{action_name} failed: Required file not found: {filepath}")
    if os.path.getsize(filepath) == 0:
        logger.error(f"File check failed for: {filepath}")
        raise EOFError(f"{action_name} failed: File is empty: {filepath}")
    logger.info(f"File check passed for: {filepath}")
# --- End ---

def create_symlink(source, dest,logger):
    if not os.path.exists(dest):
        os.symlink(source, dest)
        logger.info(f"Created symlink: {source} -> {dest}")
    else:
        logger.warning(f"Symlink already exists: {dest}")

def filter_links_by_utgs(flag, utg_file, input_file, output_file,logger):
    temp_file = f"temp_utg_file_{flag}.txt"
    try:
        command = f"cut -f1 {utg_file} > {temp_file}"
        execute_command(command, f"Failed to create temp file from {utg_file}", logger)
        # --- Check temp file ---
        check_file_exists_and_not_empty(temp_file, logger, f"Creating temporary UTG list for {flag}")
        # --- End  ---
        
        command = (
            f"awk -F \"[,\\t ]\" 'NR==FNR{{lines[$1]; next}} ($1 in lines && $2 in lines)' "
            f"{temp_file} {input_file} > {output_file}"
        )
        execute_command(command, f"Failed to filter links for {output_file}", logger)

        # --- Check output file ---
        check_file_exists_and_not_empty(output_file, logger, f"Filtering links for {flag}")
        # --- End ---
        
        command = f"rm {temp_file}"
        execute_command(command, f"Failed to remove temp file {temp_file}", logger)
        
        command = f"sed -i '1isource,target,links' {output_file}"
        execute_command(command, f"Failed to add header to {output_file}", logger)
    
    except Exception as e:
        logger.error(f"An error occurred in filter_links_by_utgs: {e}")
        raise

def adjust_r_and_cluster(initial_r, min_r, max_r, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger):

    r = initial_r
    check_cluster_dict = defaultdict(list)
    while min_r <= r <= max_r:
        # logger.info(f"Running clustering with r={r}")
        check_cluster_dict, max_group_allele_value = multilevel_cluster(csv_file, cluster_output, r , "check", utg_file, partig_file, hap_number)
        num_clusters = len(check_cluster_dict)
        # logger.info(f"Number of clusters: {num_clusters} (Target: {hap_number})")

        if num_clusters == hap_number:
            # logger.info(f"Optimal r found: {r}")
            return r 
        elif num_clusters > hap_number:
            max_r = r - step
        else:
            min_r = r + step
        r = (min_r + max_r) / 2

    return None

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len


def read_collapse_num(collapse_num_file):

    collapse_num_dict = defaultdict()
    with open(collapse_num_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0][1] == "#":
                continue
            try:
                collapse_num_dict[line[0]] = int(line[1])
            except:
                collapse_num_dict[line[0]] = 1
    return collapse_num_dict

def get_avg_uncollapse_num(REFile, collapse_num_file, hap_number):

    ctg_RE_len = read_REs(REFile)
    collapse_num_dict = read_collapse_num(collapse_num_file)

    uncollapse_avg = sum([ list_[1] for ctg, list_ in ctg_RE_len.items() if collapse_num_dict[ctg] < 2]) / int(hap_number)

    return uncollapse_avg



def multiple_adjust_r_and_cluster(initial_r, min_r, max_r, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger):


    section_length = float((max_r - min_r) / 8)

    section_1 = adjust_r_and_cluster((min_r + min_r+section_length)/2, min_r, min_r+section_length, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_2 = adjust_r_and_cluster((min_r+section_length + min_r+section_length*2)/2, min_r+section_length, min_r+section_length*2, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_3 = adjust_r_and_cluster((min_r+section_length*2 + min_r+section_length*3)/2, min_r+section_length*2, min_r+section_length*3, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_4 = adjust_r_and_cluster((min_r+section_length*3 + min_r+section_length*4)/2, min_r+section_length*3, min_r+section_length*4, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)

    section_5 = adjust_r_and_cluster((min_r+section_length*4 + min_r+section_length*5)/2, min_r+section_length*4, min_r+section_length*5, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_6 = adjust_r_and_cluster((min_r+section_length*5 + min_r+section_length*6)/2, min_r+section_length*5, min_r+section_length*6, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_7 = adjust_r_and_cluster((min_r+section_length*6 + min_r+section_length*7)/2, min_r+section_length*6, min_r+section_length*7, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_8 = adjust_r_and_cluster((min_r+section_length*7 + max_r)/2, min_r+section_length*7, max_r, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)


    if section_1 or section_2 or section_3 or section_4 or section_5 or section_6 or section_7 or section_8:
        r_list = [section_1, section_2, section_3, section_4, section_5, section_6, section_7, section_8]
        r_filter_list = [ r for r in r_list if r and r > 0]
        r = min(r_filter_list, key=lambda x: abs(x - 1))
        return r
    raise


def filter_edges_by_density(chr_num, HiC_file, utg_rescue_file, filter_HiC_file, logger, step=0.5):
    # --- Check input files ---
    check_file_exists_and_not_empty(utg_rescue_file, logger, f"Reading UTGs for Chr{chr_num}")
    check_file_exists_and_not_empty(HiC_file, logger, f"Reading Hi-C links for Chr{chr_num}")

    nodes = pd.read_csv(utg_rescue_file, header=None)[0].tolist()
    edges = pd.read_csv(HiC_file, sep=',', header=0)
    edges.columns = ['source', 'target', 'links']

    num_nodes = len(nodes)
    num_edges = len(edges)
    density = (2 * num_edges) / (num_nodes * (num_nodes - 1)) if num_nodes > 1 else 0
    logger.info(f"Chr{chr_num} hic graph info : edges -> {num_edges}\t nodes -> {num_nodes}\t density -> {density}")

    if  num_nodes < (num_edges / 50):
        logger.info(f"Chr{chr_num} HiC signal filtering...")
        threshold = 0.5
        
        while True:
            filtered_edges = edges[edges['links'] > threshold]
            filtered_num_edges = len(filtered_edges)

            logger.info(f"Chr{chr_num} HiC signal filtering : threshold -> {threshold:.1f}\t edges: -> {filtered_num_edges}")

            if filtered_num_edges < (num_nodes * 50):
                break

            threshold += step

        if filtered_edges.empty:
             logger.warning(f"Chr{chr_num} HiC signal filtering resulted in an empty edge set. Using all edges instead.")
             edges.to_csv(filter_HiC_file, sep=',', header=False, index=False)
             return 0
        else:
             filtered_edges.to_csv(filter_HiC_file, sep=',', header=False, index=False)
             # --- Check output file ---
             check_file_exists_and_not_empty(filter_HiC_file, logger, f"Filtering edges by density for Chr{chr_num}")
             return threshold
    else:
        edges.to_csv(filter_HiC_file, sep=',', header=False, index=False)
        # --- Check output file ---
        check_file_exists_and_not_empty(filter_HiC_file, logger, f"Copying all edges for Chr{chr_num}")
        return 0


def run_spectral_clustering_fallback(input_hic_file: str, final_cluster_output: str, hap_number: int, logger: logging.Logger) -> bool:
    """
    Runs Spectral Clustering on the HiC link file as a fallback when multilevel clustering fails.
    """
    logger.info("Starting: Spectral Clustering Fallback for haplotype Clustering.")

    try:
        try:
            edges_df = pd.read_csv(input_hic_file, sep=',', header=0)
            edges_df.columns = ['source', 'target', 'links']
            try:
                all_nodes = list(pd.concat([edges_df.iloc[:, 0], edges_df.iloc[:, 1]]).unique())
                num_nodes = len(all_nodes)
                logger.info(f"Loaded {num_nodes} contigs from {input_hic_file}.")
            except Exception as e:
                logger.error(f"Failed to load nodes from {input_hic_file}: {e}")
                return False
            
            G = nx.from_pandas_edgelist(
                edges_df,
                source='source',
                target='target',
                edge_attr='links',
                create_using=nx.Graph()
            )
            A_matrix = nx.to_numpy_array(G, nodelist=all_nodes, weight='links')
            logger.info(f"Built adjacency matrix A with shape {A_matrix.shape}.")
            
        except Exception as e:
            logger.error(f"Failed to build adjacency matrix from {input_hic_file}: {e}")
            return False

        logger.info(f"Running SpectralClustering with k={hap_number}.")

        from sklearn.cluster import SpectralClustering
        sc = SpectralClustering(
            n_clusters=hap_number,
            affinity='precomputed',
            assign_labels='kmeans',
            random_state=42,
            n_init=20
        )
        clusters = sc.fit_predict(A_matrix)

        results = pd.DataFrame({
            'Node': all_nodes,
            'Cluster_ID': clusters
        })
        summary_dict = results.groupby('Cluster_ID')['Node'].apply(list).to_dict()

        group_counter = 0
        with open(final_cluster_output, 'w') as file:
            for idx, utgs in summary_dict.items():
                if not utgs:
                    continue
                group_counter += 1
                file.write(f"group{group_counter}\t{len(utgs)}\t{' '.join(utgs)}\n")
        
        logger.info(f"Completed Spectral Clustering. Results written to {final_cluster_output}.")
        final_line_count = len(results)
        if final_line_count == num_nodes:
            return True
        else:
            logger.warning(f"Spectral Clustering output line count ({final_line_count}) does not match node count ({num_nodes}).")
            return False

    except Exception as e:
        logger.error(f"An unexpected error occurred during Spectral Clustering: {e}")
        return False


def process_chromosome(chr_num, args, pwd, partig_file,logger):
    script_path = os.path.abspath(sys.path[0])
    try:

        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        os.makedirs(chr_dir, exist_ok=True)

        os.chdir(chr_dir)
        logger.info(f"Processing chromosome {chr_num} in {chr_dir}")

        # Create symlinks with ThreadPoolExecutor
        param_files = {
            "chr_cluster_file": f"../{args.chr_cluster_file}",
            "digraph_file": f"../{args.digraph_file}",
            "subgraph_file": f"../{args.subgraph_file}",
            "collapse_num_file": f"../{args.collapse_num_file}",
            "HiC_file": f"../{args.HiC_file}",
            "RE_file": f"../{args.RE_file}",
            "partig_file": f"../{partig_file}"
        }
        # "merge_partig_file": f"../merge.partig.csv"
        no_expand_partig_file = partig_file
        if args.expand:
            origin_partig_file = "merge.partig.csv"
            create_symlink(f"../merge.partig.csv", origin_partig_file, logger)
        else:
            origin_partig_file = partig_file
        
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            futures = []
            for name, filepath in param_files.items():
                dest = os.path.join(chr_dir, os.path.basename(filepath))
                futures.append(executor.submit(create_symlink, filepath, dest, logger))
            
            for future in futures:
                future.result()


        utg_file = f"{args.output_prefix}.chr{chr_num}.utgs.txt"
        create_symlink(f"../group{chr_num}.txt", utg_file, logger)
        # --- Check utg_file ---
        check_file_exists_and_not_empty(utg_file, logger, f"Checking UTG list for Chr{chr_num}")

        if args.rescue:
            utg_rescue_file = f"{args.output_prefix}.chr{chr_num}.utgs.rescue.txt"
            create_symlink(f"../group{chr_num}_rescue.txt", utg_rescue_file, logger)
        else:
            utg_rescue_file = utg_file
        
        check_file_exists_and_not_empty(utg_rescue_file, logger, f"Checking UTG rescue list for Chr{chr_num}")

        
        # check utg number
        try:
            df = pd.read_csv(utg_rescue_file, sep='\t', skipinitialspace=True)
            df_len = len(df)
        except Exception as e:
            raise ValueError(f"Chr{chr_num}: Loaded utg file error...")

        if df_len < int(args.hap_number):
            logger.error(
            f"Chr{chr_num}: The number of Unitigs loaded from '{utg_rescue_file}' ({df_len}) "
            f"is less than the required haplotype number ({int(args.hap_number)}). "
            f"Cannot proceed with clustering. Please check if the chromosome clustering results are incorrect.")
            raise ValueError(f"The number of Unitigs on chromosome {chr_num} is less than the number of haplotypes...") 

        # Process files with ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            # Process links file
            links_file = f"{args.output_prefix}.chr{chr_num}.links.all.nor.csv"
            links_future = executor.submit(filter_links_by_utgs,str(chr_num)+"_links", utg_rescue_file, args.HiC_file, links_file, logger)
            
            # Process partig file
            partig_file_chr = f"{args.output_prefix}.chr{chr_num}.partig.csv"
            partig_future = executor.submit(filter_links_by_utgs,str(chr_num)+"_partig", utg_rescue_file, origin_partig_file, partig_file_chr, logger)
            
            links_future.result()
            partig_future.result()

        # filter HiC links
        filter_edge_file = f"{args.output_prefix}.chr{chr_num}.links.nor.csv"
        cut_value = filter_edges_by_density(chr_num, links_file, utg_rescue_file, filter_edge_file, logger=logger, step=0.5)

        filtered_links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.filterAllele.csv"
        # --- Check partig_file_chr ---
        check_file_exists_and_not_empty(partig_file_chr, logger, f"Checking Partig file for filtering allele links for Chr{chr_num}")
        # --- End ---

        command = (
            f"awk -F '[, \\t]' 'NR==FNR{{lines[$1\",\"$2];next}}!($1\",\"$2 in lines)' "
            f"{partig_file_chr} {filter_edge_file} > {filtered_links_file}"
        )
        execute_command(command, f"Failed to filter allele links for {filtered_links_file}",logger)
        check_file_exists_and_not_empty(filtered_links_file, logger, f"Checking filtered links file for Chr{chr_num}")
        execute_command(f"sed '1isource,target,links' -i {args.HiC_file}", f"Failed to add header to {args.HiC_file}",logger)
        execute_command(f"sed '1isource,target,links' -i {filtered_links_file}", f"Failed to add header to {filtered_links_file}",logger)


        # # The knee is used to filter HiC signals
        best_knee = find_best_knee(filtered_links_file, f"{args.output_prefix}.chr{chr_num}.knees")
        cut_value_step = 1

        if best_knee < 5:
            cut_value_max = 5
        elif best_knee > 15:
            cut_value_max = 15
        else:
            cut_value_max = best_knee

        cluster_output = f"{args.output_prefix}.chr{chr_num}.cluster.txt"

        avg_uncollapse_num = get_avg_uncollapse_num(utg_rescue_file, args.collapse_num_file, args.hap_number)
        logger.info(f"chr_num:{chr_num}\tbest_knee:{best_knee}")

        # tuple ： (check_cluster_dict, max_group_allele_value）
        satisfy_clusters_list = list()
        check_cluster_dict = defaultdict(list)

        while cut_value < cut_value_max:

            cut_links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.filterAllele.c{float(cut_value)}.csv"
            command = (
                f"awk -F ',' '($3> {float(cut_value)})' "
                f"{filtered_links_file} > {cut_links_file}"
            )
            execute_command(command, f"Failed to filter hic links for {cut_links_file}",logger)

            try:
                check_file_exists_and_not_empty(cut_links_file, logger, f"Checking cut links file (cut={cut_value}) for Chr{chr_num}")
            except (FileNotFoundError, EOFError):
                cut_value += cut_value_step
                continue 

            try:
                # # Run louvain_nei.py
                louvain_nei_result = louvain_nei(args.collapse_num_file, utg_rescue_file, cut_links_file, partig_file_chr)
                if not louvain_nei_result:
                    logger.error(f"Chr{chr_num}: Local clustering -> louvain_nei run error, try to running Direct clustering...")
                    raise Exception("Louvain_nei failed")
                check_file_exists_and_not_empty("louvain_nei.csv", logger, f"Checking Louvain_nei output for Chr{chr_num}")
                check_cluster_dict, max_group_allele_value = multilevel_cluster("louvain_nei.csv", cluster_output, float(1), "check", utg_rescue_file, no_expand_partig_file, int(args.hap_number))
                optimal_r = multiple_adjust_r_and_cluster(
                    initial_r=1.0,
                    min_r=0.01,
                    max_r=3,
                    step=0.03,
                    cluster_output=cluster_output,
                    csv_file="louvain_nei.csv", 
                    utg_file=utg_rescue_file,
                    partig_file=partig_file_chr, 
                    hap_number=args.hap_number, logger=logger
                )
                if optimal_r:
                    check_cluster_dict, max_group_allele_value = multilevel_cluster("louvain_nei.csv", cluster_output, float(optimal_r), "check", utg_rescue_file, no_expand_partig_file, int(args.hap_number))
                else:
                    raise ValueError("Optimal r not found for Louvain clustering.")

                # Check the number of clusters
                num_clusters = len(check_cluster_dict)
                if num_clusters != args.hap_number:
                    raise
                else:
                    satisfy_clusters_list.append((check_cluster_dict, max_group_allele_value))

                if max_group_allele_value > 1000000:
                    logger.info(f"Chr{chr_num}: Local clustering -> uncollapse clusters max_group_allele_value: {max_group_allele_value}, discard... ")
                    raise
                else:
                    logger.info(f"Chr{chr_num}: Local clustering -> uncollapse clusters max_group_allele_value: {max_group_allele_value}, try to reassgin collaspe utg... ")
                    break

            except:
                chr_num_collapse_num_file = f"{args.output_prefix}.chr{chr_num}.utgs.uncollapse.txt"
                command = (
                    "awk -F '[, \\t]' 'NR==FNR{lines[$1]=$2;next}{if(lines[$1]<=1){print $0}}' "
                    f"{args.collapse_num_file} {utg_rescue_file} > {chr_num_collapse_num_file}"
                )
                execute_command(command, f"Failed to filter collapse Contig for {chr_num_collapse_num_file}",logger)
                # --- Check chr_num_collapse_num_file ---
                check_file_exists_and_not_empty(chr_num_collapse_num_file, logger, f"Checking uncollapsed UTG list for Chr{chr_num}")


                chr_num_uncollapse_hic_file = f"{args.output_prefix}.chr{chr_num}.links.uncollapse.csv"
                command = (
                    "awk -F '[, \\t]' 'NR==FNR{lines[$1];next}($1 in lines && $2 in lines)' "
                    f"{chr_num_collapse_num_file} {cut_links_file} > {chr_num_uncollapse_hic_file}"
                )
                execute_command(command, f"Failed to filter hic for {chr_num_uncollapse_hic_file}",logger)
                execute_command(f"sed '1isource,target,links' -i {chr_num_uncollapse_hic_file}", f"Failed to add header to {chr_num_uncollapse_hic_file}",logger)
                # --- Check chr_num_uncollapse_hic_file ---
                check_file_exists_and_not_empty(chr_num_uncollapse_hic_file, logger, f"Checking uncollapsed Hi-C links for Chr{chr_num}")
                
                chr_num_uncollapse_hic_cut_file = f"{args.output_prefix}.chr{chr_num}.links.uncollapse.c{float(cut_value)}.csv"
                command = (
                    f"awk -F ',' '($3> {float(cut_value)})' "
                    f"{chr_num_uncollapse_hic_file} > {chr_num_uncollapse_hic_cut_file}"
                )
                execute_command(command, f"Failed to filter hic links for {chr_num_uncollapse_hic_cut_file}",logger)
                # --- Check chr_num_uncollapse_hic_cut_file ---
                try:
                    check_file_exists_and_not_empty(chr_num_uncollapse_hic_cut_file, logger, f"Checking uncollapsed cut links file (cut={cut_value}) for Chr{chr_num}")
                except (FileNotFoundError, EOFError):
                    cut_value += cut_value_step
                    continue 

                # Adjust r and run multilevel_cluster.py
                try:
                    check_cluster_dict, max_group_allele_value = multilevel_cluster(chr_num_uncollapse_hic_cut_file, cluster_output, 1, "check", utg_rescue_file, no_expand_partig_file, args.hap_number)
                    optimal_r = multiple_adjust_r_and_cluster(
                        initial_r=1.0,
                        min_r=0.01,
                        max_r=3,
                        step=0.03,
                        cluster_output=cluster_output,
                        csv_file=chr_num_uncollapse_hic_cut_file, 
                        utg_file=utg_rescue_file,
                        partig_file=partig_file_chr, 
                        hap_number=args.hap_number, logger=logger
                    )
                    if optimal_r:
                        check_cluster_dict, max_group_allele_value = multilevel_cluster(chr_num_uncollapse_hic_cut_file, cluster_output, optimal_r, "check", utg_rescue_file, no_expand_partig_file, args.hap_number)
                    else:
                        raise ValueError("Optimal r not found for Direct clustering.")

                    # Check the number of clusters
                    num_clusters = len(check_cluster_dict)

                    if num_clusters != args.hap_number:
                        raise
                    else:
                        satisfy_clusters_list.append((check_cluster_dict, max_group_allele_value))

                    # if max_group_allele_value > min([avg_uncollapse_num/20, 1000000]):
                    if max_group_allele_value > 1000000:
                        logger.info(f"Chr{chr_num}: Direct clustering -> uncollapse clusters max_group_allele_value: {max_group_allele_value}, discard... ")
                        raise
                    else:
                        logger.info(f"Chr{chr_num}: Direct clustering -> uncollapse clusters max_group_allele_value: {max_group_allele_value}, try to reassgin collaspe utg... ")
                        break
                except:
                    cut_value += cut_value_step
                    continue

        # If the number of clusters is not satisfied
        if cut_value >= cut_value_max:
            if satisfy_clusters_list:  
                sorted_satisfy_clusters_list = sorted(satisfy_clusters_list, key=lambda x: x[1])
                logger.info(f"Chr:{chr_num}: the detection is not met, take the result with the smallest detection value ({sorted_satisfy_clusters_list[0][1]}) as the clustering result.")
                final_clusters = sorted_satisfy_clusters_list[0][0]

                with open(cluster_output, 'w') as file:
                    for idx, group in enumerate(final_clusters):
                        file.write(f"group{idx+1}\t{len(final_clusters[group])}\t{' '.join(final_clusters[group])}\n")
            else:
                # Check the number of clusters
                try:
                    df = pd.read_csv(cluster_output, sep='\t', skipinitialspace=True)
                    df_len = len(df)
                except Exception as e:
                    logger.error(f"Chr{chr_num}: Check the number of clusters: Loaded utg file error(cluster_output)...")
                    df_len = 0

                if df_len != int(args.hap_number):
                    logger.info(
                    f"Chr{chr_num}: The number of clusters loaded from '{cluster_output}' ({df_len}) "
                    f"is not euqal the required haplotype number ({int(args.hap_number)}). "
                    f"Spectral Clustering will be used to cluster the haplotype. Alternatively, You can adjust the --hap_pm parameter and re-cluster the data...")

                    # Spectral Clustering
                    run_spectral_clustering_fallback(filtered_links_file, cluster_output, int(args.hap_number), logger)

        cluster_file = cluster_output

        # Run louvain_reassign_allele.py
        find_best_isolated, isolated_threshold = True, 0
        min_variance_idx = louvain_reassign_allele(args.collapse_num_file, utg_rescue_file, filtered_links_file, cluster_file , args.RE_file, partig_file_chr, args.output_prefix, find_best_isolated, no_expand_partig_file, isolated_threshold)
        logger.info(f"Chr{chr_num}: louvain_reassign_allele find best isolated : {min_variance_idx}")
        
        # --- Check reassign.cluster.txt after initial run ---
        reassign_output_file = f"{args.output_prefix}.reassign.cluster.txt"
        check_file_exists_and_not_empty(reassign_output_file, logger, f"Checking reassign cluster file after initial run for Chr{chr_num}")

        for i in range(int(args.reassign_number)-1):
            min_variance_idx = louvain_reassign_allele(args.collapse_num_file, utg_rescue_file, filtered_links_file, reassign_output_file, args.RE_file, partig_file_chr, args.output_prefix, find_best_isolated, no_expand_partig_file, isolated_threshold)
            # --- Check reassign.cluster.txt after subsequent runs ---
            check_file_exists_and_not_empty(reassign_output_file, logger, f"Checking reassign cluster file after iteration {i+2} for Chr{chr_num}")

        return f"Successfully processed chromosome {chr_num}"
    except Exception as e:
        logger.error(f"Error processing chromosome {chr_num}: {e}")
        return f"Failed to process chromosome {chr_num}: {e}"
    finally:
        os.chdir(pwd)

def run_partig(fa_file, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
    # --- Check input FASTA file ---
    check_file_exists_and_not_empty(fa_file, logger, "Checking input FASTA file")
    execute_command(f"samtools faidx {fa_file}", "Failed to indexing the assembly FASTA",logger)
    # --- Check index file ---
    check_file_exists_and_not_empty(f"{fa_file}.fai", logger, "Checking FASTA index file")

    output_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt"
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "../bin/partig")
    try:
        with open(output_file, "w") as outfile:
            command = [
                script_path_add, fa_file,
                "-k", str(partig_k),
                "-w", str(partig_w),
                "-c", str(partig_c),
                "-m", str(partig_m)
            ]
            logger.info(f"Starting: Running Partig")
            logger.info(f"Command: {' '.join(command)}")
            subprocess.run(command, check=True, stdout=outfile, stderr=subprocess.PIPE, text=True)
            logger.info(f"Completed: Running Partig\n")
        # --- Check Partig output file ---
        check_file_exists_and_not_empty(output_file, logger, "Checking Partig raw output file")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in: Running Partig\nCommand failed: {' '.join(command)}")
        logger.error(f"Error output: {e.stderr}")
        return False


def convert_partig_output(fa_file, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
    partig_txt_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt"
    trans_partig_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.csv"
    
    # --- Check Partig raw file ---
    check_file_exists_and_not_empty(partig_txt_file, logger, "Checking Partig TXT input for conversion")
    
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "trans.partig.py")
    execute_command(
        f"python {script_path_add} -fai {fa_file}.fai -p {partig_txt_file} -r {output_prefix}.RE_counts.txt -o {trans_partig_file} -d {output_prefix}.digraph.csv", 
        "Converting Partig output to CSV",logger)
    
    # --- Check converted Partig CSV and digraph file ---
    check_file_exists_and_not_empty(trans_partig_file, logger, "Checking converted Partig CSV file")
    check_file_exists_and_not_empty(f"{output_prefix}.digraph.csv", logger, "Checking digraph output file")


def log_start(logger, script_name: str, version: str, args: argparse.Namespace):
    """Log the start of the program."""
    logger.info(f"Program started, {script_name} version: {version}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")



def main():
    args = parse_arguments()
    logger = setup_logging('cluster_hap.log')

    import multiprocessing
    import os
    pwd = os.getcwd()
    total_cores = multiprocessing.cpu_count()
    max_safe_total = int(total_cores * 0.9)
    if args.process_num is None or args.process_num <= 0:
        args.process_num = min(20, total_cores)
    if args.thread_num is None or args.thread_num <= 0:
        args.thread_num = 8

    current_total = args.process_num * args.thread_num
    if current_total > max_safe_total:
        old_process = args.process_num
        old_thread = args.thread_num
        old_total = current_total
        args.thread_num = max(1, max_safe_total // args.process_num)
        new_total = args.process_num * args.thread_num

        logger.warning(
            f"User parallel setting ({old_process} processes × {old_thread} threads) "
            f"would use {old_total} threads, exceeding {total_cores} cores. "
            f"Auto adjusted to {args.process_num} processes × {args.thread_num} threads "
            f"(total {new_total} threads, within safe limit)."
        )
    else:
        logger.info(
            f"Parallel setting accepted: {args.process_num} processes × {args.thread_num} threads "
            f"(total {current_total}/{total_cores} cores)"
        )

    os.environ["OMP_NUM_THREADS"] = str(args.thread_num)
    os.environ["MKL_NUM_THREADS"] = str(args.thread_num)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.thread_num)
    os.environ["NUMEXPR_NUM_THREADS"] = str(args.thread_num)

    logger.info(f"Final parallel config: {args.process_num} processes × {args.thread_num} threads/process")
    log_start(logger, "cluster_hap.py", "1.0.0", args)

    # # # Step 1: Run partig
    partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{args.partig_m}.txt"
    trans_partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{args.partig_m}.csv"
    if not os.path.isfile(partig_file):
        if not run_partig(args.fa_file, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix,logger):
             raise Exception("Partig run failed.")
    if not os.path.isfile(trans_partig_file):
        convert_partig_output(args.fa_file, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix,logger)
    
    # --- Check final Partig CSV file, as input for next step ---
    check_file_exists_and_not_empty(trans_partig_file, logger, "Checking final Partig CSV file before next step")

    # Step 2: Run cluster2group.py
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "cluster2group.py")
    
    # --- Check input files ---
    check_file_exists_and_not_empty(args.chr_cluster_file, logger, "Checking chr_cluster_file input")
    check_file_exists_and_not_empty(args.RE_file, logger, "Checking RE_file input")

    execute_command(
        f"python {script_path_add} -c {args.chr_cluster_file} -r {args.RE_file}",
        "Failed to run cluster2group.py",logger
    )
    # The output of cluster2group.py (group*.txt) is checked implicitly in process_chromosome

    if args.rescue:
        # --- Check input files ---
        check_file_exists_and_not_empty(args.chr_cluster_rescue_file, logger, "Checking chr_cluster_rescue_file input")
        execute_command(
            f"python {script_path_add} -c {args.chr_cluster_rescue_file} -r {args.RE_file} -m rescue",
            "Failed to run cluster2group.py",logger
        )

    script_path_add = os.path.join(script_path, "filter_expand_partig.py")
    # --- Check input files ---
    check_file_exists_and_not_empty(args.digraph_file, logger, "Checking digraph_file input")
    check_file_exists_and_not_empty(args.subgraph_file, logger, "Checking subgraph_file input")

    execute_command(
            f"python {script_path_add} -d {args.digraph_file} "
            f"-r {args.RE_file} -s {args.subgraph_file} -p {trans_partig_file}",
            "Failed to run filter_expand_partig.py",logger
    )
    # The output of filter_expand_partig.py (merge.partig.csv) is checked later in process_chromosome if args.expand is True.

    # Step 3: Process chromosomes in parallel
    with Pool(processes=args.process_num) as pool:

        # Create partial function with fixed arguments
        process_chr = functools.partial(process_chromosome, args=args, pwd=pwd, partig_file=trans_partig_file,logger=logger)
        
        # # Process chromosomes in parallel
        # results = pool.map(process_chr, range(1, args.chr_number + 1))
        # # Log results
        # for result in results:
        #     logger.info(result)
        #     if result.startswith("Failed"):
        #          raise Exception(f"Chromosome processing failed: {result}")

        try:
            for result in pool.imap_unordered(process_chr, range(1, args.chr_number + 1)):
                logger.info(result)
                
                if result.startswith("Failed"):
                    logger.error(f"Chromosome processing failed: {result}")
                    raise Exception(f"Chromosome processing failed: {result}")
        except Exception as e:
            logger.error(f"Aborting parallel processing due to error: {e}")
            pool.terminate()
            pool.join()
            raise
        else:
            logger.info("All chromosomes processed successfully")

    # Step 4: Merge .reassign.cluster.txt files with modified first column format
    merged_output = f"{args.output_prefix}.hap.cluster.txt"
    missing_files = []
    merged_content = []

    for chr_num in range(1, args.chr_number + 1):
        reassign_file = os.path.join(pwd, f"chr{chr_num}", f"{args.output_prefix}.reassign.cluster.txt")
        
        # --- Check reassign_file exists and is not empty ---
        try:
             check_file_exists_and_not_empty(reassign_file, logger, f"Checking final reassign file for Chr{chr_num} before merging")
        except (FileNotFoundError, EOFError) as e:
             logger.error(f"Missing or empty reassign file for Chr{chr_num}: {e}")
             missing_files.append(reassign_file)
             continue
        
        try:
            with open(reassign_file, 'r') as f:
                lines = f.readlines()
                modified_lines = []
                for line in lines:
                    if line.strip(): 
                        # Try tab-separated first, then comma (original code uses '\t' for writing final cluster_output)
                        parts = line.strip().split('\t') 
                        if len(parts) < 3: 
                             parts = line.strip().split(',')
                             if len(parts) < 3:
                                  logger.warning(f"Invalid line format in {reassign_file}: {line.strip()}")
                                  continue

                        group_num_raw = parts[0].strip()
                        group_num = group_num_raw.replace('group', '') if group_num_raw.startswith('group') else group_num_raw
                        parts[0] = f"Chr{chr_num}g{group_num}"
                        
                        modified_line = '\t'.join(parts)
                        modified_lines.append(modified_line)

                merged_content.extend(modified_lines)
        except Exception as e:
            logger.error(f"Failed to process {reassign_file} during merging: {e}")
            missing_files.append(reassign_file)

    if missing_files:
        logger.error(f"Missing or failed reassign.cluster.txt files: {missing_files}")
        raise FileNotFoundError(f"Missing or failed haplotypes phaseing files: {missing_files}")

    # Write merged content to output file
    if not merged_content:
        logger.error(f"Merged content is empty. Check input files.")
        raise EOFError("Merged output content is empty.")

    with open(merged_output, 'w') as f_out:
        f_out.write('\n'.join(merged_content) + '\n')
    
    # --- Check final merged file ---
    check_file_exists_and_not_empty(merged_output, logger, "Checking final merged haplotype cluster file")
    
    logger.info(f"Successfully merged haplotypes phaseing files into {merged_output}")
    logger.info(f"Haplotype clustering was successful; the next step is to scaffolding the haplotypes.")



if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # Final error log is critical
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)