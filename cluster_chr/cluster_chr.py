#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import subprocess
import argcomplete
import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csgraph
from sklearn.cluster import KMeans
from collections import defaultdict
from sklearn.cluster import SpectralClustering
from pipeline_allele import Pipeline_allele
from multilevel_cluster_v2 import Multilevel_cluster



def setup_logging(log_file: str = "cluster_chr.log") -> logging.Logger:
    logger = logging.getLogger('cluster_hap')
    logger.setLevel(logging.INFO)

    if not logger.handlers:

        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        
        formatter = logging.Formatter('%(asctime)s <%(module)s.py> [%(funcName)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        
        logger.addHandler(fh)
        logger.addHandler(ch)
    
    return logger

def log_start(logger: logging.Logger, script_name: str, version: str, args: argparse.Namespace):
    """Log the start of the program."""
    logger.info(f"Program started, {script_name} version: {version}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")

def check_file_exists_and_not_empty(file_path: str, logger: logging.Logger) -> bool:
    """Check if a file exists and has a size greater than zero."""
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    if os.path.getsize(file_path) == 0:
        logger.error(f"File is empty: {file_path}")
        return False
    logger.info(f"File check passed: {file_path}")
    return True

def run_command(command: list, description: str, logger: logging.Logger, output_file: str = None) -> bool:
    try:
        logger.info(f"Starting: {description}")
        logger.info(f"Command: {' '.join(command)}")
        
        stdout_target = open(output_file, "w") if output_file else subprocess.PIPE
        
        result = subprocess.run(
            command, 
            check=True, 
            stdout=stdout_target, 
            stderr=subprocess.PIPE, 
            text=True,
            timeout=3600
        )
        
        if output_file:
            stdout_target.close()
        
        if result.stderr:
             logger.warning(f"STDERR output for {description}:\n{result.stderr.strip()}")

        logger.info(f"Completed: {description}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in: {description}")
        logger.error(f"Command failed with return code {e.returncode}: {' '.join(command)}")
        logger.error(f"STDOUT:\n{e.stdout.strip()}")
        logger.error(f"STDERR:\n{e.stderr.strip()}")
        return False
    except FileNotFoundError:
        logger.error(f"Program not found: {command[0]}. Please check your PATH or installation.")
        return False
    except TimeoutError:
        logger.error(f"Timeout occurred for: {description}")
        return False
    except Exception as e:
        logger.error(f"An unexpected error occurred during {description}: {e}")
        return False


def index_fasta(fa_file: str, logger: logging.Logger) -> bool:
    """Index the assembly FASTA file using samtools faidx."""
    faidx_file = f"{fa_file}.fai"
    if check_file_exists_and_not_empty(faidx_file, logger):
        logger.info(f"FASTA index already exists: {faidx_file}. Skipping indexing.")
        return True
    
    return run_command(["samtools", "faidx", fa_file], "Indexing the assembly FASTA", logger)

def split_gfa(gfa: str, split_gfa_n: int, split_gfa_iter: int, output_prefix: str, logger: logging.Logger) -> bool:
    """Split the GFA file using split_GFA.py."""
    output_file = f"{output_prefix}.rmTip.split.gfa"
    if check_file_exists_and_not_empty(output_file, logger):
        logger.info(f"Split GFA output already exists: {output_file}. Skipping split.")
        return True
    
    script_path = os.path.abspath(sys.path[0])
    command = [
        "python", os.path.join(script_path, "split_GFA.py"),
        "-g", gfa,
        "-n", str(split_gfa_n),
        "-iter", str(split_gfa_iter),
        "-o", output_prefix
    ]
    return run_command(command, "Splitting the GFA file", logger)

def run_partig(fa_file: str, partig_k: int, partig_w: int, partig_c: int, partig_m: float, output_prefix: str, logger: logging.Logger) -> bool:
    """Run the Partig tool."""
    output_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m:.2f}.txt"
    if check_file_exists_and_not_empty(output_file, logger):
        logger.info(f"Partig output already exists: {output_file}. Skipping Partig run.")
        return True
    
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "../bin/partig")
    
    command = [
        script_path_add, fa_file,
        "-k", str(partig_k),
        "-w", str(partig_w),
        "-c", str(partig_c),
        "-m", str(partig_m)
    ]
    return run_command(command, "Running Partig", logger, output_file=output_file)


def convert_partig_output(fa_file: str, partig_k: int, partig_w: int, partig_c: int, partig_m: float, output_prefix: str, RE_file: str, logger: logging.Logger) -> bool:
    """Convert Partig output to CSV format."""
    input_partig_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m:.2f}.txt"
    output_csv_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m:.2f}.csv"
    
    if not check_file_exists_and_not_empty(input_partig_file, logger):
        return False
    
    if check_file_exists_and_not_empty(output_csv_file, logger):
        logger.info(f"Partig CSV output already exists: {output_csv_file}. Skipping conversion.")
        return True
        
    script_path = os.path.abspath(sys.path[0])
    command = [
        "python", os.path.join(script_path, "trans_partig.py"),
        "-fai", f"{fa_file}.fai",
        "-p", input_partig_file,
        "-o", output_csv_file,
        "-r", RE_file, 
        "-d", f"{output_prefix}.digraph.csv"
    ]
    return run_command(command, "Converting Partig output to CSV", logger)

def run_pipeline_allele(split_gfa_file: str, HiC_file: str, RE_file: str, partig_file: str, output_prefix: str, n_chr: int, logger: logging.Logger) -> bool:
    """Run the Pipeline_allele clustering step."""
    if not all(check_file_exists_and_not_empty(f, logger) for f in [split_gfa_file, HiC_file, RE_file, partig_file]):
        return False
        
    output_check_file = f"{output_prefix}.allele.cluster.expand.txt"
    if check_file_exists_and_not_empty(output_check_file, logger):
        logger.info(f"Pipeline_allele output already exists: {output_check_file}. Skipping Pipeline_allele.")
        return True

    try:
        logger.info("Starting: Running Pipeline_allele")
        # Assuming Pipeline_allele is a class/function that executes the pipeline
        Pipeline_allele(split_gfa_file, HiC_file, RE_file, partig_file, output_prefix, n_chr)
        
        if not check_file_exists_and_not_empty(output_check_file, logger):
            logger.error("Pipeline_allele completed, but output file is missing or empty.")
            return False
            
        logger.info("Completed: Running Pipeline_allele")
        return True
    except Exception as e:
        logger.error(f"Error in: Running Pipeline_allele: {e}")
        return False

def run_trans_cluster(output_prefix: str, logger: logging.Logger) -> bool:
    """Run trans_cluster.py to transform allele cluster results."""
    input_file = f"{output_prefix}.allele.cluster.expand.txt"
    if not check_file_exists_and_not_empty(input_file, logger):
        return False

    subgraphs_file = "group_ctgs_save.txt"
    script_path = os.path.abspath(sys.path[0])
    command = [
        "python", os.path.join(script_path, "trans_cluster.py"),
        "-c", input_file,
        "-s", subgraphs_file,
        "-n", output_prefix
    ]
    return run_command(command, "Running trans_cluster.py", logger)

def run_pipeline_chr(output_prefix: str, HiC_file: str, logger: logging.Logger) -> bool:
    """Run pipeline_chr.py for chromosome clustering pre-processing."""
    input_file_c = f"{output_prefix}.allele.cluster.expand.txt"
    input_file_s = "group_ctgs_save.txt"
    if not all(check_file_exists_and_not_empty(f, logger) for f in [input_file_c, HiC_file, input_file_s]):
        return False

    script_path = os.path.abspath(sys.path[0])
    command = [
        "python", os.path.join(script_path,"pipeline_chr.py"),
        "-c", input_file_c,
        "-l", HiC_file,
        "-s", input_file_s,
        "-n", output_prefix
    ]
    return run_command(command, "Running pipeline_chr.py (HiC clustering pre-processing)", logger)

def run_multilevel_cluster_optimized(input_hic_file: str, output_prefix: str, chr_number: int, logger: logging.Logger, r_min=0.01, r_max=3.0, tolerance=0.001, max_iter=50) -> bool:
    """
    Run multilevel clustering with an optimized binary search for the ratio 'r'
    to achieve the desired cluster count.
    """
    input_RE_counts = f"{output_prefix}.RE_counts.txt"
    chr_cluster_output = f"{output_prefix}.chr.cluster.ctg.txt"
    allele_cluster = f"{output_prefix}.allele.cluster.ctg.txt"
    
    if not all(check_file_exists_and_not_empty(f, logger) for f in [input_hic_file, input_RE_counts]):
        return False
        
    r = 1.0
    r_start = r_min
    r_end = r_max
    
    # for i in range(max_iter):
    while r_start <= r <= r_end:

        if r_end - r_start < tolerance:
            logger.warning(f"r adjustment range is too small ({r_end - r_start:.4f} < {tolerance}). Stopping search.")
            break
        
        logger.info(f"Trying r={r:.4f} (range: [{r_start:.4f}, {r_end:.4f}])")
        
        try:
            cluster_count = Multilevel_cluster(
                input_hic_file,
                chr_cluster_output,
                float(r),
                True,
                input_RE_counts,
                allele_cluster,
                int(chr_number)
            )
            if cluster_count is None:
                logger.error("Multilevel_cluster returned None, check inner component errors.")
                return False
                
            logger.info(f"Current cluster count: {cluster_count}")

            if cluster_count == chr_number:
                logger.info(f"Success: Desired cluster count {chr_number} achieved with r={r:.4f}.")
                return True
            elif cluster_count > chr_number:
                logger.info(f"Cluster count {cluster_count} > target {chr_number}. Decreasing r (r_end = {r:.4f}).")
                r_end = r
                r = (r_start + r_end) / 2
            else: # cluster_count < chr_number
                logger.info(f"Cluster count {cluster_count} < target {chr_number}. Increasing r (r_start = {r:.4f}).")
                r_start = r
                r = (r_start + r_end) / 2
                
        except Exception as e:
            logger.error(f"Error running Multilevel_cluster with r={r:.4f}: {e}")
            return False

    logger.error("Failed to achieve the desired cluster count after all attempts/iterations.")
    return False

def run_rescue_base_subgraph(HiC_file: str, output_prefix: str, logger: logging.Logger) -> bool:
    """Run rescue_base_subgraph.py to rescue unclustered contigs."""
    input_chr_cluster = f"{output_prefix}.chr.cluster.ctg.txt"
    input_allele_cluster = f"{output_prefix}.allele.cluster.expand.txt"
    input_s = "group_ctgs_save.txt"
    
    if not all(check_file_exists_and_not_empty(f, logger) for f in [input_chr_cluster, input_allele_cluster, input_s, HiC_file]):
        return False

    # The final output is often written to the input_chr_cluster or a new file, 
    # but we rely on the command's success.
        
    script_path = os.path.abspath(sys.path[0])
    command = [
        "python", os.path.join(script_path,"rescue_base_subgraph.py"),
        "--chr_cluster_ctg", input_chr_cluster,
        "--allele_cluster", input_allele_cluster,
        "-s", input_s,
        "-rs", "group_ctgs_filter.txt",
        "-l", HiC_file
    ]
    return run_command(command, "Running rescue_base_subgraph.py", logger)

def filter_edges_by_density(chr_num: int, HiC_file: str, group_ctgs_save: str, filter_HiC_file: str, logger: logging.Logger, filter_threshold: int = 30, step: float = 0.5) -> None:
    """Filter HiC edges based on graph density and a link threshold."""
    if not all(check_file_exists_and_not_empty(f, logger) for f in [HiC_file, group_ctgs_save]):
        # If input files are missing, we can't filter, but we might still copy the original if possible.
        logger.warning(f"Missing input for filtering, attempting to copy original HiC file.")
        try:
             import shutil
             shutil.copy(HiC_file, filter_HiC_file)
             return
        except Exception as e:
             logger.error(f"Failed to copy original HiC file to filter file: {e}")
             return

    try:
        nodes = pd.read_csv(group_ctgs_save, sep='\t' ,header=None)
        edges = pd.read_csv(HiC_file, sep=',', header=0)
        edges.columns = ['source', 'target', 'links']

        num_nodes = len(nodes)
        num_edges = len(edges)

        density = (2 * num_edges) / (num_nodes * (num_nodes - 1)) if num_nodes > 1 else 0
        logger.info(f"HiC graph info : edges -> {num_edges}, nodes -> {num_nodes}, density -> {density:.4f}")

        # Check conditions for filtering
        # The condition uses an OR logic in the original: density < 0.2 OR num_nodes < (num_edges / filter_threshold)
        # Using the simplified check from the original for robustness unless the logic is confirmed to be an AND
        
        # Original logic: if density < 0.2 and num_nodes < (num_edges / filter_threshold):
        # We'll use the original logic and ensure variables are checked
        
        should_filter = (density < 0.2 and num_nodes < (num_edges / filter_threshold)) or (num_nodes < (num_edges / filter_threshold))
        
        if should_filter:
            logger.info(f"HiC signal filtering triggered. Target filter threshold: {filter_threshold} edges per node.")
            threshold = 0.5
            
            max_threshold_search = 100 
            search_count = 0
            
            while True:
                if search_count >= max_threshold_search:
                    logger.warning(f"Exceeded max threshold search attempts ({max_threshold_search}). Stopping filtering search.")
                    break
                    
                filtered_edges = edges[edges['links'] > threshold]
                filtered_num_edges = len(filtered_edges)

                logger.info(f"HiC signal filtering : threshold -> {threshold:.1f}\t edges -> {filtered_num_edges}")

                if filtered_num_edges < (num_nodes * filter_threshold):
                    logger.info(f"Achieved target edge count ({filtered_num_edges} < {num_nodes * filter_threshold}).")
                    break

                threshold += step
                search_count += 1
                
            filtered_edges.to_csv(filter_HiC_file, sep=',', header=True, index=False)
            logger.info(f"Filtered HiC written to {filter_HiC_file}")
        else:
            edges.to_csv(filter_HiC_file, sep=',', header=True, index=False)
            logger.info(f"No filtering required. Original HiC written to {filter_HiC_file}")
            
    except Exception as e:
        logger.error(f"Error during HiC edge filtering: {e}. Writing original edges as fallback.")
        try:
             edges.to_csv(filter_HiC_file, sep=',', header=True, index=False)
        except:
             logger.error(f"Fallback to writing original HiC also failed.")



def run_spectral_clustering_fallback(input_hic_file: str, groups_file: str, final_cluster_output: str, chr_number: int, logger: logging.Logger) -> bool:
    """
    Runs Spectral Clustering on the HiC link file as a fallback when multilevel clustering fails.
    """
    logger.info("Starting: Spectral Clustering Fallback for Chromosome Clustering.")

    def read_c(c):
        cluster_dict = defaultdict(list)
        subgraph_group_dict = defaultdict(list)
        with open(c, 'r') as file:
            for line in file:
                line = line.strip().split('\t')
                for utg in line[2].split():
                    cluster_dict[line[0]].append(utg)
                    subgraph_group_dict[utg].append(line[0])
        return cluster_dict, subgraph_group_dict

    try:
        cluster_dict, subgraph_group_dict = read_c(groups_file)
        logger.info(f"Loaded file {groups_file}.")
    except Exception as e:
        logger.error(f"Failed to load nodes from {groups_file}: {e}")
        return False
    
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

        logger.info(f"Running SpectralClustering with k={chr_number}.")
        sc = SpectralClustering(
            n_clusters=chr_number,
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

        chr_cluster_dict = defaultdict(set)
        for idx, group in summary_dict.items():
            for group_member in group:
                for utg in cluster_dict.get(group_member, []):
                    chr_cluster_dict[idx].add(utg)

        group = 0
        with open(final_cluster_output, 'w') as file:
            for idx in chr_cluster_dict:
                group += 1
                utgs = list(chr_cluster_dict[idx])
                file.write(f"group{group}\t{len(utgs)}\t{' '.join(utgs)}\n")
        
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


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog='cluster_chr', description='A pipeline for chromosome clustering using Hi-C data and genome assembly.')

    base_group  = parser.add_argument_group('>>> Parameters for basic data')
    base_group.add_argument("-f", "--fa_file", metavar='\b', required=True, help="Path to the assembly FASTA file.")
    base_group.add_argument("-r", "--RE_file", metavar='\b', required=True, help="Path to the restriction enzyme file.")

    hic_group  = parser.add_argument_group('>>> Parameters for HiC data alignment')
    hic_group.add_argument("-l", "--HiC_file", metavar='\b', required=True, help="Path to the Hi-C data file.")

    output_group  = parser.add_argument_group('>>> Parameter for the prefix of the result file')
    output_group.add_argument("-op", "--output_prefix", metavar='\b', required=True, help="Prefix for output files.")

    genome_group  = parser.add_argument_group('>>> Parameters of chromosome numbers')
    genome_group.add_argument("-n_chr", "--chr_number", type=int, metavar='\b', required=True, help="Desired number of clusters corresponding to chromosomes.")

    partig_group  = parser.add_argument_group('>>> Parameter for partig')
    partig_group.add_argument("-pk", "--partig_k", metavar='\b', type=int, default=17, help="K-mer size for Partig. Default: 17.")
    partig_group.add_argument("-pw", "--partig_w", metavar='\b', type=int, default=17, help="Minimizer window size for Partig. Default: 17.")
    partig_group.add_argument("-pc", "--partig_c", metavar='\b', type=int, default=60, help="Max occurrance for Partig. Default: 60.")
    partig_group.add_argument("-pm", "--partig_m", metavar='\b', type=float, default=0.95, help="Mini K-mer similarity for Partig. Default: 0.95.")

    split_GFA_group  = parser.add_argument_group('>>> Parameter for split GFA')
    split_GFA_group.add_argument("-g", "--gfa", metavar='\b', required=True, help="Path to the GFA file.")
    split_GFA_group.add_argument("-n", "--split_gfa_n", metavar='\b', type=int, default=2, help="Number of common neighbors when splitting GFA. Default: 2.")
    split_GFA_group.add_argument("-i", "--split_gfa_iter", metavar='\b', type=int, default=3, help="Number of iterations when splitting GFA. Default: 3.")

    clustering_group  = parser.add_argument_group('>>> Parameter for clustering')
    clustering_group.add_argument("-r_max", "--r_max", metavar='\b', default=3, help="Maximum value of parameter R during Louvain clustering.")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    
    # Robustness: Input file existence checks
    for file_arg in ['fa_file', 'RE_file', 'HiC_file', 'gfa']:
        file_path = getattr(args, file_arg)
        if not os.path.exists(file_path):
            sys.exit(f"Error: Required input file not found: {file_path}")
            
    # Robustness: Integer/Float validation (already handled by argparse type, but check logical values)
    if args.chr_number <= 0:
        sys.exit(f"Error: Chromosome number (-n_chr) must be a positive integer.")
    if args.partig_k <= 0 or args.partig_w <= 0:
        sys.exit(f"Error: Partig k-mer/window size (-pk, -pw) must be positive.")
    if not (0.80 <= args.partig_m <= 1.0):
        sys.exit(f"Error: Partig k-mer similarity (-pm) must be between 0.80 and 1.0.")

    return args

def main():
    args = parse_arguments()
    pwd = os.getcwd()
    logger = setup_logging('cluster_chr.log')
    log_start(logger, "cluster_chr.py", "1.0.0", args)
    successful_clustering = False

    # --- Step 0: Pre-processing and Initial GFA Split ---
    if not index_fasta(args.fa_file, logger):
        logger.error("Indexing FASTA failed. Exiting.")
        return

    if not split_gfa(args.gfa, args.split_gfa_n, args.split_gfa_iter, args.output_prefix, logger):
        logger.error("Split GFA failed. Exiting.")
        return

    # --- Step 1: Running partig ---
    min_partig_m = 0.8
    if args.partig_m < min_partig_m:
        logger.error("Your parameter (partig_m) is too small. Please select a value greater than or equal to 0.8 and less than or equal to 0.95.. Exiting.")
        return

    partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{min_partig_m:.2f}.txt"
    trans_partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{min_partig_m:.2f}.csv"
    
    if not check_file_exists_and_not_empty(partig_file, logger):
        if not run_partig(args.fa_file, args.partig_k, args.partig_w, args.partig_c, min_partig_m, args.output_prefix, logger):
            logger.error("Run partig failed. Exiting.")
            return
    if not check_file_exists_and_not_empty(trans_partig_file, logger):
        if not convert_partig_output(args.fa_file, args.partig_k, args.partig_w, args.partig_c, min_partig_m, args.output_prefix, args.RE_file, logger):
            logger.error("Conversion partig output failed. Exiting.")
            return

    vals = [x / 100 for x in range(int(args.partig_m*100), 75, -5)]

    for partig_m_iter in vals:
        partig_m_iter = float(partig_m_iter)
        logger.info(f"The current parameter partig_m is {partig_m_iter}...")
        trans_partig_df = pd.read_csv(trans_partig_file)
        filtered_trans_partig_df = trans_partig_df[trans_partig_df.iloc[:, 2] >= partig_m_iter]

        trans_partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{partig_m_iter:.2f}.csv"
        filtered_trans_partig_df.to_csv(trans_partig_file, index=False)

        # --- Step 2: Allele Clustering Pipeline ---
        split_gfa_file = f"{args.output_prefix}.rmTip.split.gfa"
        if not run_pipeline_allele(split_gfa_file, args.HiC_file, args.RE_file, trans_partig_file, args.output_prefix, args.chr_number, logger):
            logger.error("Run allele_cluster pipeline failed. Exiting.")
            return

        # --- Step 3: Transform Allele Cluster Results ---
        if not run_trans_cluster(args.output_prefix, logger):
            logger.error("Transform cluster results failed. Exiting.")
            return

        # --- Step 4: Chromosome Clustering Pre-processing ---
        if not run_pipeline_chr(args.output_prefix, args.HiC_file, logger):
            logger.error("Run chromosome clustering pre-processing failed. Exiting.")
            return
        
        # --- Step 5: Filter Allele HiC Data ---
        filter_threshold = 30
        allele_hic_file, filter_HiC_file = f"{args.output_prefix}.allele.hic.csv", f"{args.output_prefix}.allele.hic.filter.csv"
        groups_file = "group_ctgs_save.txt"
        
        # Ensure the required HiC file for filtering exists after run_pipeline_chr
        if not check_file_exists_and_not_empty(allele_hic_file, logger):
            logger.error(f"Required HiC file for filtering not found: {allele_hic_file}. Exiting.")
            return

        # --- Step 6: Multilevel Chromosome Clustering ---
        max_filter_attempts = 6
        successful_clustering = False
        
        for attempt in range(max_filter_attempts):
            if attempt > 0:
                input_hic_file = filter_HiC_file
                output_hic_file = filter_HiC_file
            else:
                input_hic_file = allele_hic_file
                output_hic_file = filter_HiC_file
            current_threshold = filter_threshold - (attempt * 5)
            if current_threshold < 5: # Set a minimum practical threshold
                logger.warning("Minimum filter threshold reached. Stopping threshold search.")
                break
                
            logger.info(f"Attempt {attempt + 1}: Chromosome clustering uses the filter_threshold : {current_threshold}.")
            
            filter_edges_by_density(
                args.chr_number, 
                input_hic_file, 
                groups_file, 
                output_hic_file, 
                logger, 
                filter_threshold=current_threshold, 
                step=0.5
            )
            
            # Check if the filtered HiC file is valid before proceeding
            if not check_file_exists_and_not_empty(filter_HiC_file, logger):
                logger.info(f"Filtered HiC file is invalid: {filter_HiC_file}. Cannot proceed to clustering with this threshold.")
                filter_HiC_file = allele_hic_file
                continue
                
            # Run optimized multilevel cluster
            if run_multilevel_cluster_optimized(filter_HiC_file, args.output_prefix, int(args.chr_number), logger, r_min=0.01, r_max=float(args.r_max)):
                successful_clustering = True
                break
            else:
                logger.warning(f"Run Chromosome multilevel_cluster failed with filter_threshold : {current_threshold}. Trying next threshold.")
        
        if successful_clustering:
            break
    
    final_cluster_file = f"{args.output_prefix}.chr.cluster.ctg.txt"
    if not successful_clustering:
        logger.error("Louvain failed to cluster to the provided number of chromosomes, and the clustering subgraphs failed. Spectral Clustering will be used to cluster the chromosomes. Alternatively, you can use a reference genome to cluster the chromosomes yourself, and then use the subsequent steps of GPhase.")
        
        # run Spectral Clustering
        groups_file = f"{args.output_prefix}.allele.cluster.ctg.txt"
        final_cluster_file = f"{args.output_prefix}.chr.cluster.ctg.txt"
        if run_spectral_clustering_fallback(allele_hic_file, groups_file, final_cluster_file, int(args.chr_number), logger):
             successful_clustering = True
             logger.info("Spectral Clustering fallback completed successfully.")
        else:
             logger.error("Spectral Clustering fallback also failed. Exiting pipeline.")
             return

    # Final check on cluster count
    if not check_file_exists_and_not_empty(final_cluster_file, logger):
        logger.error("Final chromosome cluster file is missing or empty after successful clustering flag.")
        return
        
    try:
        with open(final_cluster_file, "r") as f:
            line_count = sum(1 for _ in f)
    except Exception as e:
        logger.error(f"Error reading final cluster file {final_cluster_file}: {e}")
        return

    if line_count != int(args.chr_number):
        logger.error(f"Clustering Chromosome Error: Incorrect number of clusters ({line_count} found, {args.chr_number} desired).")
        # Log this as an error but might proceed if the downstream step can handle it, 
        # though typically this means the pipeline failed. We exit for high robustness.
        return

    # --- Step 7: Rescue Base Subgraph ---
    if not run_rescue_base_subgraph(args.HiC_file, args.output_prefix, logger):
        logger.error("Rescue base subgraph failed.")
        return
    
    if check_file_exists_and_not_empty("rescue.cluster.ctg.txt", logger):
        logger.info(f"The clustering results of the {args.chr_number} chromosomes in gphase are stored in file: rescue.cluster.ctg.txt .")

    logger.info("Chromosome clustering was successful; haplotype clustering will be performed next.")

if __name__ == "__main__":
    main()