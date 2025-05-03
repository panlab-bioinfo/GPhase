#!/usr/bin/env python3

import os
import subprocess
import argparse
import logging
import argcomplete
import numpy as np
import pandas as pd
import sys
import random
from multilevel_cluster_v2 import Multilevel_cluster
from argcomplete.completers import FilesCompleter

def setup_logging(log_file: str = "cluster_chr.log") -> logging.Logger:
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

def log_start(logger: logging.Logger, script_name: str, version: str, args: argparse.Namespace):
    """Log the start of the program."""
    logger.info(f"Program started, {script_name} version: {version}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")

def run_command(command, description, logger):
    try:
        logger.info(f"Starting: {description}")
        logger.info(f"Command: {' '.join(command)}")
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logger.info(f"Completed: {description}\n")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in: {description}\n")
        logger.error(f"Command failed: {' '.join(command)}")
        return None

def index_fasta(fa_file, logger):
    flag = run_command(["samtools", "faidx", fa_file], "Indexing the assembly FASTA",logger)
    if flag:
        return True
    return False 


def split_gfa(gfa, split_gfa_n, split_gfa_iter, output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    flag = run_command([
        "python", os.path.join(script_path, "split_GFA.py"),
        "-g", gfa,
        "-n", str(split_gfa_n),
        "-iter", str(split_gfa_iter),
        "-o", output_prefix
    ], "Splitting the GFA file",logger)
    if flag:
        return True
    return False 


def run_partig(fa_file, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
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
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in: Running Partig\nCommand failed: {' '.join(command)}")
        logger.error(f"Error output: {e.stderr}")
        return False


def convert_partig_output(fa_file, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
    script_path = os.path.abspath(sys.path[0])
    flag = run_command([
        "python", os.path.join(script_path, "trans.partig.py"),
        "-fai", f"{fa_file}.fai",
        "-p", f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt",
        "-o", f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.csv"
    ], "Converting Partig output to CSV",logger)
    if flag:
        return True
    return False 

def run_pipeline_allele(split_gfa_file, HiC_file, RE_file, partig_file, output_prefix,logger):
    script_path = os.path.abspath(sys.path[0])
    flag = run_command([
        "python", os.path.join(script_path, "pipeline.allele.py"),
        "-g", split_gfa_file,
        "-l", HiC_file,
        "-r", RE_file,
        "-a", partig_file,
        "-n", output_prefix
    ], "Running pipeline.allele.py", logger)
    if flag:
        return True
    return False 

def run_trans_cluster(output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    flag = run_command([
        "python", os.path.join(script_path, "trans_cluster.py"),
        "-c", f"{output_prefix}.allele.cluster.expand.txt",
        "-s", "group_ctgs_save.txt",
        "-n", output_prefix
    ], "Running trans_cluster.py", logger)
    if flag:
        return True
    return False 

def run_pipeline_chr(output_prefix, HiC_file, logger):
    script_path = os.path.abspath(sys.path[0])
    flag = run_command([
        "python", os.path.join(script_path,"pipeline.chr.py"),
        "-c", f"{output_prefix}.allele.cluster.expand.txt",
        "-l", HiC_file,
        "-s", "group_ctgs_save.txt",
        "-n", output_prefix
    ], "Running pipeline.chr.py", logger)
    if flag:
        return True
    return False 


def run_multilevel_cluster_v1(output_prefix, chr_number, logger):
    r_min, r_max = 0.01, 30.0
    r = 1
    while r_min <= r <= r_max:
        logger.info(f"Running multilevel clustering with r={r}")
        script_path = os.path.abspath(sys.path[0])


        cluster_count = Multilevel_cluster(f"{output_prefix}.allele.hic.csv", f"{output_prefix}.chr.cluster.ctg.txt", str(r), True, f"{output_prefix}.RE_counts.txt", f"{output_prefix}.allele.cluster.ctg.txt", chr_number)

        if cluster_count is None:
            return r

        logger.info(f"Current cluster count: {cluster_count}")

        if cluster_count == chr_number:
            logger.info(f"Desired cluster count {chr_number} achieved with r={r}.")
            return True
        elif cluster_count > chr_number:
            r_max = r
            r = (r_min + r) / 2
        else:
            r_min = r
            r = (r + r_max) / 2

        if abs(r_max - r_min) < 0.01:
            logger.warning("r adjustment range is too small to continue.")
            break

    logger.error("Failed to achieve the desired cluster count within the r range.")
    return 0

def run_multilevel_cluster(output_prefix, chr_number, logger, max_attempts=50):
    r_min_global, r_max_global = 0.01, 30.0

    for attempt in range(max_attempts):
        # 每次尝试从全局范围中随机选择一个初始r
        r = random.uniform(r_min_global, r_max_global)
        r_min, r_max = r_min_global, r_max_global

        logger.info(f"Attempt {attempt + 1}: Trying initial r={r}")

        while r_min <= r <= r_max:
            logger.info(f"Running multilevel clustering with r={r}")
            script_path = os.path.abspath(sys.path[0])

            cluster_count = Multilevel_cluster(
                f"{output_prefix}.allele.hic.filter.csv",
                f"{output_prefix}.chr.cluster.ctg.txt",
                str(r),
                True,
                f"{output_prefix}.RE_counts.txt",
                f"{output_prefix}.allele.cluster.ctg.txt",
                chr_number
            )

            if cluster_count is None:
                logger.warning("Multilevel_cluster returned None, exiting current attempt.")
                break

            logger.info(f"Current cluster count: {cluster_count}")

            if cluster_count == chr_number:
                logger.info(f"Success: Desired cluster count {chr_number} achieved with r={r}.")
                return True
            elif cluster_count > chr_number:
                r_max = r
                r = (r_min + r) / 2
            else:
                r_min = r
                r = (r + r_max) / 2

            if abs(r_max - r_min) < 0.01:
                logger.warning("r adjustment range is too small to continue in this attempt.")
                break

    logger.error("Failed to achieve the desired cluster count after all attempts.")
    return False

def run_rescue_base_subgraph(HiC_file, output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    flag = run_command([
        "python", os.path.join(script_path,"rescue_base_subgraph.py"),
        "--chr_cluster_ctg", f"{output_prefix}.chr.cluster.ctg.txt",
        "--allele_cluster", f"{output_prefix}.allele.cluster.expand.txt",
        "-s", "group_ctgs_save.txt",
        "-rs", "group_ctgs_filter.txt",
        "-l", HiC_file
    ], "Running rescue_base_subgraph.py", logger)
    if flag:
        return True
    return False 

def filter_edges_by_density(chr_num, HiC_file, group_ctgs_save, filter_HiC_file, logger, step=0.5):

    nodes = pd.read_csv(group_ctgs_save, sep='\t' ,header=None)
    edges = pd.read_csv(HiC_file, sep=',', header=0)
    edges.columns = ['source', 'target', 'links']

    num_nodes = len(nodes)
    num_edges = len(edges)
    density = (2 * num_edges) / (num_nodes * (num_nodes - 1)) if num_nodes > 1 else 0
    logger.info(f"Chr{chr_num} hic graph info : edges -> {num_edges}\t nodes -> {num_nodes}\t density -> {density}")

    if density < 0.2 and num_nodes < (num_edges / 30):
        logger.info(f"Chr{chr_num} HiC signal filtering...")
        threshold = 0.5
        
        while True:
            filtered_edges = edges[edges['links'] > threshold]
            filtered_num_edges = len(filtered_edges)

            logger.info(f"Chr{chr_num} HiC signal filtering : c -> {threshold:.1f}\t edges: -> {filtered_num_edges}")

            if filtered_num_edges < (num_nodes * 30):
                break

            threshold += step

        filtered_edges.to_csv(filter_HiC_file, sep=',', header=True, index=False)
    else:
        edges.to_csv(filter_HiC_file, sep=',', header=True, index=False)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog='cluster_chr')

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
    split_GFA_group.add_argument("-n", "--split_gfa_n", metavar='\b', type=int, default=5, help="Number of common neighbors when splitting GFA. Default: 5.")
    split_GFA_group.add_argument("-i", "--split_gfa_iter", metavar='\b', type=int, default=3, help="Number of iterations when splitting GFA. Default: 3.")

    return parser.parse_args()

def main():
    args = parse_arguments()
    argcomplete.autocomplete(args)
    pwd = os.getcwd()
    logger = setup_logging('cluster_chr.log')
    log_start(logger, "cluster_chr.py", "1.0.0", args)

    if not index_fasta(args.fa_file, logger):
        logger.error("Split GFA Error: An error occurred while splitting the GFA.")
        return


    if not split_gfa(args.gfa, args.split_gfa_n, args.split_gfa_iter, args.output_prefix, logger):
        logger.error("Split GFA Error: An error occurred while splitting the GFA.")
        return

    split_gfa_file = f"{args.output_prefix}.rmTip.split.gfa"

    if not run_partig(args.fa_file, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix,logger):
        logger.error("Run partig Error: An error occurred while running Partig.")
        return

    if not convert_partig_output(args.fa_file, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix, logger):
        logger.error("Conversion partig Error: in Partig output to CSV.")
        return

    partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{args.partig_m}.csv"

    if not run_pipeline_allele(split_gfa_file, args.HiC_file, args.RE_file, partig_file, args.output_prefix, logger):
        logger.error("Run allele_cluster Error: An error occurred while running allele_cluster.")
        return

    if not run_trans_cluster(args.output_prefix, logger):
        logger.error("Transform Error: An error occurred while converting allele_cluster results.")
        return

    if not run_pipeline_chr(args.output_prefix, args.HiC_file, logger):
        logger.error("Run hic_cluster Error: An error occurred while running hic_cluster.")
        return
    
    # filter allele HiC
    allele_hic_file = f"{args.output_prefix}.allele.hic.csv"
    groups_file = "group_ctgs_save.txt"
    filter_HiC_file = f"{args.output_prefix}.allele.hic.filter.csv"
    filter_edges_by_density(args.chr_number, allele_hic_file, groups_file, filter_HiC_file, logger, step=0.5)

    if not run_multilevel_cluster(args.output_prefix, args.chr_number, logger):
        logger.error("Run multilevel_cluster Error: An error occurred while running multilevel_cluster.")
        return


    if not run_rescue_base_subgraph(args.HiC_file, args.output_prefix, logger):
        logger.error("Rescue Error: An error occurred while rescue base subgraph.")
        return

if __name__ == "__main__":
    main()
