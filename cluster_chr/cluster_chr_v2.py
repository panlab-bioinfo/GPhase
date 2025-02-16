# PYTHON_ARGCOMPLETE_OK

import os
import subprocess
import argparse
import logging
import argcomplete
import sys
from argcomplete.completers import FilesCompleter

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
    
    # Formatter with timestamp and script name (based on your example format)
    formatter = logging.Formatter('%(asctime)s <%(module)s.py> [%(funcName)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def run_command(command, description, logger):
    try:
        logger.info(f"Starting: {description}")
        logger.info(f"Command: {' '.join(command)}")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logger.info(f"Running result: {result.stdout}")
        logger.info(f"Completed: {description}\n")
        return 1
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in: {description}\nCommand failed: {' '.join(command)}")
        logger.error(f"Error output: {e.stderr}")
        return None

def index_fasta(asm_fa, logger):
    run_command(["samtools", "faidx", asm_fa], "Indexing the assembly FASTA",logger)

def split_gfa(gfa, split_gfa_n, split_gfa_iter, output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path, "split_GFA.py"),
        "-g", gfa,
        "-n", str(split_gfa_n),
        "-iter", str(split_gfa_iter),
        "-o", output_prefix
    ], "Splitting the GFA file",logger)


def run_partig(asm_fa, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
    output_file = f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt"
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "../bin/partig")
    try:
        with open(output_file, "w") as outfile:
            command = [
                script_path_add, asm_fa,
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


def convert_partig_output(asm_fa, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path, "trans.partig.py"),
        "-fai", f"{asm_fa}.fai",
        "-p", f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt",
        "-o", f"{output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.csv"
    ], "Converting Partig output to CSV",logger)

def run_pipeline_allele(split_gfa_file, HiC_file, RE_file, partig_file, output_prefix,logger):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path, "pipeline.allele.py"),
        "-g", split_gfa_file,
        "-l", HiC_file,
        "-r", RE_file,
        "-a", partig_file,
        "-n", output_prefix
    ], "Running pipeline.allele.py", logger)

def run_trans_cluster(output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path, "trans_cluster.py"),
        "-c", f"{output_prefix}.allele.cluster.expand.txt",
        "-s", "group_ctgs_save.txt",
        "-n", output_prefix
    ], "Running trans_cluster.py", logger)

def run_pipeline_chr(output_prefix, HiC_file):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path,"pipeline.chr.py"),
        "-c", f"{output_prefix}.allele.cluster.expand.txt",
        "-l", HiC_file,
        "-s", "group_ctgs_save.txt",
        "-n", output_prefix
    ], "Running pipeline.chr.py", logger)

def get_cluster_count(output_prefix, logger):
    cluster_file = f"{output_prefix}.chr.cluster.txt"
    if not os.path.exists(cluster_file):
        logger.error(f"Cluster file {cluster_file} does not exist.")
        return None
    with open(cluster_file, 'r') as f:
        clusters = set(line.strip().split()[0] for line in f if line.strip())
    return len(clusters)

def run_multilevel_cluster(output_prefix, chr_number, logger):
    r_min, r_max = 0.01, 5.0
    r = 1.0
    while r_min <= r <= r_max:
        logger.info(f"Running multilevel clustering with r={r}")
        script_path = os.path.abspath(sys.path[0])
        result = run_command([
            "python", os.path.join(script_path,"multilevel_cluster.py"),
            "-c", f"{output_prefix}.allele.hic.csv",
            "-o", f"{output_prefix}.chr.cluster.txt",
            "-r", str(r)
        ], "Running multilevel_cluster.py", logger)

        if result is None:
            return False

        cluster_count = get_cluster_count(output_prefix)
        if cluster_count is None:
            return False

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
    return False

def run_trans_allele_cluster(output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path,"trans_allele_cluster.py"),
        "-c1", f"{output_prefix}.chr.cluster.txt",
        "-c2", f"{output_prefix}.allele.cluster.ctg.txt",
        "-n", output_prefix
    ], "Running trans_allele_cluster.py", logger)

def run_rescue_base_subgraph(HiC_file, output_prefix, logger):
    script_path = os.path.abspath(sys.path[0])
    run_command([
        "python", os.path.join(script_path,"rescue_base_subgraph.py"),
        "--chr_cluster_ctg", f"{output_prefix}.chr.cluster.ctg.txt",
        "--allele_cluster", f"{output_prefix}.allele.cluster.expand.txt",
        "-s", "group_ctgs_save.txt",
        "-rs", "group_ctgs_filter.txt",
        "-l", HiC_file
    ], "Running rescue_base_subgraph.py", logger)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pipeline for clustering chromosomes.")

    # Required arguments
    parser.add_argument("-a", "--asm_fa", required=True, help="Path to the assembly FASTA file.")
    parser.add_argument("-g", "--gfa", required=True, help="Path to the GFA file.")
    parser.add_argument("-l", "--HiC_file", required=True, help="Path to the Hi-C data file.")
    parser.add_argument("-e", "--RE_file", required=True, help="Path to the restriction enzyme file.")
    parser.add_argument("-op", "--output_prefix", required=True, help="Prefix for output files.")
    parser.add_argument("-cn", "--chr_number", type=int, required=True, help="Desired number of clusters corresponding to chromosomes.")

    # Optional arguments with default values
    parser.add_argument("-k", "--partig_k", type=int, default=17, help="K-mer size for Partig. Default: 17.")
    parser.add_argument("-w", "--partig_w", type=int, default=17, help="Minimizer window size for Partig. Default: 17.")
    parser.add_argument("-c", "--partig_c", type=int, default=60, help="Max occurrance for Partig. Default: 60.")
    parser.add_argument("-m", "--partig_m", type=float, default=0.95, help="Mini k-mer similarity for Partig. Default: 0.95.")
    parser.add_argument("-n", "--split_gfa_n", type=int, default=2, help="Number of common neighbors when splitting GFA. Default: 2.")
    parser.add_argument("-i", "--split_gfa_iter", type=int, default=3, help="Number of iterations when splitting GFA. Default: 3.")

    return parser.parse_args()

def main():

    logger = setup_logging('cluster_chr.log')
    args = parser.parse_args()
    argcomplete.autocomplete(parser)
    pwd = os.getcwd()

    if not index_fasta(args.asm_fa, logger):
        logger.error("Get fasta index Error: Samtools gets fasta file index error.")
        return

    if not split_gfa(args.gfa, args.split_gfa_n, args.split_gfa_iter, args.output_prefix, logger):
        logger.error("Split GFA Error: An error occurred while splitting the GFA.")
        return

    split_gfa_file = f"{args.output_prefix}.rmTip.split.gfa"

    if not run_partig(args.asm_fa, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix,logger):
        logger.error("Run partig Error: An error occurred while running Partig.")
        return

    if not convert_partig_output(args.asm_fa, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix, logger):
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

    if not run_multilevel_cluster(args.output_prefix, args.chr_number, logger):
        logger.error("Run multilevel_cluster Error: An error occurred while running multilevel_cluster.")
        return

    if not run_trans_allele_cluster(args.output_prefix, logger):
        logger.error("Transform Error: An error occurred while converting hic_cluster results.")
        return

    if not run_rescue_base_subgraph(args.HiC_file, args.output_prefix, logger):
        logger.error("Rescue Error: An error occurred while rescue base subgraph.")
        return

if __name__ == "__main__":
    main()
