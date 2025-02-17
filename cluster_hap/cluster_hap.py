import os
import argparse
import subprocess
import logging
import argcomplete
from argcomplete.completers import FilesCompleter
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor
import functools
import sys

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
    performance_group.add_argument("--process_num", metavar='\b',type=int, default=10,
                       help=f"Number of parallel processes (default: 10)")
    performance_group.add_argument("--thread_num", metavar='\b',type=int,default=10,
                       help=f"Number of threads per process (default: 10)")

    Optional_group  = parser.add_argument_group('>>> Optional parameters')
    Optional_group.add_argument('--correct', action='store_true', help='correct the cluster')
    Optional_group.add_argument('--isolated_threshold', default=5,help='<int>Detect whether the intensity of the hic signal is an outlier')

    return parser.parse_args()

def execute_command(command, error_message,logger):
    try:
        logger.info(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True, executable='/bin/bash')
        logger.info(f"Running command successfully: {command}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(command)}")
        logger.error(f"Error: {e}")
    except Exception as e:
        logger.error(f"{error_message}")
        logger.error(f"Error executing command: {e}")
        return False


def create_symlink(source, dest,logger):
    if not os.path.exists(dest):
        os.symlink(source, dest)
        logging.info(f"Created symlink: {source} -> {dest}")
    else:
        logging.warning(f"Symlink already exists: {dest}")

def filter_links_by_utgs(flag, utg_file, input_file, output_file,logger):
    temp_file = f"temp_utg_file_{flag}.txt"
    try:
        command = f"cut -f1 {utg_file} > {temp_file}"
        execute_command(command, f"Failed to create temp file from {utg_file}", logger)
        
        command = (
            f"awk -F \"[,\\t ]\" 'NR==FNR{{lines[$1]; next}} ($1 in lines && $2 in lines)' "
            f"{temp_file} {input_file} > {output_file}"
        )
        execute_command(command, f"Failed to filter links for {output_file}", logger)
        
        command = f"rm {temp_file}"
        execute_command(command, f"Failed to remove temp file {temp_file}", logger)
        
        command = f"sed -i '1isource,target,links' {output_file}"
        execute_command(command, f"Failed to add header to {output_file}", logger)
    
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise

def adjust_r_and_cluster(initial_r, min_r, max_r, step, cluster_file, cluster_command, hap_number,logger):
    r = initial_r
    while min_r <= r <= max_r:
        logger.info(f"Running clustering with r={r}")
        execute_command(cluster_command.format(r=r), f"Failed to run clustering with r={r}",logger)
        
        # Check the number of clusters
        with open(cluster_file, 'r') as f:
            clusters = set(line.split(',')[0].strip() for line in f if line.strip())
        num_clusters = len(clusters)
        logger.info(f"Number of clusters: {num_clusters} (Target: {hap_number})")

        if num_clusters == hap_number:
            logger.info(f"Optimal r found: {r}")
            return r
        elif num_clusters > hap_number:
            max_r = r - step
        else:
            min_r = r + step

        r = (min_r + max_r) / 2

    raise ValueError(f"Failed to find optimal r within range {min_r}-{max_r}")

def process_chromosome(chr_num, args, pwd, partig_file,logger):
    try:

        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        os.makedirs(chr_dir, exist_ok=True)

        os.chdir(chr_dir)
        logging.info(f"Processing chromosome {chr_num} in {chr_dir}")

        # Create symlinks with ThreadPoolExecutor
        param_files = {
            "chr_cluster_file": f"../{args.chr_cluster_file}",
            "digraph_file": f"../{args.digraph_file}",
            "subgraph_file": f"../{args.subgraph_file}",
            "collapse_num_file": f"../{args.collapse_num_file}",
            "HiC_file": f"../{args.HiC_file}",
            "RE_file": f"../{args.RE_file}",
            "partig_file": f"../{partig_file}",
            "merge_partig_file": f"../merge.partig.csv"
        }

        merge_partig_file = "merge.partig.csv"
        
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            futures = []
            for name, filepath in param_files.items():
                dest = os.path.join(chr_dir, os.path.basename(filepath))
                futures.append(executor.submit(create_symlink, filepath, dest, logger))
            
            # Wait for all symlinks to be created
            for future in futures:
                future.result()


        utg_file = f"{args.output_prefix}.chr{chr_num}.utgs.txt"
        utg_rescue_file = f"{args.output_prefix}.chr{chr_num}.utgs.rescue.txt"
        create_symlink(f"../group{chr_num}.txt", utg_file, logger)
        create_symlink(f"../group{chr_num}_rescue.txt", utg_rescue_file, logger)

        # Process files with ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            # Process links file
            links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.csv"
            links_future = executor.submit(filter_links_by_utgs,str(chr_num)+"_links", utg_file, args.HiC_file, links_file, logger)
            
            # Process partig file
            partig_file = f"{args.output_prefix}.chr{chr_num}.partig.csv"
            partig_future = executor.submit(filter_links_by_utgs,str(chr_num)+"_partig", utg_file, merge_partig_file, partig_file, logger)
            
            # Wait for both operations to complete
            links_future.result()
            partig_future.result()


        filtered_links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.filterAllele.csv"
        command = (
            f"awk -F '[, \\t]' 'NR==FNR{{lines[$1\",\"$2];next}}!($1\",\"$2 in lines)' "
            f"{partig_file} {links_file} > {filtered_links_file}"
        )
        execute_command(command, f"Failed to filter allele links for {filtered_links_file}",logger)
        execute_command(f"sed '1isource,target,links' -i {filtered_links_file}", f"Failed to add header to {filtered_links_file}",logger)

        # Run louvain_nei.py
        script_path = os.path.abspath(sys.path[0])
        script_path_add = os.path.join(script_path, "louvain_nei.py")
        execute_command(
            f"python {script_path_add} -c {args.collapse_num_file} -chr {utg_file} "
            f"-l {filtered_links_file} -a {partig_file}",
            "Failed to run louvain_nei.py",logger
        )

        # Adjust r and run multilevel_cluster.py
        script_path_add = os.path.join(script_path, "multilevel_cluster.py")
        cluster_output = f"{args.output_prefix}.chr{chr_num}.cluster.txt"
        cluster_command = (
            f"python {script_path_add} -c louvain_nei.csv -o {cluster_output} -r {{r}}"
        )
        try:
            optimal_r = adjust_r_and_cluster(
                initial_r=1.0,
                min_r=0.01,
                max_r=5.0,
                step=0.01,
                cluster_file=cluster_output,
                cluster_command=cluster_command,
                hap_number=args.hap_number, logger=logger
            )
            logger.info(f"Optimal r for chromosome {chr_num}: {optimal_r}")
        except ValueError as e:
            logger.error(f"Failed to determine optimal r for chromosome {chr_num}: {e}")
            return


        if args.correct:
            script_path_add = os.path.join(script_path, "cluster_correct.py")
            execute_command(
                f"python {script_path_add} -cn {args.collapse_num_file} "
                f"-chr {utg_file} -l {filtered_links_file} -a {partig_file} -c {cluster_output} -n_hap {args.hap_number}",
                "Failed to run cluster_correct.py",logger
            )
            cluster_file = "cor.cluster.txt"
        
        else:
            cluster_file = cluster_output

        # Run louvain_reassign_allele.py
        script_path_add = os.path.join(script_path, "louvain_reassign_allele.py")
        execute_command(
            f"python {script_path_add} -c {args.collapse_num_file} "
            f"-chr {utg_rescue_file} -l {links_file} -r {args.RE_file} -a {merge_partig_file} "
            f"--clusters {cluster_file} --isolated_threshold {args.isolated_threshold}",
            "Failed to run louvain_reassign_allele.py",logger
        )

        
        return f"Successfully processed chromosome {chr_num}"
    except Exception as e:
        logger.error(f"Error processing chromosome {chr_num}: {e}")
        return f"Failed to process chromosome {chr_num}: {e}"

def run_partig(fa_file, partig_k, partig_w, partig_c, partig_m, output_prefix,logger):
    execute_command(f"samtools faidx {fa_file}", "Failed to indexing the assembly FASTA",logger)
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
    script_path_add = os.path.join(script_path, "trans.partig.py")
    execute_command(
        f"python {script_path_add} -fai {fa_file}.fai -p {output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt -o {output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.csv", "Converting Partig output to CSV",logger)

def log_start(logger, script_name: str, version: str, args: argparse.Namespace):
    """Log the start of the program."""
    logger.info(f"Program started, {script_name} version: {version}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")



def main():
    logger = setup_logging('cluster_hap.log')
    args = parse_arguments()
    pwd = os.getcwd()

    log_start(logger, "cluster_hap.py", "1.0.0", args)

    # Step 1: Run partig
    run_partig(args.fa_file, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix,logger)
    convert_partig_output(args.fa_file, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix,logger)
    partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{args.partig_m}.csv"

    # Step 2: Run cluster2group.py
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "cluster2group.py")
    execute_command(
        f"python {script_path_add} -c {args.chr_cluster_file} -r {args.RE_file}",
        "Failed to run cluster2group.py",logger
    )
    execute_command(
        f"python {script_path_add} -c {args.chr_cluster_rescue_file} -r {args.RE_file} -m rescue",
        "Failed to run cluster2group.py",logger
    )

    script_path_add = os.path.join(script_path, "filter_expand_partig.py")
    execute_command(
            f"python {script_path_add} -d {args.digraph_file} "
            f"-r {args.RE_file} -s {args.subgraph_file} -p {partig_file}",
            "Failed to run filter_expand_partig.py",logger
    )

    # Step 3: Process chromosomes in parallel
    with Pool(processes=args.process_num) as pool:

        # Create partial function with fixed arguments
        process_chr = functools.partial(process_chromosome, args=args, pwd=pwd, partig_file=partig_file,logger=logger)
        
        # Process chromosomes in parallel
        results = pool.map(process_chr, range(1, args.chr_number + 1))
        
        # Log results
        for result in results:
            logging.info(result)



if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")