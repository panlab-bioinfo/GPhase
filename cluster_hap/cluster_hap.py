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

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def execute_command(command, error_message):
    try:
        subprocess.run( command, shell=True, check=True, executable='/bin/bash')
        logging.info(f"Command executed successfully: {command}")
    except subprocess.CalledProcessError:
        logging.error(error_message)
        raise

def create_symlink(source, dest):
    if not os.path.exists(dest):
        os.symlink(source, dest)
        logging.info(f"Created symlink: {source} -> {dest}")
    else:
        logging.warning(f"Symlink already exists: {dest}")

def filter_links_by_utgs(utg_file, input_file, output_file):
    command = (
        f"awk -F \"[,\\t ]\" \"NR==FNR{{lines[\$1];next}}(\$1 in lines && \$2 in lines)\" "
        f"<(cut -f1 {utg_file}) {input_file} > {output_file}"
    )
    execute_command(command, f"Failed to filter links for {output_file}")
    execute_command(f"sed '1isource,target,links' -i {output_file}", f"Failed to add header to {output_file}")

def adjust_r_and_cluster(initial_r, min_r, max_r, step, cluster_file, cluster_command, hap_number):
    r = initial_r
    while min_r <= r <= max_r:
        logging.info(f"Running clustering with r={r}")
        execute_command(cluster_command.format(r=r), f"Failed to run clustering with r={r}")
        
        # Check the number of clusters
        with open(cluster_file, 'r') as f:
            clusters = set(line.split(',')[0].strip() for line in f if line.strip())
        num_clusters = len(clusters)
        logging.info(f"Number of clusters: {num_clusters} (Target: {hap_number})")

        if num_clusters == hap_number:
            logging.info(f"Optimal r found: {r}")
            return r
        elif num_clusters > hap_number:
            max_r = r - step
        else:
            min_r = r + step

        r = (min_r + max_r) / 2

    raise ValueError(f"Failed to find optimal r within range {min_r}-{max_r}")

def process_chromosome(chr_num, args, pwd, partig_file):
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
                futures.append(executor.submit(create_symlink, filepath, dest))
            
            # Wait for all symlinks to be created
            for future in futures:
                future.result()


        utg_file = f"{args.output_prefix}.chr{chr_num}.utgs.txt"
        utg_rescue_file = f"{args.output_prefix}.chr{chr_num}.utgs.rescue.txt"
        create_symlink(f"../group{chr_num}.txt", utg_file)
        create_symlink(f"../group{chr_num}_rescue.txt", utg_rescue_file)

        # Process files with ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            # Process links file
            links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.csv"
            links_future = executor.submit(filter_links_by_utgs, utg_file, args.HiC_file, links_file)
            
            # Process partig file
            partig_file = f"{args.output_prefix}.chr{chr_num}.partig.csv"
            partig_future = executor.submit(filter_links_by_utgs, utg_file, merge_partig_file, partig_file)
            
            # Wait for both operations to complete
            links_future.result()
            partig_future.result()


        filtered_links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.filterAllele.csv"
        command = (
            f"awk -F '[, \\t]' 'NR==FNR{{lines[$1\",\"$2];next}}!($1\",\"$2 in lines)' "
            f"{partig_file} {links_file} > {filtered_links_file}"
        )
        execute_command(command, f"Failed to filter allele links for {filtered_links_file}")
        execute_command(f"sed '1isource,target,links' -i {filtered_links_file}", f"Failed to add header to {filtered_links_file}")

        # Run louvain_nei.py
        script_path = os.path.abspath(sys.path[0])
        script_path_add = os.path.join(script_path, "louvain_nei.py")
        execute_command(
            f"python {script_path_add} -c {args.collapse_num_file} -chr {utg_file} "
            f"-l {filtered_links_file} -a {partig_file}",
            "Failed to run louvain_nei.py"
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
                hap_number=args.hap_number
            )
            logging.info(f"Optimal r for chromosome {chr_num}: {optimal_r}")
        except ValueError as e:
            logging.error(f"Failed to determine optimal r for chromosome {chr_num}: {e}")
            return


        if args.correct:
            script_path_add = os.path.join(script_path, "cluster_correct.py")
            execute_command(
                f"python {script_path_add} -cn {args.collapse_num_file} "
                f"-chr {utg_file} -l {filtered_links_file} -a {partig_file} -c {cluster_output} -n_hap {args.hap_number}",
                "Failed to run cluster_correct.py"
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
            "Failed to run louvain_reassign_allele.py"
        )

        
        return f"Successfully processed chromosome {chr_num}"
    except Exception as e:
        logging.error(f"Error processing chromosome {chr_num}: {e}")
        return f"Failed to process chromosome {chr_num}: {e}"

def run_partig(asm_fa, partig_k, partig_w, partig_c, partig_m, output_prefix):
    execute_command(f"samtools faidx {asm_fa}, Indexing the assembly FASTA")
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
            logging.info(f"Starting: Running Partig")
            logging.info(f"Command: {' '.join(command)}")
            subprocess.run(command, check=True, stdout=outfile, stderr=subprocess.PIPE, text=True)
            logging.info(f"Completed: Running Partig\n")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in: Running Partig\nCommand failed: {' '.join(command)}")
        logging.error(f"Error output: {e.stderr}")
        return False


def convert_partig_output(asm_fa, partig_k, partig_w, partig_c, partig_m, output_prefix):
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "trans.partig.py")
    return execute_command(
        f"python {script_path_add} -fai {asm_fa}.fai -p {output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.txt -o {output_prefix}.partig.{partig_k}_{partig_w}_{partig_c}_{partig_m}.csv", "Converting Partig output to CSV")

def main(args):
    setup_logging()
    pwd = os.getcwd()

    # Log configuration
    logging.info(f"Running with {args.process_num} processes and {args.thread_num} threads per process")

    # Step 1: Run partig
    run_partig(args.asm_fa, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix)
    convert_partig_output(args.asm_fa, args.partig_k, args.partig_w, args.partig_c, args.partig_m, args.output_prefix)
    partig_file = f"{args.output_prefix}.partig.{args.partig_k}_{args.partig_w}_{args.partig_c}_{args.partig_m}.csv"

    # Step 2: Run cluster2group.py
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "cluster2group.py")
    execute_command(
        f"python {script_path_add} -c {args.chr_cluster_file} -r {args.RE_file}",
        "Failed to run cluster2group.py"
    )
    execute_command(
        f"python {script_path_add} -c {args.chr_cluster_rescue_file} -r {args.RE_file} -m rescue",
        "Failed to run cluster2group.py"
    )

    script_path_add = os.path.join(script_path, "filter_expand_partig.py")
    execute_command(
            f"python {script_path_add} -d {args.digraph_file} "
            f"-r {args.RE_file} -s {args.subgraph_file} -p {partig_file}",
            "Failed to run filter_expand_partig.py"
    )

    # Step 3: Process chromosomes in parallel
    with Pool(processes=args.process_num) as pool:

        # Create partial function with fixed arguments
        process_chr = functools.partial(process_chromosome, args=args, pwd=pwd, partig_file=partig_file)
        
        # Process chromosomes in parallel
        results = pool.map(process_chr, range(1, args.chr_number + 1))
        
        # Log results
        for result in results:
            logging.info(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallel Hi-C data processing pipeline")

    # Required arguments
    parser.add_argument("-a", "--asm_fa", required=True, help="Path to the assembly FASTA file.")
    parser.add_argument("-c", "--chr_cluster_file", required=True, help="Chromosome cluster file")
    parser.add_argument("-cr", "--chr_cluster_rescue_file", required=True, help="Chromosome cluster file")
    parser.add_argument("-d", "--digraph_file", required=True, help="Directed graph file")
    parser.add_argument("-s", "--subgraph_file", required=True, help="Subgraph file")
    parser.add_argument("-cn", "--collapse_num_file", required=True, help="Collapse number file")
    parser.add_argument("-l", "--HiC_file", required=True, help="Hi-C file")
    parser.add_argument("-r", "--RE_file", required=True, help="Restriction enzyme file")
    parser.add_argument('--correct', action='store_true', help='correct the cluster')
    parser.add_argument('--isolated_threshold', default=5,help='<int>Detect whether the intensity of the hic signal is an outlier')
    parser.add_argument("-pk", "--partig_k", type=int, default=17, help="K-mer size for Partig. Default: 17.")
    parser.add_argument("-pw", "--partig_w", type=int, default=17, help="Minimizer window size for Partig. Default: 17.")
    parser.add_argument("-pc", "--partig_c", type=int, default=60, help="Max occurrance for Partig. Default: 60.")
    parser.add_argument("-pm", "--partig_m", type=float, default=0.6, help="Mini k-mer similarity for Partig. Default: 0.6.")
    parser.add_argument("-op", "--output_prefix", required=True, help="Output prefix")
    parser.add_argument("-n_chr", "--chr_number", type=int, required=True, help="Number of chromosomes")
    parser.add_argument("-n_hap", "--hap_number", type=int, required=True, help="Number of haplotypes")

    # Optional arguments for parallel processing
    parser.add_argument("--process_num", type=int, default=10,
                       help=f"Number of parallel processes (default: 10)")
    parser.add_argument("--thread_num", type=int,default=10,
                       help=f"Number of threads per process (default: 10)")

    args = parser.parse_args()
    argcomplete.autocomplete(parser)


    try:
        main(args)
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")