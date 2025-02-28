#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple, Optional, Union

def setup_logging(log_file: str = "scaffold.log") -> logging.Logger:
    """Configure logging to both file and console."""
    logger = logging.getLogger('genome_scaffolding')
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
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(prog='scaffold_hap')

    base_group  = parser.add_argument_group('>>> Parameters for basic data')
    base_group.add_argument("-f", "--fa_file", metavar='\b', required=True,help="Fasta file")
    base_group.add_argument("-r", "--RE_file", metavar='\b', required=True,help="Path to RE file")

    cluster_group  = parser.add_argument_group('>>> Parameters of the file for the haplotype clustering results')
    cluster_group.add_argument("-CHP", "--cluster_hap_path",metavar='\b', required=True, help="Path to cluster hap directory")

    graph_group  = parser.add_argument_group('>>> Parameters related to subgraphs')
    graph_group.add_argument("-s", "--subgraph_file",metavar='\b', required=True, help="Subgraph resulting from the GFA split")
    graph_group.add_argument("-g", "--gfa_file", metavar='\b', required=True,help="GFA file")
    graph_group.add_argument("-d", "--digraph_file", metavar='\b', required=True,help="Directed graph resulting from GFA conversion")

    hic_group  = parser.add_argument_group('>>> Parameters for HiC data alignment')
    hic_group.add_argument("--map_file_type",metavar='\b', default="pa5", choices=['bam', 'pa5'], type=str.lower,help="HiC mapping file type")
    hic_group.add_argument("-m", "--map_file", metavar='\b', required=True,help="HiC mapping file")
    hic_group.add_argument("-l", "--HiC_file", metavar='\b', required=True,help="HiC links file")

    genome_group  = parser.add_argument_group('>>> Parameters of chromosome and haplotype numbers')
    genome_group.add_argument("-n_chr", "--chr_number", metavar='\b', required=True,type=int, help="Number of chromosomes")
    genome_group.add_argument("-n_hap", "--hap_number", metavar='\b', required=True,type=int, help="Number of haplotypes")

    output_group  = parser.add_argument_group('>>> Parameter for the prefix of the result file')
    output_group.add_argument("-op", "--output_prefix", metavar='\b', required=True, help="output_group")

    performance_group  = parser.add_argument_group('>>> Parameters for performance')
    performance_group .add_argument("-t", "--thread_number", metavar='\b', type=int, default=1, help="Number of parallel processes")
    
    return parser.parse_args()

def run_command(cmd: List[str], logger, cwd: Optional[str] = None, shell: bool = False) -> bool:
    """Execute a shell command."""
    try:
        logger.info(f"Running command: {cmd}")
        subprocess.run(cmd, cwd=cwd, shell=shell, check=True)
        logger.info(f"Running command successfully: {cmd}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Error: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error executing command: {str(e)}")
        return False

def process_chromosome(pwd: str, chr_num: int, args: argparse.Namespace,logger) -> bool:
    """Process a single chromosome directory."""
    logger.info(f"Processing chromosome {chr_num}...")
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "cluster2group.py")
    try:
        # Create chromosome directory
        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        os.makedirs(chr_dir, exist_ok=True)
        os.chdir(chr_dir)

        # Create symbolic links
        os.symlink(os.path.join("../", f"{args.RE_file}"), f"{args.output_prefix}.RE_counts.txt")
        os.symlink(
            os.path.join("../", args.cluster_hap_path, f"chr{chr_num}", f"{args.output_prefix}.reassign.cluster.txt"),
            f"{args.output_prefix}.reassign.cluster.txt"
        )

        # Run cluster2group.sh
        cmd = ["python", script_path_add, "-c", f"{args.output_prefix}.reassign.cluster.txt", "-r", f"{args.output_prefix}.RE_counts.txt"]
        if not run_command(cmd, logger):
            return False

        logger.info(f"Chromosome {chr_num} processing completed.")
        return True
    except Exception as e:
        logger.error(f"Error processing chromosome {chr_num}: {str(e)}")
        return False

def process_haplotype(pwd: str, chr_num: int, hap_num: int, args: argparse.Namespace,logger) -> bool:
    """Process a single haplotype."""
    logger.info(f"Processing haplotype {hap_num} of chromosome {chr_num}...")
    try:
        # Change to chromosome directory
        os.chdir(pwd)
        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        hap_dir = os.path.join(chr_dir, f"chr{chr_num}g{hap_num}")
        os.makedirs(hap_dir, exist_ok=True)
        os.chdir(hap_dir)

        # Create symbolic links
        src_files = [
            f"chr{chr_num}/group{hap_num}.txt",
            args.subgraph_file,
            args.gfa_file,
            args.digraph_file,
            args.RE_file,
            args.map_file,
            args.HiC_file,
            args.fa_file
        ]
        for src in src_files:
            try:
                os.symlink(os.path.join(pwd, src), os.path.basename(src))
            except FileExistsError:
                pass

        # Run scaffold.subgraph.py
        script_path = os.path.abspath(sys.path[0])
        script_path_add = os.path.join(script_path, "get_subgraph_scaffold.py")

        cmd = [
            "python", script_path_add,
            "-c", f"group{hap_num}.txt",
            "-subgraph", args.subgraph_file,
            "-graph", args.gfa_file,
            "-digraph", args.digraph_file,
            "-r", args.RE_file,
            "-l", args.HiC_file
        ]

        if not run_command(cmd,logger):
            return False

        # split asm.fa
        tmp_file = f"tmp_{hap_num}.txt"
        cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
        if not run_command([cut_command], shell=True, logger=logger):
            return False
        # seqkit 
        seqkit_cmd = f"seqkit grep -f {tmp_file} {args.fa_file} > chr{chr_num}g{hap_num}.fa"
        if not run_command([seqkit_cmd], shell=True, logger=logger):
            return False
        cleanup_command = f"rm {tmp_file}"
        if not run_command([cleanup_command], shell=True, logger=logger):
            return False

        # Create fasta index
        if not run_command(["samtools", "faidx", f"chr{chr_num}g{hap_num}.fa"], logger=logger):
            return False

        # Process bam or pairs file
        if args.map_file_type == "bam":
            command = f"cut -f1 group{hap_num}.txt | xargs -I T echo T | tr '\n' ' ' "
            cut_output = subprocess.check_output(command, shell=True).decode('utf-8').strip()

            # Use it in the samtools command
            samtools_cmd = f"samtools view -b {args.map_file} {cut_output} -o chr{chr_num}g{hap_num}.bam"
            if not run_command([samtools_cmd], shell=True, logger=logger):
                return False

            # Run YAHS
            if not run_command([
                "yahs",
                f"chr{chr_num}g{hap_num}.fa",
                f"chr{chr_num}g{hap_num}.bam",
                "-a", "subgraphGroup.agp",
                "-o", "subgraph_yahs"
            ], logger=logger):
                return False

        if args.map_file_type == "pa5":

            tmp_file = f"tmp_{hap_num}.txt"
            cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
            if not run_command([cut_command], shell=True,logger=logger):
                return False

            # 使用 awk 命令
            awk_cmd = f"""awk 'NR==FNR{{lines[$1];next}}{{if(NF < 8 && $2 in lines){{print $0;next}}if($2 in lines && $4 in lines){{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5}}}}' {tmp_file} {args.map_file} > chr{chr_num}g{hap_num}.pairs"""
            if not run_command([awk_cmd], shell=True, logger=logger):
                return False

            cleanup_command = f"rm {tmp_file}"
            if not run_command([cleanup_command], shell=True, logger=logger):
                return False

            # Run YAHS
            if not run_command([
                "yahs",
                f"chr{chr_num}g{hap_num}.fa",
                f"chr{chr_num}g{hap_num}.pairs",
                "-a", "subgraphGroup.agp",
                "-o", "subgraph_yahs",
                "--file-type", "pa5"
            ], logger=logger):
                return False

        # Update scaffold names
        for file_pattern, sed_cmd in [
            ("subgraph_yahs_scaffolds_final.agp", f"s/^/chr{chr_num}g{hap_num}_/g"),
            ("subgraph_yahs_scaffolds_final.fa", f"s/>/>chr{chr_num}g{hap_num}_/g")
        ]:
            if not run_command(["sed", "-i", sed_cmd, file_pattern], logger=logger):
                return False

        # Set up scaffold_HapHiC_sort directory
        scaffold_dir = os.path.join(hap_dir, "scaffold_HapHiC_sort")
        os.makedirs(scaffold_dir, exist_ok=True)
        os.chdir(scaffold_dir)

        # Create links in scaffold directory
        links_list = [f"chr{chr_num}g{hap_num}/subgraph_yahs.bin", \
                    f"chr{chr_num}g{hap_num}/chr{chr_num}g{hap_num}.fa", \
                    f"chr{chr_num}g{hap_num}/subgraph_yahs_scaffolds_final.agp", \
                    args.RE_file]

        if args.map_file_type == "bam":
            links_list.append(f"chr{chr_num}g{hap_num}/chr{chr_num}g{hap_num}.bam")
        elif args.map_file_type == "pa5":
            links_list.append(f"chr{chr_num}g{hap_num}/chr{chr_num}g{hap_num}.pairs")

        for src in links_list:
            try:
                os.symlink(os.path.join(pwd, f"chr{chr_num}", src), os.path.basename(src))
            except FileExistsError:
                pass

        # Run get_data_HapHiC_sort.py
        script_path_add = os.path.join(script_path, "get_data_HapHiC_sort.py")
        if args.map_file_type == "bam":
            if not run_command([
                "python", f"{script_path_add}",
                "-m", f"chr{chr_num}g{hap_num}.bam",
                "-a", "subgraph_yahs_scaffolds_final.agp",
                "-r", args.RE_file,
                "-o", f"{args.output_prefix}.chr{chr_num}g{hap_num}"
            ],logger=logger):
                return False
        if args.map_file_type == "pa5":
            if not run_command([
                "python", f"{script_path_add}",
                "-m", f"chr{chr_num}g{hap_num}.pairs",
                "-a", "subgraph_yahs_scaffolds_final.agp",
                "-r", args.RE_file,
                "-o", f"{args.output_prefix}.chr{chr_num}g{hap_num}"
            ], logger=logger):
                return False

        logger.info(f"Haplotype {hap_num} of chromosome {chr_num} processing completed.")
        return True
    except Exception as e:
        logger.error(f"Error processing chr{chr_num}g{hap_num}: {str(e)}")
        return False

def perform_final_merge(pwd: str, args: argparse.Namespace,logger) -> bool:
    """Perform final merging steps."""
    logger.info(f"Starting final merge steps...")
    script_path = os.path.abspath(sys.path[0])
    try:
        os.chdir(pwd)
        haphic_dir = os.path.join(pwd, "HapHiC_sort")
        os.makedirs(os.path.join(haphic_dir, "split_clms"), exist_ok=True)
        os.makedirs(os.path.join(haphic_dir, "groups_REs"), exist_ok=True)
        os.chdir(os.path.join(haphic_dir, "split_clms"))

        # Link CLM files
        clm_files = glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort", f"{args.output_prefix}.*chr*g*.clm"))
        for clm_file in clm_files:
            try:
                os.symlink(clm_file, os.path.basename(clm_file))
            except FileExistsError:
                pass

        # Link RE files
        os.chdir(os.path.join(haphic_dir, "groups_REs"))
        re_files = glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort/", f"{args.output_prefix}.chr*g*scaffold*txt"))
        for re_file in re_files:
            try:
                os.symlink(re_file, os.path.basename(re_file))
            except FileExistsError:
                pass

        os.chdir(haphic_dir)

        # Merge AGP files
        with open(f"{args.output_prefix}.merge.agp", 'w') as outfile:
            for agp_file in glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort/subgraph_yahs_scaffolds_final.agp")):
                with open(agp_file) as infile:
                    outfile.write(infile.read())


        # Merge HT files
        with open(f"{args.output_prefix}.merge.HT.pkl", 'wb') as outfile:
            for ht_file in glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort/*HT*")):
                with open(ht_file, 'rb') as infile:
                    outfile.write(infile.read())

        # Merge scaffold files
        with open(f"{args.output_prefix}.merge.fa", 'w') as outfile:
            for scaffold_file in glob.glob(os.path.join(pwd, "chr*/chr*/subgraph_yahs_scaffolds_final.fa")):
                with open(scaffold_file) as infile:
                    outfile.write(infile.read())

        # # Merge pairs files
        # with open(f"{args.output_prefix}.merge.map.bin", 'wb') as outfile:
        #     for pairs_file in glob.glob(os.path.join(pwd, "chr*/chr*/subgraph_yahs.bin")):
        #         with open(pairs_file, 'rb') as infile:
        #             outfile.write(infile.read())

        # Create and execute final commands script
        script_path_add = os.path.join(script_path, "../src/HapHiC/haphic")
        script_content = f"""
            cd {haphic_dir}
            {script_path_add} sort {args.output_prefix}.merge.fa {args.output_prefix}.merge.HT.pkl split_clms/ groups_REs/* 
            {script_path_add} build {args.output_prefix}.merge.fa {args.output_prefix}.merge.fa {args.output_prefix}.merge.map.bin final_tours/*tour
            seqkit sort scaffolds.fa > scaffolds.sort.fa
            samtools faidx scaffolds.sort.fa

        """

        with tempfile.NamedTemporaryFile(delete=False) as temp_script:
            temp_script.write(script_content.encode())
            temp_script_path = temp_script.name

        cmd = ["bash", temp_script_path]
        if not run_command(cmd, logger=logger):
            return False
        logger.info("HapHiC sort completed.")

        # get final agp (using agptools)
        os.chdir(haphic_dir)
        os.makedirs(os.path.join(haphic_dir, "final_agp"), exist_ok=True)
        os.chdir(os.path.join(haphic_dir, "final_agp"))
        # Loop through i and j values
        for i in range(1, int(args.chr_number) + 1):
            for j in range(1, int(args.hap_number) + 1):
                # Step 1: Extract rows from scaffolds.agp
                subprocess.run(f"grep chr{i}g{j} ../scaffolds.agp > chr{i}g{j}.agp", shell=True)

                # Step 2: Process and convert the data from .agp file into .joins.txt
                with open(f"chr{i}g{j}.agp", 'r') as infile, open(f"chr{i}g{j}.joins.txt", 'w') as outfile:

                    first_scaffold = True
                    for line in infile:
                        cols = line.strip().split()

                        if first_scaffold:
                            first_scaffold = False
                        elif cols[8] == "+" or cols[8] == "-":
                            outfile.write(",")
                           
                        if cols[8] == "+":
                            outfile.write(f"{cols[5]}")
                        elif cols[8] == "-":
                            outfile.write(f"-{cols[5]}")

                # Step 4: Join using agptools
                subprocess.run(f"agptools join chr{i}g{j}.joins.txt <(grep 'chr{i}g{j}' ../{args.output_prefix}.merge.agp) | "
                            f"awk -v i={i} -v j={j} 'BEGIN{{OFS=\"\\t\"}}{{$1=\"Chr\"i\"g\"j; print $0}}' >> {args.output_prefix}.final.agp", shell=True)
        
        return True
    except Exception as e:
        logger.error(f"Error during final merge: {str(e)}")
        return False

def log_start(logger: logging.Logger, script_name: str, version: str, args: argparse.Namespace):
    """Log the start of the program."""
    logger.info(f"Program started, {script_name} version: {version}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")

def main():
    logger = setup_logging('scaffold_hap.log')  # Set up logging
    args = parse_arguments()
    pwd = os.getcwd()
    
    # Log program start
    log_start(logger, "scaffold_hap.py", "1.0.0", args)
    
    # Process chromosomes
    for i in range(1, args.chr_number + 1):
        logger.info(f"Processing chromosome {i}...")
        if not process_chromosome(pwd, i, args ,logger):
            logger.error(f"Error processing chromosome {i}")
            sys.exit(1)
    
    # Process haplotypes in parallel
    with ProcessPoolExecutor(max_workers=args.thread_number) as executor:
        futures = []
        for i in range(1, args.chr_number + 1):
            for j in range(1, args.hap_number + 1):
                futures.append(executor.submit(process_haplotype, pwd, i, j, args,logger))
        
        for future in as_completed(futures):
            if not future.result():
                logger.error("Error processing haplotype.")
                sys.exit(1)

    # Final merge
    logger.info("Performing final merge steps...")
    if not perform_final_merge(pwd, args,logger):
        logger.error("Error in final merge steps")
        sys.exit(1)
    
    logger.info("Program completed successfully.")

if __name__ == "__main__":
    main()
