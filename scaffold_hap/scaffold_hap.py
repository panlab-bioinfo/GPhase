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
from get_subgraph_scaffold_v2 import Get_subgraph_scaffold_v2
from get_data_HapHiC_sort import Get_data_HapHiC_sort

def setup_logging(log_file: str = "scaffold.log") -> logging.Logger:
    """Configure logging to both file and console."""
    logger = logging.getLogger('genome_scaffolding')
    logger.setLevel(logging.INFO)
    
    # File handler
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
    performance_group .add_argument("-t", "--thread_number", metavar='\b', type=int, default=12, help="Number of parallel processes, default: 12")

    yahs_group  = parser.add_argument_group('>>> Parameters for YaHS')
    yahs_group .add_argument("--no_contig_ec", action='store_true', help="do not do contig error correction")
    yahs_group .add_argument("--no_scaffold_ec", action='store_true', help="do not do scaffold error correction")

    haphic_group  = parser.add_argument_group('>>> Parameters for HapHiC sort')
    haphic_group .add_argument("--min_len", metavar='\b',type=int, default=100, help="minimum scaffold length(kb), default: 100")
    haphic_group .add_argument("--mutprob", metavar='\b', type=float, default=0.6, help="mutation probability in the genetic algorithm, default: 0.6")
    haphic_group .add_argument("--ngen", metavar='\b', type=int, default=20000, help="number of generations for convergence, default: 20000")
    haphic_group .add_argument("--npop", metavar='\b', type=int, default=200, help="mopulation size, default: 200")
    haphic_group .add_argument("--processes", metavar='\b', type=int, default=32, help="processes for fast sorting and ALLHiC optimization, default: 32")

    
    return parser.parse_args()

def run_command(cmd: List[str], logger, cwd: Optional[str] = None, shell: bool = False, stdout=None) -> bool:
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

#
def process_one_Contig(file_path):
    
    dir_path = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    output_path = os.path.dirname(dir_path)
    new_file_name = file_name.replace(".txt", ".tour")
    final_output_path = os.path.join(output_path, new_file_name)
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

        if len(lines) == 2:
            contig_id = lines[1].strip().split()[0]
            with open(final_output_path, 'w', encoding='utf-8') as out:
                out.write(">INIT\n")
                out.write(contig_id + "+\n")
                out.write(">FLIPWHOLE1\n")
                out.write(contig_id + "+\n")
                out.write(">FLIPONE1\n")
                out.write(contig_id + "+\n")
            os.symlink(final_output_path, "../" + new_file_name)  
            return new_file_name
        else:
            os.symlink(file_path, file_name)  
            return None


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

        # Run get subgraph sort result  
        script_path = os.path.abspath(sys.path[0])
        logger.info(f"Chromosome {chr_num} hap {hap_num} run Get_subgraph_scaffold_v2...")
        logger.info(f"Chromosome {chr_num} hap {hap_num} run Get_subgraph_scaffold_v2: {args.gfa_file} {args.RE_file} {args.digraph_file} group{hap_num}.txt {args.subgraph_file}")
        Get_subgraph_scaffold_v2(args.gfa_file, args.RE_file, args.digraph_file,f"group{hap_num}.txt", args.subgraph_file)

        logger.info(f"Chromosome {chr_num} hap {hap_num} run Get_subgraph_scaffold_v2 done...")

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

        if args.map_file_type == "pa5":

            tmp_file = f"tmp_{hap_num}.txt"
            cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
            if not run_command([cut_command], shell=True,logger=logger):
                return False

            awk_cmd = f"""awk 'NR==FNR{{lines[$1];next}}{{if(NF < 8 && $2 in lines){{print $0;next}}if($2 in lines && $4 in lines){{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5}}}}' {tmp_file} {args.map_file} > chr{chr_num}g{hap_num}.pairs"""
            if not run_command([awk_cmd], shell=True, logger=logger):
                return False

            cleanup_command = f"rm {tmp_file}"
            if not run_command([cleanup_command], shell=True, logger=logger):
                return False

            # Run YAHS
            contig_ec = "--no-contig-ec" if args.no_contig_ec else ""
            scaffold_ec = "--no-scaffold-ec" if args.no_scaffold_ec else ""
            if not run_command([
                "yahs",
                f"chr{chr_num}g{hap_num}.fa",
                f"chr{chr_num}g{hap_num}.pairs",
                "-a", "subgraphGroup.agp",
                "-o", "subgraph_yahs",
                "--file-type", "pa5", f"{contig_ec}", f"{scaffold_ec}"
            ], logger=logger):
                return False

        # Update scaffold names
        for file_pattern, sed_cmd in [
            ("subgraph_yahs_scaffolds_final.agp", f"s/^/chr{chr_num}g{hap_num}_/g"),
            ("subgraph_yahs_scaffolds_final.fa", f"s/>/>chr{chr_num}g{hap_num}_/g"),
            ("subgraphGroup.agp", f"s/^/chr{chr_num}g{hap_num}_/g"),
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
                    f"chr{chr_num}g{hap_num}/subgraphGroup.agp", \
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
        if args.map_file_type == "pa5":
            flag = Get_data_HapHiC_sort(f"chr{chr_num}g{hap_num}.pairs", "pa5",  "subgraph_yahs_scaffolds_final.agp", args.RE_file, f"{args.output_prefix}.chr{chr_num}g{hap_num}", args.min_len)
            if not flag:
                logger.error(f"Error : Chr{chr_num}g{hap_num} get data for HapHiC -> {str(e)}")
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
    one_contigs_list = list()
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
                re = process_one_Contig(re_file)
                if re:
                    one_contigs_list.append(re)
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


        # HapHiC sort & build
        script_path_add = os.path.join(script_path, "../src/HapHiC/haphic")
        script_content_1 = f"""
            cd {haphic_dir}
            {script_path_add} sort {args.output_prefix}.merge.fa {args.output_prefix}.merge.HT.pkl split_clms/ groups_REs/* --skip_fast_sort --mutprob {args.mutprob} --ngen {args.ngen} --npop {args.npop}  --processes {args.processes}
        """
        script_content_2 = f"""
            cd {haphic_dir}
            {script_path_add} build {args.output_prefix}.merge.fa {args.output_prefix}.merge.fa {args.output_prefix}.merge.fa final_tours/*tour
            seqkit sort scaffolds.fa > scaffolds.sort.fa
            samtools faidx scaffolds.sort.fa
        """

        with tempfile.NamedTemporaryFile(delete=False) as temp_script:
            temp_script.write(script_content_1.encode())
            temp_script_path = temp_script.name

        cmd = ["bash", temp_script_path]
        if not run_command(cmd, logger=logger):
            return False
        logger.info("HapHiC sort completed.")

        # copy ont Contig scaffold
        os.chdir(os.path.join(haphic_dir, "final_tours"))
        for group in one_contigs_list:
            os.symlink("../" + group, group)  

        os.chdir(haphic_dir)
        with tempfile.NamedTemporaryFile(delete=False) as temp_script:
            temp_script.write(script_content_2.encode())
            temp_script_path = temp_script.name

        cmd = ["bash", temp_script_path]
        if not run_command(cmd, logger=logger):
            return False
        logger.info("HapHiC build completed.")

        # get final agp (using agptools)
        os.chdir(haphic_dir)
        os.makedirs(os.path.join(haphic_dir, "final_agp"), exist_ok=True)
        os.chdir(os.path.join(haphic_dir, "final_agp"))

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
                subprocess.run(f"agptools join chr{i}g{j}.joins.txt <(grep 'chr{i}g{j}' ../{args.output_prefix}.merge.agp) -n 100 -e proximity_ligation | "
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
    args = parse_arguments()
    logger = setup_logging('scaffold_hap.log') 
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
