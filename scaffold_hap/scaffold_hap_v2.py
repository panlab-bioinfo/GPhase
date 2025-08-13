#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import subprocess
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Optional, Union
from trans_pairs import Trans_pairs
from get_subgraph_scaffold import Get_subgraph_scaffold
from get_data_HapHiC_sort import Get_data_HapHiC_sort
from rescue_base_graph import Rescue_base_graph

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'cluster_chr')))
from get_RE import Get_RE


def setup_logging(log_file: str = "scaffold.log") -> logging.Logger:
    """Configure logging to both file and console."""
    logger = logging.getLogger('genome_scaffolding')
    logger.setLevel(logging.INFO)

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
    base_group.add_argument("-e", "--enzyme_site", metavar='\b', default="GATC",help="Restriction enzyme file.")

    cluster_group  = parser.add_argument_group('>>> Parameters of the file for the haplotype clustering results')
    cluster_group.add_argument("-CHP", "--cluster_hap_path",metavar='\b', required=True, help="Path to cluster hap directory")

    graph_group  = parser.add_argument_group('>>> Parameters related to subgraphs')
    graph_group.add_argument("-s", "--subgraph_file",metavar='\b', required=True, help="Subgraph resulting from the GFA split")
    graph_group.add_argument("-g", "--gfa_file", metavar='\b', required=True,help="GFA file")
    graph_group.add_argument("-d", "--digraph_file", metavar='\b', required=True,help="Directed graph resulting from GFA conversion")

    hic_group  = parser.add_argument_group('>>> Parameters for HiC data alignment')
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
    haphic_group .add_argument("--mutprob", metavar='\b', type=float, default=0.3, help="mutation probability in the genetic algorithm, default: 0.3")
    haphic_group .add_argument("--ngen", metavar='\b', type=int, default=5000, help="number of generations for convergence, default: 5000")
    haphic_group .add_argument("--npop", metavar='\b', type=int, default=100, help="mopulation size, default: 100")
    haphic_group .add_argument("--processes", metavar='\b', type=int, default=32, help="processes for fast sorting and ALLHiC optimization, default: 32")

    
    return parser.parse_args()

def trans_agp(agp_file, original_agp, output_file):

    scaffold_map = defaultdict(list)

    with open(agp_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            scaffold = fields[0]
            if len(fields) >= 9 and "scaffold" in fields[5]:
                target = fields[5]
                ori = fields[8]
                scaffold_map[scaffold].append((ori, target))

    with open(output_file, "w") as out:
        for i, (scaffold, pairs) in enumerate(scaffold_map.items(), 1):
            if not pairs:
                continue

            ori, target = pairs[0]

            if len(pairs) == 1:
                with open(original_agp) as orig_f:
                    for line in orig_f:
                        if line.startswith("#") or line.strip() == "":
                            out.write(line)
                            continue
                        fields = line.strip().split('\t')
                        if fields[0] == target:
                            fields[0] = scaffold
                            out.write('\t'.join(fields) + "\n")
            else:

                join_string = ",".join([ori + tgt for ori, tgt in pairs]) + "\t" + scaffold

                with open("joins.txt", "w") as jf:
                    jf.write(join_string)

                cmd = ["agptools", "join", "joins.txt", original_agp, "-n100"]
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                for line in result.stdout.splitlines():
                    if scaffold in line:
                        out.write(line + "\n")

                os.remove("joins.txt")
                        

def split_pairs(hap_num, chr_num, map_file, logger):

    tmp_file = f"tmp_{hap_num}.txt"
    cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
    if not run_command([cut_command], shell=True, logger=logger):
        return False

    awk_cmd = f"""awk 'NR==FNR{{lines[$1];next}}
    {{
        if($2!=before_utg){{
            before_utg=$2;
            before_utg_in_set=($2 in lines);
        }}
        if(before_utg_in_set){{
            if(before_utg_in_set && ($4 in lines)){{
                print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6;
            }}
        }}
    }}' {tmp_file} {map_file} > chr{chr_num}g{hap_num}.pairs"""
    if not run_command([awk_cmd], shell=True, logger=logger):
        return False

    cleanup_command = f"rm {tmp_file}"
    if not run_command([cleanup_command], shell=True, logger=logger):
        return False


def run_command(cmd: List[str], logger, cwd: Optional[str] = None, shell: bool = False, stdout=None) -> bool:
    """Execute a shell command."""
    try:
        # logger.info(f"Running command: {cmd}")
        subprocess.run(cmd, cwd=cwd, shell=shell, check=True)
        # logger.info(f"Running command successfully: {cmd}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Error: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error executing command: {str(e)}")
        return False

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

def sort_file(input_file):

    with open(input_file, 'r') as f:
        lines = f.readlines()

    lines.sort()
    with open(input_file, 'w') as f:
        f.writelines(lines)


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
        logger.info(f"Chromosome {chr_num} hap {hap_num} run Get_subgraph_scaffold...")
        logger.info(f"Chromosome {chr_num} hap {hap_num} run Get_subgraph_scaffold: {args.gfa_file} {args.RE_file} {args.digraph_file} group{hap_num}.txt {args.subgraph_file}")

        # Get_subgraph_scaffold_v2(args.gfa_file, args.RE_file, args.digraph_file,f"group{hap_num}.txt", args.subgraph_file)
        Get_subgraph_scaffold(args.digraph_file, args.RE_file, args.HiC_file, f"group{hap_num}.txt", args.subgraph_file, args.gfa_file)

        logger.info(f"Chromosome {chr_num} hap {hap_num} run Get_subgraph_scaffold done...")

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

        # split pairs
        split_pairs(hap_num, chr_num, args.map_file, logger)

        # Run YAHS
        if not run_command([
            "yahs",
            f"chr{chr_num}g{hap_num}.fa",
            f"chr{chr_num}g{hap_num}.pairs",
            "-a", "subgraph_sort.agp",
            "-o", "yahs_iter_1",
            "--file-type", "pa5",  "-q0"
        ], logger=logger):
            return False

        # trans pairs
        Trans_pairs("yahs_iter_1_scaffolds_final.agp", f"chr{chr_num}g{hap_num}.pairs", "yahs_iter_2.pairs") 

        # Create fasta index
        if not run_command(["samtools", "faidx", f"yahs_iter_1_scaffolds_final.fa"], logger=logger):
            return False

        # Run YAHS_iter_2
        if not run_command([
            "yahs",
            "yahs_iter_1_scaffolds_final.fa",
            "yahs_iter_2.pairs",
            "-o", "yahs_iter_2",
            "--file-type", "pa5", "--no-contig-ec", "--no-scaffold-ec", "-q0"
        ], logger=logger):
            return False

        Get_RE("yahs_iter_2_scaffolds_final.fa", "yahs_iter_2", f"{args.enzyme_site}")
        
        # Update scaffold names
        for file_pattern, sed_cmd in [
            ("yahs_iter_2_scaffolds_final.agp", f"s/^/chr{chr_num}g{hap_num}_/g"),
            ("yahs_iter_2_scaffolds_final.fa", f"s/>/>chr{chr_num}g{hap_num}_/g"),
        ]:
            if not run_command(["sed", "-i", sed_cmd, file_pattern], logger=logger):
                return False

        # trans agp
        trans_agp("yahs_iter_2_scaffolds_final.agp", "yahs_iter_1_scaffolds_final.agp", "yahs_iter_2_scaffolds_final_trans.agp")


        # Run get_data_HapHiC_sort.py
        try:
            flag = Get_data_HapHiC_sort(
                "yahs_iter_2.pairs",
                "yahs_iter_2_scaffolds_final.agp",
                "yahs_iter_2.RE_counts.txt",
                f"{args.output_prefix}.chr{chr_num}g{hap_num}",
                args.min_len
            )
            if not flag:
                logger.error(f"Error: Chr{chr_num}g{hap_num} get data for HapHiC")
                return False
        except Exception as e:
            logger.error(f"Error: Chr{chr_num}g{hap_num} get data for HapHiC -> {str(e)}")
            return False


        logger.info(f"Haplotype {hap_num} of chromosome {chr_num} processing completed.")
        return True
    except Exception as e:
        logger.error(f"Error processing chr{chr_num}g{hap_num}: {str(e)}")
        return False

def haphic_sort(pwd: str, args: argparse.Namespace,logger) -> bool:

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
        clm_files = glob.glob(os.path.join(pwd, "chr*/chr*/", f"{args.output_prefix}.*chr*g*.clm"))
        for clm_file in clm_files:
            try:
                os.symlink(clm_file, os.path.basename(clm_file))
            except FileExistsError:
                pass

        # Link RE files
        os.chdir(os.path.join(haphic_dir, "groups_REs"))
        re_files = glob.glob(os.path.join(pwd, "chr*/chr*/", f"{args.output_prefix}.chr*g*scaffold*txt"))
        for re_file in re_files:
            try:
                re = process_one_Contig(re_file)
                if re:
                    one_contigs_list.append(re)
            except FileExistsError:
                pass

        os.chdir(haphic_dir)


        # Merge HT files
        with open(f"{args.output_prefix}.merge.HT.pkl", 'wb') as outfile:
            for ht_file in glob.glob(os.path.join(pwd, "chr*/chr*/*HT*")):
                with open(ht_file, 'rb') as infile:
                    outfile.write(infile.read())

        # Merge scaffold files
        with open(f"{args.output_prefix}.merge.fa", 'w') as outfile:
            for scaffold_file in glob.glob(os.path.join(pwd, "chr*/chr*/yahs_iter_2_scaffolds_final.fa")):
                with open(scaffold_file) as infile:
                    outfile.write(infile.read())

        # Merge trans agp for get final agp
        with open(f"{args.output_prefix}.merge.agp", 'w') as outfile:
            for scaffold_file in glob.glob(os.path.join(pwd, "chr*/chr*/yahs_iter_2_scaffolds_final_trans.agp")):
                with open(scaffold_file) as infile:
                    outfile.write(infile.read())


        # HapHiC sort & build
        script_path_add = os.path.join(script_path, "../src/HapHiC/haphic")
        script_content_1 = f"""
            cd {haphic_dir}
            {script_path_add} sort {args.output_prefix}.merge.fa {args.output_prefix}.merge.HT.pkl split_clms/ groups_REs/*  --mutprob {args.mutprob} --ngen {args.ngen} --npop {args.npop}  --processes {args.processes}
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


        try:
            # get final agp
            trans_agp("scaffolds.agp", f"{args.output_prefix}.merge.agp", "../gphase_final.agp")

            # rescue and connect utg base graph
            os.chdir(os.path.join(haphic_dir, "../"))
            Rescue_base_graph(args.digraph_file, "gphase_final.agp", args.gfa_file, args.RE_file, args.fa_file)
            logger.info("GPhase rescue completed.")

            # sort file
            # sort_file("gphase_final.agp")
            # sort_file("gphase_final_rescue.agp")

            # # get rescue fasta
            # haphic_utils_dir = os.path.join(script_path, "../src/HapHiC/utils")
            # cmd = [f"{haphic_utils_dir}/agp_to_fasta", "gphase_final_rescue.agp", f"{args.fa_file}"]
            # with open("gphase_final_rescue.fasta", "w") as outfile:
            #     if subprocess.run(cmd, stdout=outfile).returncode != 0:
            #         return False

            # rename scaffolds.sort.fa
            src_file = os.path.join(haphic_dir, "scaffolds.sort.fa")
            dst_file = "gphase_final.fasta"
            shutil.copy(src_file, dst_file)
        
        except:
            logger.info("Final conversion error.")
            return False

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
    logger.info("haphic sort and final rescue...")
    if not haphic_sort(pwd, args,logger):
        logger.error("Error in haphic sort and final rescue")
        sys.exit(1)
    
    logger.info("Program completed successfully.")

if __name__ == "__main__":
    main()