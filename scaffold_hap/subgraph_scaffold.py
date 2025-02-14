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
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Genome scaffolding pipeline")
    
    parser.add_argument("-r", "--RE_file", required=True,help="Path to RE file")
    parser.add_argument("-s", "--group_ctgs_file",required=True, help="Group contigs file")
    parser.add_argument("-CHP", "--cluster_hap_path",required=True, help="Path to cluster hap directory")
    parser.add_argument("-g", "--gfa_file", required=True,help="GFA file")
    parser.add_argument("-d", "--digraph_file", required=True,help="Digraph file")
    # parser.add_argument("-p", "--paris_file",required=True, help="Paris file")
    parser.add_argument("--map_file_type",default="bam", choices=['bam', 'pa5'], type=str.lower,help="hic mapping file type")
    parser.add_argument("-m", "--map_file", required=True,help="hic mapping file")
    parser.add_argument("-l", "--HiC_file", required=True,help="HiC file")
    parser.add_argument("-f", "--fa_file", required=True,help="Fasta file")
    parser.add_argument("-nc", "--n_chr", required=True,type=int, help="Number of chromosomes")
    parser.add_argument("-nh", "--n_hap", required=True,type=int, help="Number of haplotypes")
    parser.add_argument("-op", "--op", required=True, help="Operation prefix")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of parallel processes")
    
    return parser.parse_args()

def run_command(cmd: List[str], cwd: Optional[str] = None, shell: bool = False) -> bool:
    """Execute a shell command."""
    try:
        subprocess.run(cmd, cwd=cwd, shell=shell, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}")
        print(f"Error: {e}")
        return False
    except Exception as e:
        print(f"Error executing command: {e}")
        return False

def process_chromosome(pwd: str, chr_num: int, args: argparse.Namespace) -> bool:
    """Process a single chromosome directory."""
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "cluster2group.py")
    try:
        # Create chromosome directory
        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        os.makedirs(chr_dir, exist_ok=True)
        os.chdir(chr_dir)

        # Create symbolic links
        os.symlink(os.path.join("../", f"{args.RE_file}"), f"{args.op}.RE_counts.txt")
        os.symlink(
            os.path.join("../", args.cluster_hap_path, f"chr{chr_num}", "reassign_collapse.cluster.txt"),
            "reassign_collapse.cluster.txt"
        )

        # Run cluster2group.sh
        cmd = ["python", script_path_add, "-c", "reassign_collapse.cluster.txt", "-r", f"{args.op}.RE_counts.txt"]
        if not run_command(cmd):
            return False

        return True
    except Exception as e:
        print(f"Error processing chromosome {chr_num}: {e}")
        return False

def process_haplotype(pwd: str, chr_num: int, hap_num: int, args: argparse.Namespace) -> bool:
    """Process a single haplotype."""
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
            args.group_ctgs_file,
            args.gfa_file,
            args.digraph_file,
            args.RE_file,
            args.map_file,
            args.map_file+".bai",
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
        script_path_add = os.path.join(script_path, "scaffold.subgraph.py")

        cmd = [
            "python", script_path_add,
            "-c", f"group{hap_num}.txt",
            "-subgraph", args.group_ctgs_file,
            "-graph", args.gfa_file,
            "-digraph", args.digraph_file,
            "-r", args.RE_file,
            "-l", args.HiC_file
        ]

        if not run_command(cmd):
            return False
        # split asm.fa
        tmp_file = f"tmp_{hap_num}.txt"
        cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
        if not run_command([cut_command], shell=True):
            return False
        # seqkit 
        seqkit_cmd = f"seqkit grep -f {tmp_file} {args.fa_file} > chr{chr_num}g{hap_num}.fa"
        if not run_command([seqkit_cmd], shell=True):
            return False
        cleanup_command = f"rm {tmp_file}"
        if not run_command([cleanup_command], shell=True):
            return False

        # Create fasta index
        if not run_command(["samtools", "faidx", f"chr{chr_num}g{hap_num}.fa"]):
            return False

        # Process bam or pairs file
        if args.map_file_type == "bam":
            command = f"cut -f1 group{hap_num}.txt | xargs -I T echo T | tr '\n' ' ' "
            cut_output = subprocess.check_output(command, shell=True).decode('utf-8').strip()

            # Use it in the samtools command
            samtools_cmd = f"samtools view -b {args.map_file} {cut_output} -o chr{chr_num}g{hap_num}.bam"
            if not run_command([samtools_cmd], shell=True):
                return False

            # Run YAHS
            if not run_command([
                "yahs",
                f"chr{chr_num}g{hap_num}.fa",
                f"chr{chr_num}g{hap_num}.bam",
                "-a", "subgraphGroup.agp",
                "-o", "subgraph_yahs"
            ]):
                return False

        if args.map_file_type == "pa5":

            tmp_file = f"tmp_{hap_num}.txt"
            cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
            if not run_command([cut_command], shell=True):
                return False

            # 使用 awk 命令
            awk_cmd = f"""awk 'NR==FNR{{lines[$1];next}}{{if(NF < 8 && $2 in lines){{print $0;next}}if($2 in lines && $4 in lines){{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5}}}}' {tmp_file} {args.map_file} > chr{chr_num}g{hap_num}.pairs"""
            if not run_command([awk_cmd], shell=True):
                return False

            cleanup_command = f"rm {tmp_file}"
            if not run_command([cleanup_command], shell=True):
                return False

            # Run YAHS
            if not run_command([
                "yahs",
                f"chr{chr_num}g{hap_num}.fa",
                f"chr{chr_num}g{hap_num}.pairs",
                "-a", "subgraphGroup.agp",
                "-o", "subgraph_yahs",
                "--file-type", "pa5"
            ]):
                return False

        # Update scaffold names
        for file_pattern, sed_cmd in [
            ("subgraph_yahs_scaffolds_final.agp", f"s/^/chr{chr_num}g{hap_num}_/g"),
            ("subgraph_yahs_scaffolds_final.fa", f"s/>/>chr{chr_num}g{hap_num}_/g")
        ]:
            if not run_command(["sed", "-i", sed_cmd, file_pattern]):
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
                "-o", f"{args.op}.chr{chr_num}g{hap_num}"
            ]):
                return False
        if args.map_file_type == "pa5":
            if not run_command([
                "python", f"{script_path_add}",
                "-m", f"chr{chr_num}g{hap_num}.pairs",
                "-a", "subgraph_yahs_scaffolds_final.agp",
                "-r", args.RE_file,
                "-o", f"{args.op}.chr{chr_num}g{hap_num}"
            ]):
                return False

        return True
    except Exception as e:
        print(f"Error processing chr{chr_num}g{hap_num}: {e}")
        return False

def perform_final_merge(pwd: str, args: argparse.Namespace) -> bool:
    """Perform final merging steps."""
    try:
        os.chdir(pwd)
        haphic_dir = os.path.join(pwd, "HapHiC_sort")
        os.makedirs(os.path.join(haphic_dir, "split_clms"), exist_ok=True)
        os.makedirs(os.path.join(haphic_dir, "groups_REs"), exist_ok=True)
        os.chdir(os.path.join(haphic_dir, "split_clms"))

        # Link CLM files
        clm_files = glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort", f"{args.op}.*chr*g*.clm"))
        for clm_file in clm_files:
            try:
                os.symlink(clm_file, os.path.basename(clm_file))
            except FileExistsError:
                pass

        # Link RE files
        os.chdir(os.path.join(haphic_dir, "groups_REs"))
        re_files = glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort", f"{args.op}.chr*g*scaffold*txt"))
        for re_file in re_files:
            try:
                os.symlink(re_file, os.path.basename(re_file))
            except FileExistsError:
                pass

        os.chdir(haphic_dir)

        # Merge HT files
        with open(f"{args.op}.merge.HT.pkl", 'wb') as outfile:
            for ht_file in glob.glob(os.path.join(pwd, "chr*/chr*/scaffold_HapHiC_sort/*HT*")):
                with open(ht_file, 'rb') as infile:
                    outfile.write(infile.read())

        # Merge scaffold files
        with open(f"{args.op}.merge.fa", 'w') as outfile:
            for scaffold_file in glob.glob(os.path.join(pwd, "chr*/chr*/subgraph_yahs_scaffolds_final.fa")):
                with open(scaffold_file) as infile:
                    outfile.write(infile.read())

        # Merge pairs files
        with open(f"{args.op}.merge.map.bin", 'wb') as outfile:
            for pairs_file in glob.glob(os.path.join(pwd, "chr*/chr*/subgraph_yahs.bin")):
                with open(pairs_file, 'rb') as infile:
                    outfile.write(infile.read())

        # Create and execute final commands script
        script_path_add = os.path.join(script_path, "../src/HapHiC/haphic")
        script_content = f"""
            cd {haphic_dir}
            {script_path_add} sort {args.op}.merge.fa {args.op}.merge.HT.pkl split_clms/ groups_REs/* --quick_view
            {script_path_add} build {args.op}.merge.fa {args.op}.merge.fa {args.op}.merge.pairs.bin final_tours/*tour
            seqkit sort scaffolds.fa > scaffolds.sort.fa
            samtools faidx scaffolds.sort.fa
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as script_file:
            script_file.write(script_content)
            script_path = script_file.name

        if not run_command(["bash", script_path]):
            return False

        os.unlink(script_path)
        return True
    except Exception as e:
        print(f"Error in final merge steps: {e}")
        return False

def main():
    args = parse_arguments()
    pwd = os.getcwd()

    # Process chromosomes
    for i in range(1, args.n_chr + 1):
        if not process_chromosome(pwd, i, args):
            sys.exit(1)

    # Process haplotypes in parallel
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for i in range(1, args.n_chr + 1):
            for j in range(1, args.n_hap + 1):
                futures.append(executor.submit(process_haplotype, pwd, i, j, args))

        for future in as_completed(futures):
            if not future.result():
                sys.exit(1)

    # Perform final merge
    if not perform_final_merge(pwd, args):
        sys.exit(1)

if __name__ == "__main__":
    main()