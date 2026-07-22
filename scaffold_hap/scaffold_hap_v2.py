#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import re
import subprocess
import shutil
import sys
import tempfile
from pathlib import Path
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

# Output filenames
INTERMEDIATE_UNITIG_AGP = "gphase_final.agp"
FINAL_UNITIG_SCAFFOLD_FASTA = "gphase_final.unitig.scaffold.fasta"
RESCUE_RESCUE_AGP = "gphase_final_rescue.agp"
FINAL_RESCUE_AGP = "gphase_final.unitig.rescue.agp"
RESCUE_CONTIG_AGP = "gphase_final_contig.agp"
FINAL_CONTIG_AGP = "gphase_final.contig.scaffold.agp"
RESCUE_CONTIG_FASTA = "gphase_final_contig.fasta"
FINAL_CONTIG_FASTA = "gphase_final.contig.fasta"
RESCUE_CTG2UTG = "gphase_final_ctg2utg.txt"
FINAL_CTG2UTG = "gphase_final.ctg2utg.txt"
FINAL_CONTIG_SCAFFOLD_FASTA = "gphase_final.contig.scaffold.fasta"
RENAME_COLLAPSE_PREFIX = "gphase_final.unitig.scaffold"
FINAL_UNITIG_SCAFFOLD_AGP = f"{RENAME_COLLAPSE_PREFIX}.agp"
FINAL_UNITIG_FASTA = "gphase_final.unitig.fasta"


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
    output_group.add_argument("-op", "--output_prefix", metavar='\b', required=True, help="Output file prefix")

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

def check_file_exists_and_not_empty(file_path: Union[str, Path], logger: logging.Logger, action_name: str, min_size: int = 0) -> None:
    """
    Check if a file exists and meets a minimum size. Raise an error if it fails.
    :param file_path: Path to the file to check.
    :param logger: The logger instance to use for reporting.
    :param action_name: Descriptive name of the action that was expected to create the file.
    :param min_size: Minimum required size in bytes (default 0).
    :raises FileNotFoundError: If the file does not exist.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        logger.error(f"File check failed: Required file not found: {file_path}")
        raise FileNotFoundError(f"{action_name} failed: Required file not found: {file_path}")
    
    file_size = file_path.stat().st_size
    if file_size <= min_size:
        logger.error(f"File check failed: File is too small (size={file_size} bytes, min_required={min_size} bytes): {file_path}")
        raise EOFError(f"{action_name} failed: File is too small: {file_path}")
        
    logger.info(f"File check passed for: {file_path}")

def trans_agp(agp_file, original_agp, output_file, logger):
    """Translates an AGP file by merging scaffolds using agptools, includes robustness checks."""
    scaffold_map = defaultdict(list)
    action_name = f"AGP translation: {agp_file} to {output_file}"

    # Parse new AGP to get mapping
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            scaffold = fields[0]
            if len(fields) >= 9 and "scaffold" in fields[5]:
                target = fields[5]
                ori = fields[8]
                scaffold_map[scaffold].append((ori, target))

    def flip_ori(o):
        if o == '+':
            return '-'
        if o == '-':
            return '+'
        return o

    def emit_target_lines(scaffold, target_lines, reverse):
        """Rename target AGP rows to scaffold; if reverse, flip order and W orientations."""
        lines = [list(f) for f in target_lines]
        if reverse:
            lines.reverse()
            for fields in lines:
                if fields[4].upper() == 'W':
                    fields[8] = flip_ori(fields[8])
        pos = 1
        out_rows = []
        for part_num, fields in enumerate(lines, 1):
            fields[0] = scaffold
            fields[3] = str(part_num)
            ctype = fields[4].upper()
            if ctype == 'W':
                length = int(fields[7]) - int(fields[6]) + 1
            elif ctype in ('N', 'U'):
                length = int(fields[5])
            else:
                length = int(fields[2]) - int(fields[1]) + 1
            fields[1] = str(pos)
            fields[2] = str(pos + length - 1)
            pos += length
            out_rows.append('\t'.join(fields) + '\n')
        return out_rows

    # Write translated AGP
    with open(output_file, "w") as out:
        for scaffold, pairs in scaffold_map.items():
            if not pairs:
                continue

            ori, target = pairs[0]
            if len(scaffold_map[scaffold]) == 1:
                # Direct rewrite for single components; honor parent orientation
                target_lines = []
                with open(original_agp, 'r') as orig_f:
                    for line in orig_f:
                        if line.startswith("#") or line.strip() == "":
                            continue
                        fields = line.strip().split('\t')
                        if fields[0] == target:
                            target_lines.append(fields)
                for row in emit_target_lines(scaffold, target_lines, reverse=(ori == '-')):
                    out.write(row)
            else:
                # Use agptools join for multiple components
                join_string = ",".join([ori + tgt for ori, tgt in pairs]) + "\t" + scaffold

                with open("joins.txt", "w") as jf:
                    jf.write(join_string)
                # Run agptools
                cmd = ["agptools", "join", "joins.txt", original_agp, "-n100", "-tscaffold", "-eproximity_ligation"]
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

                for line in result.stdout.splitlines():
                    if scaffold in line:
                        out.write(line + "\n")

                if result.stderr:
                    logger.warning(f"agptools join stderr for {scaffold}: {result.stderr.strip()}")
                os.remove("joins.txt")

    # Check if the final output file exists and is non-empty
    check_file_exists_and_not_empty(output_file, logger, action_name, min_size=10)
    return True

def split_pairs(hap_num, chr_num, map_file, logger):
    """Splits HiC pairs file based on the clustered contig list."""
    tmp_file = f"tmp_{hap_num}.txt"
    output_pairs_file = f"chr{chr_num}g{hap_num}.pairs"
    action_name = f"Split pairs for chr{chr_num}g{hap_num}"
    
    # 1. cut command to get contig list
    cut_command = f"cut -f1 group{hap_num}.txt > {tmp_file}"
    if not run_command([cut_command], shell=True, logger=logger, check_files=[tmp_file], action_prefix=action_name):
        return False

    with open(tmp_file, "r") as f:
        tmp_file_line_count = len(f.readlines())

    if tmp_file_line_count == 1:
        with open(output_pairs_file, "w") as f:
            pass 
        return True

    # 2. awk command to filter pairs
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
    }}' {tmp_file} {map_file} > {output_pairs_file}"""
    
    if not run_command([awk_cmd], shell=True, logger=logger, action_prefix=action_name):
        return False

    # 3. cleanup command
    cleanup_command = f"rm {tmp_file}"
    # Cleanup failure is non-fatal
    if not run_command([cleanup_command], shell=True, logger=logger, action_prefix=action_name):
        logger.warning(f"Failed to clean up file: {tmp_file}")
        
    return True


def run_command(cmd: List[str], logger, cwd: Optional[str] = None, shell: bool = False, stdout=None, check_files: Optional[List[str]] = None, action_prefix: str = "Command execution") -> bool:
    """Execute a shell command and optionally check for output file existence."""
    cmd_str = ' '.join(cmd) if not shell else cmd[0]
    
    try:
        # logger.info(f"Running command: {cmd_str}")
        subprocess.run(cmd, cwd=cwd, shell=shell, check=True, stdout=stdout)
        # logger.info(f"Running command successfully: {cmd_str}")

        if check_files:
            for file in check_files:
                # check_file_exists_and_not_empty raises exception on failure
                check_file_exists_and_not_empty(Path(cwd or '.') / file, logger, f"{action_prefix} -> Command: {cmd_str}")
        
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {cmd_str}")
        logger.error(f"Error: {str(e)}")
        return False
    # Catch exceptions raised by the file check
    except (FileNotFoundError, EOFError) as e:
        logger.error(f"Robustness check failed after command: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error executing command: {str(e)}")
        return False

def process_one_Contig(file_path):
    """Processes a single-contig RE counts file to generate a .tour file for HapHiC."""
    dir_path = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    output_path = os.path.dirname(dir_path)
    new_file_name = file_name.replace(".txt", ".tour")
    final_output_path = os.path.join(output_path, new_file_name)
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

        if len(lines) == 2:
            os.symlink(file_path, file_name) 
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
    """Sorts an AGP file first by component and then by start position."""
    lines_to_sort = []
    header_lines = []

    # Read file and distinguish header lines
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    scaffold = parts[0]
                    # Assuming component start is in the second column (AGP standard)
                    try:
                        start = int(parts[1])
                        lines_to_sort.append((scaffold, start, line))
                    except ValueError:
                        # Handle lines where the second column is not an integer (e.g., gaps)
                        lines_to_sort.append((scaffold, -1, line))
                else:
                    # Keep lines that don't fit the expected format as is
                    lines_to_sort.append(("", 0, line))

    # Sort by scaffold name and then by start position
    def sk(name: str):
        m = re.search(r'chr(\d+)g(\d+)', name)
        nums = (int(m.group(1)), int(m.group(2))) if m else (10**9, 0)
        if re.search(r'\.chr\d+g\d+\.scaffold$', name):
            tier = 0         
        elif name.startswith('chr'):
            tier = 1          # chr
        else:
            tier = 2          # utg
        return (tier, nums, name)

    lines_to_sort.sort(key=lambda x: (sk(x[0]), x[1]))

    # Write back to file
    with open(input_file, 'w') as f:
        for hl in header_lines:
            f.write(hl)
        for _, _, line in lines_to_sort:
            f.write(line)


def normalize_agp_n_gaps(agp_file):
    """Convert AGP gap rows with type N to U and set linkage evidence to proximity_ligation.

    AGP gap line: object, beg, end, part, type(N/U), gap_len, gap_type, linkage,
    linkage_evidence. Change type N -> U and last column -> proximity_ligation.
    """
    converted = 0
    lines = []
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                lines.append(line)
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 5 and parts[4].upper() == 'N':
                parts[4] = 'U'
                if len(parts) >= 9:
                    parts[-1] = 'proximity_ligation'
                else:
                    parts.append('proximity_ligation')
                converted += 1
                lines.append('\t'.join(parts) + '\n')
            else:
                lines.append(line if line.endswith('\n') else line + '\n')

    with open(agp_file, 'w') as f:
        f.writelines(lines)
    return converted


def process_chromosome(pwd: str, chr_num: int, args: argparse.Namespace,logger) -> bool:
    """Process a single chromosome directory: symlink files and run cluster2group."""
    logger.info(f"Processing chromosome {chr_num}...")
    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path, "cluster2group.py")
    action_name = f"Process chr{chr_num}"
    
    try:
        # Create chromosome directory
        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        os.makedirs(chr_dir, exist_ok=True)
        os.chdir(chr_dir)

        # Create symbolic links and check existence
        re_counts_file = f"{args.output_prefix}.RE_counts.txt"
        cluster_file = f"{args.output_prefix}.reassign.cluster.txt"
        
        os.symlink(os.path.join("../", f"{args.RE_file}"), re_counts_file)
        os.symlink(
            os.path.join("../", args.cluster_hap_path, f"chr{chr_num}", f"{args.output_prefix}.reassign.cluster.txt"),
            cluster_file
        )
        check_file_exists_and_not_empty(re_counts_file, logger, action_name)
        check_file_exists_and_not_empty(cluster_file, logger, action_name, min_size=10)


        # Run cluster2group.sh
        group_file = f"group1.txt" # Expect at least one group file
        cmd = ["python", script_path_add, "-c", cluster_file, "-r", re_counts_file]
        if not run_command(cmd, logger, check_files=[group_file], action_prefix="cluster2group"):
            return False

        logger.info(f"Chromosome {chr_num} processing completed.")
        return True
    except (FileNotFoundError, EOFError, Exception) as e:
        logger.error(f"Error processing chromosome {chr_num}: {str(e)}")
        return False
    finally:
        os.chdir(pwd)


def process_haplotype(pwd: str, chr_num: int, hap_num: int, args: argparse.Namespace,logger) -> bool:
    """Process a single haplotype: run subgraph, yahs, and initial HapHiC preparation."""
    logger.info(f"Processing haplotype {hap_num} of chromosome {chr_num}...")
    
    fa_file = f"chr{chr_num}g{hap_num}.fa"
    fa_fai = f"{fa_file}.fai"
    
    pairs_file = f"chr{chr_num}g{hap_num}.pairs"
    yahs1_fa = "yahs_iter_1_scaffolds_final.fa"
    yahs1_agp = "yahs_iter_1_scaffolds_final.agp"

    yahs2_fa = "yahs_iter_2_scaffolds_final.fa"
    yahs2_agp = "yahs_iter_2_scaffolds_final.agp"
    yahs2_pairs = "yahs_iter_2.pairs"
    yahs2_re = "yahs_iter_2.RE_counts.txt"
    yahs2_trans_agp = "yahs_iter_2_scaffolds_final_trans.agp"
    subgraph_agp = "subgraph_sort.agp"
    group_file = f"group{hap_num}.txt"
    action_name = f"Process chr{chr_num}g{hap_num}"

    try:
        # Change to haplotype directory
        os.chdir(pwd)
        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        hap_dir = os.path.join(chr_dir, f"chr{chr_num}g{hap_num}")
        os.makedirs(hap_dir, exist_ok=True)
        os.chdir(hap_dir)

        # Symlink creation and checks
        src_files = [
            f"chr{chr_num}/{group_file}", args.subgraph_file, args.gfa_file, args.digraph_file,
            args.RE_file, args.map_file, args.HiC_file, args.fa_file
        ]
        for src in src_files:
            try:
                os.symlink(os.path.join(pwd, src), os.path.basename(src))
            except FileExistsError:
                pass
        check_file_exists_and_not_empty(group_file, logger, action_name, min_size=1)


        # Run Get_subgraph_scaffold
        logger.info(f"{action_name} run Get_subgraph_scaffold...")
        try:
            Get_subgraph_scaffold(args.digraph_file, args.RE_file, args.HiC_file, group_file, args.subgraph_file, args.gfa_file)
            # script_path = os.path.abspath(sys.path[0])
            # logger.info(f"script_path:{script_path}")
            # command = [
            #     "python", os.path.join(script_path, "get_subgraph_scaffold.py"),
            #     "-c", f"{group_file}",
            #     "-l", f"{args.HiC_file}",
            #     "-subgraph", f"{args.subgraph_file}",
            #     "-graph", f"{args.gfa_file}", 
            #     "-digraph", f"{args.digraph_file}", 
            #     "-r", f"{args.RE_file}"
            # ]

            # if not run_command([command], shell=True, logger=logger, action_prefix=action_name):
            #     logger.error(f"Run get_subgraph_scaffold.py error...")
            #     return False
            check_file_exists_and_not_empty(subgraph_agp, logger, "Get_subgraph_scaffold execution", min_size=10)
        except Exception as e:
            logger.error(f"Error in Get_subgraph_scaffold: {str(e)}")
            return False
        logger.info(f"{action_name} run Get_subgraph_scaffold done.")


        # Split FASTA (seqkit)
        tmp_file = f"tmp_{hap_num}.txt"
        cut_command = f"cut -f1 {group_file} > {tmp_file}"
        if not run_command([cut_command], shell=True, logger=logger, check_files=[tmp_file], action_prefix=action_name):
            return False
        seqkit_cmd = f"seqkit grep -f {tmp_file} {args.fa_file} > {fa_file}"
        if not run_command([seqkit_cmd], shell=True, logger=logger, check_files=[fa_file], action_prefix=action_name):
            return False
        cleanup_command = f"rm {tmp_file}"
        run_command([cleanup_command], shell=True, logger=logger, action_prefix=action_name)

        # Create fasta index
        if not run_command(["samtools", "faidx", fa_file], logger=logger, check_files=[fa_fai], action_prefix=action_name):
            return False

        # Split pairs (awk/cut)
        if not split_pairs(hap_num, chr_num, args.map_file, logger):
            return False
        
        # Run YAHS iter 1
        yahs1_cmd = [
            "yahs", fa_file, pairs_file, "-a", subgraph_agp, "-o", "yahs_iter_1",
            "--file-type", "pa5", "-q0"
        ]
        if not run_command(yahs1_cmd, logger=logger, check_files=[yahs1_fa, yahs1_agp], action_prefix="YAHS iter 1"):
            return False

        # Trans pairs (Trans_pairs)
        try:
            with open(yahs1_agp, "r") as f:
                yahs1_agp_line_count = len(f.readlines())

            if yahs1_agp_line_count == 1:
                with open(yahs2_pairs, "w") as f:
                    pass 
            else:
                Trans_pairs(yahs1_agp, pairs_file, yahs2_pairs)
                check_file_exists_and_not_empty(yahs2_pairs, logger, "Trans_pairs script execution", min_size=10)
        except Exception as e:
            logger.error(f"Error in Trans_pairs: {str(e)}")
            return False

        # Create fasta index for iter 1 scaffolds
        if not run_command(["samtools", "faidx", yahs1_fa], logger=logger, check_files=[f"{yahs1_fa}.fai"], action_prefix=action_name):
            return False

        # Run YAHS_iter_2
        yahs2_cmd = [
            "yahs", yahs1_fa, yahs2_pairs, "-o", "yahs_iter_2",
            "--file-type", "pa5", "--no-contig-ec", "--no-scaffold-ec", "-q0"
        ]
        if not run_command(yahs2_cmd, logger=logger, check_files=[yahs2_fa, yahs2_agp], action_prefix="YAHS iter 2"):
            return False

        # Get_RE
        try:
            Get_RE(yahs2_fa, "yahs_iter_2", f"{args.enzyme_site}")
            check_file_exists_and_not_empty(yahs2_re, logger, "Get_RE script execution", min_size=10)
        except Exception as e:
            logger.error(f"Error in Get_RE: {str(e)}")
            return False
        
        # Update scaffold names (sed commands)
        for file_pattern, sed_cmd in [
            (yahs2_agp, f"s/^/chr{chr_num}g{hap_num}_/g"),
            (yahs2_fa, f"s/>/>chr{chr_num}g{hap_num}_/g"),
        ]:
            if not run_command(["sed", "-i", sed_cmd, file_pattern], logger=logger, action_prefix="Sed scaffold rename"):
                return False

        # trans agp
        if not trans_agp(yahs2_agp, yahs1_agp, yahs2_trans_agp, logger):
            return False

        # Run get_data_HapHiC_sort.py
        prefix = f"{args.output_prefix}.chr{chr_num}g{hap_num}"
        haphic_output_files = [
            f"{prefix}.clm", f"{prefix}.HT.pkl", f"{prefix}.txt"
        ]
        
        try:
            flag = Get_data_HapHiC_sort(yahs2_pairs, yahs2_agp, yahs2_re, prefix, args.min_len)
            if not flag:
                logger.error(f"Error: {action_name} failed to get data for HapHiC.")
                return False
        except Exception as e:
            logger.error(f"Error: {action_name} get data for HapHiC -> {str(e)}")
            return False

        logger.info(f"Haplotype {hap_num} of chromosome {chr_num} processing completed.")
        return True
    except (FileNotFoundError, EOFError, Exception) as e:
        logger.error(f"Error processing {action_name}: {str(e)}")
        return False
    finally:
        os.chdir(pwd)

def haphic_sort(pwd: str, args: argparse.Namespace,logger) -> bool:
    """Perform final merging steps: merge files, run HapHiC sort/build, and final rescue."""
    logger.info(f"Starting final merge steps...")
    script_path = os.path.abspath(sys.path[0])
    one_contigs_list = list()
    action_name = "HapHiC sort and final rescue"

    # Define key output files
    merge_fa = f"{args.output_prefix}.merge.fa"
    merge_ht = f"{args.output_prefix}.merge.HT.pkl"
    merge_agp = f"{args.output_prefix}.merge.agp"
    scaffolds_fa = "scaffolds.fa"
    scaffolds_sort_fa = "scaffolds.sort.fa"
    final_agp = INTERMEDIATE_UNITIG_AGP
    final_fasta = FINAL_UNITIG_SCAFFOLD_FASTA
    
    try:
        os.chdir(pwd)
        haphic_dir = os.path.join(pwd, "HapHiC_sort")
        os.makedirs(os.path.join(haphic_dir, "split_clms"), exist_ok=True)
        os.makedirs(os.path.join(haphic_dir, "groups_REs"), exist_ok=True)
        
        # Link CLM files
        os.chdir(os.path.join(haphic_dir, "split_clms"))
        clm_files = glob.glob(os.path.join(pwd, "chr*/chr*/", f"{args.output_prefix}.*chr*g*.clm"))
        for clm_file in clm_files:
            try:
                os.symlink(clm_file, os.path.basename(clm_file))
            except FileExistsError:
                pass

        # Link RE files and create tour for singletons
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
        with open(merge_ht, 'wb') as outfile:
            for ht_file in glob.glob(os.path.join(pwd, "chr*/chr*/*HT*")):
                with open(ht_file, 'rb') as infile:
                    outfile.write(infile.read())
        check_file_exists_and_not_empty(merge_ht, logger, "Merging HT files", min_size=100)

        # Merge scaffold files
        with open(merge_fa, 'w') as outfile:
            for scaffold_file in glob.glob(os.path.join(pwd, "chr*/chr*/yahs_iter_2_scaffolds_final.fa")):
                with open(scaffold_file) as infile:
                    outfile.write(infile.read())
        check_file_exists_and_not_empty(merge_fa, logger, "Merging FASTA files", min_size=100)

        # Merge trans agp for get final agp
        with open(merge_agp, 'w') as outfile:
            for scaffold_file in glob.glob(os.path.join(pwd, "chr*/chr*/yahs_iter_2_scaffolds_final_trans.agp")):
                with open(scaffold_file) as infile:
                    outfile.write(infile.read())
        check_file_exists_and_not_empty(merge_agp, logger, "Merging AGP files", min_size=100)


        # HapHiC sort
        script_path_add = os.path.join(script_path, "../src/HapHiC/haphic")
        script_content_1 = f"""
            cd {haphic_dir}
            {script_path_add} sort {merge_fa} {merge_ht} split_clms/ groups_REs/* --mutprob {args.mutprob} --ngen {args.ngen} --npop {args.npop}  --processes {args.processes} 
        """
        with tempfile.NamedTemporaryFile(delete=False) as temp_script:
            temp_script.write(script_content_1.encode())
            temp_script_path = temp_script.name

        if not run_command(["bash", temp_script_path], logger=logger, action_prefix="HapHiC sort"):
            return False
        # Check for final_tours directory (minimal check for existence)
        check_file_exists_and_not_empty(os.path.join(haphic_dir, "final_tours"), logger, "HapHiC sort output", min_size=0)
        logger.info("HapHiC sort completed.")

        # Copy singleton tour files
        for group in one_contigs_list:
            if not os.path.exists(group):
                os.symlink("../" + group, group)  

        os.chdir(haphic_dir)
        
        # HapHiC build
        script_content_2 = f"""
            cd {haphic_dir}
            {script_path_add} build {merge_fa} {merge_fa} {merge_fa} *tour
            seqkit sort {scaffolds_fa} > {scaffolds_sort_fa}
            samtools faidx {scaffolds_sort_fa}
        """
        with tempfile.NamedTemporaryFile(delete=False) as temp_script:
            temp_script.write(script_content_2.encode())
            temp_script_path = temp_script.name

        build_check_files = [scaffolds_fa, scaffolds_sort_fa, f"{scaffolds_sort_fa}.fai"]
        if not run_command(["bash", temp_script_path], logger=logger, check_files=build_check_files, action_prefix="HapHiC build"):
            return False
        logger.info("HapHiC build completed.")

        # Get final agp
        if not trans_agp("scaffolds.agp", merge_agp, final_agp, logger):
            logger.error("Final AGP translation failed.")
            return False

        # Post-processing and Rescue steps
        os.chdir(os.path.join(haphic_dir, "../"))

        # Sort files
        sort_file(f"HapHiC_sort/{final_agp}")
        Path(f"HapHiC_sort/{final_agp}").rename(f"{final_agp}")

        # Rename scaffolds.sort.fa
        src_file = os.path.join(haphic_dir, scaffolds_sort_fa)
        shutil.copy(src_file, final_fasta)
        check_file_exists_and_not_empty(final_fasta, logger, "Final FASTA copy", min_size=100)

        return True
    except (FileNotFoundError, EOFError, Exception) as e:
        logger.error(f"Error during final merge: {str(e)}")
        return False
    finally:
        os.chdir(pwd)

def log_start(logger: logging.Logger, script_name: str, version: str, args: argparse.Namespace):
    """Log the start of the program."""
    logger.info(f"Program started, {script_name} version: {version}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def fix_reversed_singleton_utg_scaffolds(agp_path, scaffold_fasta, utg_fasta, logger):
    """
    For scaffolds that contain exactly one W component with orientation '+':
    if the scaffold sequence equals the reverse complement of the source utg
    (a YaHS/trans_agp orientation sync bug), reverse-complement the FASTA entry
    back so it matches the original utg direction. AGP is left unchanged.
    """
    scaffold_components = defaultdict(list)
    with open(agp_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            scaffold_components[parts[0]].append(parts)

    # scaffold_id -> (utg_id, contig_start, contig_end)
    singleton_plus = {}
    for scaffold_id, comps in scaffold_components.items():
        w_comps = [c for c in comps if c[4].upper() == "W"]
        if len(w_comps) != 1:
            continue
        w = w_comps[0]
        if w[8] != "+":
            continue
        try:
            start, end = int(w[6]), int(w[7])
        except ValueError:
            continue
        singleton_plus[scaffold_id] = (w[5], start, end)

    if not singleton_plus:
        logger.info("No singleton '+' scaffolds to check for reversed orientation.")
        return 0

    needed_utgs = {utg_id for utg_id, _, _ in singleton_plus.values()}
    utg_seq = {}
    current_id = None
    current_seq = []
    with open(utg_fasta, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_id in needed_utgs:
                    utg_seq[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                if current_id in needed_utgs:
                    current_seq.append(line)
        if current_id and current_id in needed_utgs:
            utg_seq[current_id] = "".join(current_seq)

    records = []
    current_id = None
    current_seq = []
    with open(scaffold_fasta, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            records.append((current_id, "".join(current_seq)))

    fixed = 0
    new_records = []
    for scaffold_id, seq in records:
        info = singleton_plus.get(scaffold_id)
        if not info:
            new_records.append((scaffold_id, seq))
            continue
        utg_id, start, end = info
        src = utg_seq.get(utg_id)
        if src is None:
            logger.warning(
                f"Singleton scaffold {scaffold_id}: utg {utg_id} not found in {utg_fasta}"
            )
            new_records.append((scaffold_id, seq))
            continue
        expected = src[start - 1 : end]
        if len(seq) == len(expected) and seq == _reverse_complement(expected) and seq != expected:
            new_records.append((scaffold_id, _reverse_complement(seq)))
            fixed += 1
            logger.info(
                f"Reversed singleton scaffold {scaffold_id} back to match {utg_id} (+)"
            )
        else:
            new_records.append((scaffold_id, seq))

    if fixed:
        with open(scaffold_fasta, "w") as out:
            for scaffold_id, seq in new_records:
                out.write(f">{scaffold_id}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + "\n")

    logger.info(f"Fixed {fixed} reversed singleton-utg scaffold(s) in {scaffold_fasta}")
    return fixed


def append_unplaced_utgs(agp_path, all_utg_fasta, target_fasta):
    """
    Identify the Unitigs present in the input FASTA but absent from the final AGP,
    and append them to the end of both the final AGP file and the final FASTA.
    """
    used_utgs = set()
    if os.path.exists(agp_path):
        with open(agp_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) >= 9 and parts[4].upper() not in ['N', 'U']:
                    used_utgs.add(parts[5].strip())

    def append_to_files(utg_id, sequence, f_agp, f_fa):
        length = len(sequence)
        f_agp.write(f"{utg_id}\t1\t{length}\t1\tW\t{utg_id}\t1\t{length}\t+\n")
        f_fa.write(f">{utg_id}\n")
        for i in range(0, length, 80):
            f_fa.write(sequence[i:i+80] + "\n")

    count = 0
    with open(agp_path, 'a') as f_agp, open(target_fasta, 'a') as f_fa:
        current_id = None
        current_seq = []
        
        with open(all_utg_fasta, 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_id not in used_utgs:
                        append_to_files(current_id, "".join(current_seq), f_agp, f_fa)
                        count += 1
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id and current_id not in used_utgs:
                append_to_files(current_id, "".join(current_seq), f_agp, f_fa)
                count += 1
                
    return count

def main():
    args = parse_arguments()
    logger = setup_logging('scaffold_hap.log') 
    pwd = os.getcwd()
    
    # Log program start
    log_start(logger, "scaffold_hap.py", "1.0.0", args)
    
    # Process chromosomes
    for i in range(1, args.chr_number + 1):
        logger.info(f"Starting processing for chromosome {i}...")
        if not process_chromosome(pwd, i, args ,logger):
            logger.error(f"Error processing chromosome {i}. Exiting.")
            sys.exit(1)
    
    # Process haplotypes in parallel
    logger.info(f"Starting parallel haplotype processing with {args.thread_number} threads...")
    with ProcessPoolExecutor(max_workers=args.thread_number) as executor:
        futures = []
        for i in range(1, args.chr_number + 1):
            for j in range(1, args.hap_number + 1):
                futures.append(executor.submit(process_haplotype, pwd, i, j, args,logger))
        
        for future in as_completed(futures):
            try:
                if not future.result():
                    logger.error("Error processing a haplotype. Exiting.")
                    sys.exit(1)
            except Exception as e:
                logger.error(f"An exception occurred during parallel haplotype processing: {e}")
                sys.exit(1)


    # haphic sort & build
    logger.info("Starting HapHiC sort...")
    if not haphic_sort(pwd, args,logger):
        logger.error("Error in HapHiC sort and final rescue. Exiting.")
        sys.exit(1)

    # add unitigs discarded during the clustering stage
    added_num = append_unplaced_utgs(INTERMEDIATE_UNITIG_AGP, args.fa_file, FINAL_UNITIG_SCAFFOLD_FASTA)
    logger.info(f"A total of {added_num} missing unitigs have been added to the AGP and FASTA files.")

    # Singleton utg scaffolds: if FASTA is RC of source utg while AGP says '+', flip FASTA back
    fixed_num = fix_reversed_singleton_utg_scaffolds(
        INTERMEDIATE_UNITIG_AGP,
        FINAL_UNITIG_SCAFFOLD_FASTA,
        args.fa_file,
        logger,
    )
    logger.info(f"Singleton orientation fix done ({fixed_num} scaffold(s) corrected).")

    # GPhase rescue using assembly graph 
    logger.info("Starting final rescue...")
    try:
        # Rescue and connect utg base graph
        try:
            Rescue_base_graph(args.digraph_file, INTERMEDIATE_UNITIG_AGP, args.gfa_file, args.RE_file, args.fa_file)
            check_file_exists_and_not_empty(RESCUE_RESCUE_AGP, logger, "Rescue_base_graph execution", min_size=100)
            check_file_exists_and_not_empty(RESCUE_CONTIG_AGP, logger, "Rescue_base_graph execution", min_size=100)
            check_file_exists_and_not_empty(RESCUE_CONTIG_FASTA, logger, "Rescue_base_graph execution", min_size=100)
            check_file_exists_and_not_empty(RESCUE_CTG2UTG, logger, "Rescue_base_graph execution", min_size=10)

            Path(RESCUE_RESCUE_AGP).rename(FINAL_RESCUE_AGP)
            Path(RESCUE_CONTIG_AGP).rename(FINAL_CONTIG_AGP)
            Path(RESCUE_CONTIG_FASTA).rename(FINAL_CONTIG_FASTA)
            Path(RESCUE_CTG2UTG).rename(FINAL_CTG2UTG)
            check_file_exists_and_not_empty(FINAL_RESCUE_AGP, logger, "Rename rescue AGP", min_size=100)
            check_file_exists_and_not_empty(FINAL_CONTIG_AGP, logger, "Rename contig AGP", min_size=100)
            check_file_exists_and_not_empty(FINAL_CONTIG_FASTA, logger, "Rename contig fasta", min_size=100)
            check_file_exists_and_not_empty(FINAL_CTG2UTG, logger, "Rename ctg2utg", min_size=10)
            logger.info("GPhase rescue completed.")
        except Exception as e:
            logger.error(f"Error in Rescue_base_graph: {str(e)}")
            sys.exit(1)

        script_path = os.path.abspath(sys.path[0])
        haphic_utils_dir = os.path.join(script_path, "../src/HapHiC/utils")

        # Rename duplicated collapsed unitigs/contigs for Hi-C heatmap
        logger.info("Starting rename_collapse for duplicated unitigs...")
        rename_script = os.path.join(script_path, "rename_collapse_agp_pairs_fasta.py")
        rename_unitig_fa = f"{RENAME_COLLAPSE_PREFIX}.fa"
        rename_rescue_agp = f"{RENAME_COLLAPSE_PREFIX}.rescue.agp"
        rename_contig_agp = f"{RENAME_COLLAPSE_PREFIX}.contig.agp"
        rename_contig_fa = f"{RENAME_COLLAPSE_PREFIX}.contig.fa"
        rename_cmd = [
            "python", rename_script,
            INTERMEDIATE_UNITIG_AGP,
            args.fa_file,
            RENAME_COLLAPSE_PREFIX,
            "--no-hic",
            "--rescue-agp", FINAL_RESCUE_AGP,
            "--contig-agp", FINAL_CONTIG_AGP,
            "--contig-fasta", FINAL_CONTIG_FASTA,
        ]
        rename_outputs = [
            FINAL_UNITIG_SCAFFOLD_AGP,
            rename_unitig_fa,
            rename_rescue_agp,
            rename_contig_agp,
            rename_contig_fa,
        ]
        if not run_command(
            rename_cmd,
            logger=logger,
            check_files=rename_outputs,
            action_prefix="rename_collapse_agp_pairs_fasta",
        ):
            logger.error("rename_collapse_agp_pairs_fasta failed.")
            sys.exit(1)

        Path(rename_unitig_fa).replace(FINAL_UNITIG_FASTA)
        Path(rename_rescue_agp).replace(FINAL_RESCUE_AGP)
        Path(rename_contig_agp).replace(FINAL_CONTIG_AGP)
        Path(rename_contig_fa).replace(FINAL_CONTIG_FASTA)
        check_file_exists_and_not_empty(FINAL_UNITIG_FASTA, logger, "Rename unitig fasta", min_size=100)
        check_file_exists_and_not_empty(FINAL_RESCUE_AGP, logger, "Rename rescue AGP", min_size=100)
        check_file_exists_and_not_empty(FINAL_CONTIG_AGP, logger, "Rename contig AGP", min_size=100)
        check_file_exists_and_not_empty(FINAL_CONTIG_FASTA, logger, "Rename contig fasta", min_size=100)
        logger.info("rename_collapse completed.")

        sort_file(FINAL_RESCUE_AGP)
        sort_file(FINAL_CONTIG_AGP)
        sort_file(FINAL_UNITIG_SCAFFOLD_AGP)

        # Normalize gap rows: N -> U, linkage evidence -> proximity_ligation
        logger.info("Normalizing AGP N gaps to U / proximity_ligation...")
        for agp_path in (FINAL_UNITIG_SCAFFOLD_AGP, FINAL_RESCUE_AGP, FINAL_CONTIG_AGP):
            n_converted = normalize_agp_n_gaps(agp_path)
            logger.info(f"Converted {n_converted} N-gap row(s) in {agp_path}")

        # Build contig-level scaffold FASTA from renamed contig AGP/FASTA
        for agp, fa, out in (
            (FINAL_CONTIG_AGP, FINAL_CONTIG_FASTA, FINAL_CONTIG_SCAFFOLD_FASTA),
            (FINAL_UNITIG_SCAFFOLD_AGP, FINAL_UNITIG_FASTA, FINAL_UNITIG_SCAFFOLD_FASTA),
        ):
            cmd = [f"{haphic_utils_dir}/agp_to_fasta", agp, fa]
            with open(out, "w") as outfile:
                if subprocess.run(cmd, stdout=outfile).returncode != 0:
                    logger.error(f"agp_to_fasta failed for {out}")
                    sys.exit(1)
            check_file_exists_and_not_empty(out, logger, "agp_to_fasta execution", min_size=100)

        Path(INTERMEDIATE_UNITIG_AGP).unlink(missing_ok=True)
        logger.info(f"Removed intermediate file: {INTERMEDIATE_UNITIG_AGP}")

    except Exception as e:
        logger.error(f"An exception occurred during rescue using assembly graph: {e}")
        sys.exit(1)

    logger.info("GPhase completed successfully.")

if __name__ == "__main__":
    main()