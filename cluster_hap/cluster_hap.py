#!/usr/bin/env python3

import os
import argparse
import subprocess
import logging
import pandas as pd
import argcomplete
import networkx as nx
from argcomplete.completers import FilesCompleter
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor
import functools
import sys
from collections import defaultdict
from find_knees import find_best_knee
from multilevel_cluster_v2 import multilevel_cluster
from louvain_reassign_allele import louvain_reassign_allele
from louvain_nei import louvain_nei


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
    performance_group.add_argument("--process_num", metavar='\b',type=int, default=20,
                       help=f"Number of parallel processes (default: 20)")
    performance_group.add_argument("--thread_num", metavar='\b',type=int,default=20,
                       help=f"Number of threads per process (default: 20)")
    
    reassign_group = parser.add_argument_group('>>> Parameters for reassign step')
    reassign_group.add_argument("--reassign_number", metavar='\b',type=int, default=1,
                       help=f"Number of reassign step (default: 1)")

    Optional_group  = parser.add_argument_group('>>> Optional parameters')
    Optional_group.add_argument('--correct', action='store_true', help='correct the cluster')
    Optional_group.add_argument('--expand', action='store_true', help='expand the allele')
    Optional_group.add_argument('--rescue', action='store_true', help='rescue filtered subgraph ')


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
        logger.info(f"Created symlink: {source} -> {dest}")
    else:
        logger.warning(f"Symlink already exists: {dest}")

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

def adjust_r_and_cluster(initial_r, min_r, max_r, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger):
    r = initial_r
    check_cluster_dict = defaultdict(list)
    while min_r <= r <= max_r:
        logger.info(f"Running clustering with r={r}")
        check_cluster_dict, max_group_allele_value = multilevel_cluster(csv_file, cluster_output, r , "check", utg_file, partig_file, hap_number)

        num_clusters = len(check_cluster_dict)
        logger.info(f"Number of clusters: {num_clusters} (Target: {hap_number})")

        if num_clusters == hap_number:
            logger.info(f"Optimal r found: {r}")
            return r
        elif num_clusters > hap_number:
            max_r = r - step
        else:
            min_r = r + step

        r = (min_r + max_r) / 2

    return 0

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len


def read_collapse_num(collapse_num_file):

    collapse_num_dict = defaultdict()
    with open(collapse_num_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0].startswith("utg") or line[0].startswith("utig") :
                try:
                    collapse_num_dict[line[0]] = int(line[1])
                except:
                    collapse_num_dict[line[0]] = 1
    return collapse_num_dict

def get_avg_uncollapse_num(REFile, collapse_num_file, hap_number):

    ctg_RE_len = read_REs(REFile)
    collapse_num_dict = read_collapse_num(collapse_num_file)

    uncollapse_avg = sum([ list_[1] for ctg, list_ in ctg_RE_len.items() if collapse_num_dict[ctg] < 2]) / int(hap_number)

    return uncollapse_avg



def multiple_adjust_r_and_cluster(initial_r, min_r, max_r, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger):

    section_length = float((max_r - min_r) / 8)

    section_1 = adjust_r_and_cluster((min_r + min_r+section_length)/2, min_r, min_r+section_length, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_2 = adjust_r_and_cluster((min_r+section_length + min_r+section_length*2)/2, min_r+section_length, min_r+section_length*2, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_3 = adjust_r_and_cluster((min_r+section_length*2 + min_r+section_length*3)/2, min_r+section_length*2, min_r+section_length*3, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_4 = adjust_r_and_cluster((min_r+section_length*3 + min_r+section_length*4)/2, min_r+section_length*3, min_r+section_length*4, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)

    section_5 = adjust_r_and_cluster((min_r+section_length*4 + min_r+section_length*5)/2, min_r+section_length*4, min_r+section_length*5, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_6 = adjust_r_and_cluster((min_r+section_length*5 + min_r+section_length*6)/2, min_r+section_length*5, min_r+section_length*6, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_7 = adjust_r_and_cluster((min_r+section_length*6 + min_r+section_length*7)/2, min_r+section_length*6, min_r+section_length*7, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)
    section_8 = adjust_r_and_cluster((min_r+section_length*7 + max_r)/2, min_r+section_length*7, max_r, step, cluster_output, csv_file, utg_file, partig_file, hap_number,logger)


    if section_1 or section_2 or section_3 or section_4 or section_5 or section_6 or section_7 or section_8:
        list_ = [section_1, section_2, section_3, section_4, section_5, section_6, section_7, section_8]
        r = max([ r for r in list_ if r > 0])
        return r
    raise 


def filter_edges_by_density(chr_num, HiC_file, utg_rescue_file, filter_HiC_file, logger, step=0.5):

    nodes = pd.read_csv(utg_rescue_file, header=None)[0].tolist()
    edges = pd.read_csv(HiC_file, sep=',', header=0)
    edges.columns = ['source', 'target', 'links']

    num_nodes = len(nodes)
    num_edges = len(edges)
    density = (2 * num_edges) / (num_nodes * (num_nodes - 1)) if num_nodes > 1 else 0
    logger.info(f"Chr{chr_num} hic graph info : edges -> {num_edges}\t nodes -> {num_nodes}\t density -> {density}")

    if  num_nodes < (num_edges / 50):
        logger.info(f"Chr{chr_num} HiC signal filtering...")
        threshold = 0.5
        
        while True:
            filtered_edges = edges[edges['links'] > threshold]
            filtered_num_edges = len(filtered_edges)

            logger.info(f"Chr{chr_num} HiC signal filtering : threshold -> {threshold:.1f}\t edges: -> {filtered_num_edges}")

            if filtered_num_edges < (num_nodes * 50):
                break

            threshold += step

        filtered_edges.to_csv(filter_HiC_file, sep=',', header=False, index=False)
        return threshold
    else:
        edges.to_csv(filter_HiC_file, sep=',', header=False, index=False)
        return 0



def process_chromosome(chr_num, args, pwd, partig_file,logger):
    script_path = os.path.abspath(sys.path[0])
    try:

        chr_dir = os.path.join(pwd, f"chr{chr_num}")
        os.makedirs(chr_dir, exist_ok=True)

        os.chdir(chr_dir)
        logger.info(f"Processing chromosome {chr_num} in {chr_dir}")

        # Create symlinks with ThreadPoolExecutor
        param_files = {
            "chr_cluster_file": f"../{args.chr_cluster_file}",
            "digraph_file": f"../{args.digraph_file}",
            "subgraph_file": f"../{args.subgraph_file}",
            "collapse_num_file": f"../{args.collapse_num_file}",
            "HiC_file": f"../{args.HiC_file}",
            "RE_file": f"../{args.RE_file}",
            "partig_file": f"../{partig_file}"
        }
        # "merge_partig_file": f"../merge.partig.csv"
        no_expand_partig_file = partig_file
        if args.expand:
            origin_partig_file = "merge.partig.csv"
            create_symlink(f"../merge.partig.csv", origin_partig_file, logger)
        else:
            origin_partig_file = partig_file
        
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            futures = []
            for name, filepath in param_files.items():
                dest = os.path.join(chr_dir, os.path.basename(filepath))
                futures.append(executor.submit(create_symlink, filepath, dest, logger))
            
            for future in futures:
                future.result()


        utg_file = f"{args.output_prefix}.chr{chr_num}.utgs.txt"
        create_symlink(f"../group{chr_num}.txt", utg_file, logger)

        if args.rescue:
            utg_rescue_file = f"{args.output_prefix}.chr{chr_num}.utgs.rescue.txt"
            create_symlink(f"../group{chr_num}_rescue.txt", utg_rescue_file, logger)
        else:
            utg_rescue_file = utg_file


        # Process files with ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.thread_num) as executor:
            # Process links file
            links_file = f"{args.output_prefix}.chr{chr_num}.links.all.nor.csv"
            links_future = executor.submit(filter_links_by_utgs,str(chr_num)+"_links", utg_rescue_file, args.HiC_file, links_file, logger)
            
            # Process partig file
            partig_file = f"{args.output_prefix}.chr{chr_num}.partig.csv"
            partig_future = executor.submit(filter_links_by_utgs,str(chr_num)+"_partig", utg_rescue_file, origin_partig_file, partig_file, logger)
            
            # Wait for both operations to complete
            links_future.result()
            partig_future.result()

        # filter HiC links
        filter_edge_file = f"{args.output_prefix}.chr{chr_num}.links.nor.csv"
        cut_value = filter_edges_by_density(chr_num, links_file, utg_rescue_file, filter_edge_file, logger=logger, step=0.5)

        filtered_links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.filterAllele.csv"
        command = (
            f"awk -F '[, \\t]' 'NR==FNR{{lines[$1\",\"$2];next}}!($1\",\"$2 in lines)' "
            f"{partig_file} {filter_edge_file} > {filtered_links_file}"
        )
        execute_command(command, f"Failed to filter allele links for {filtered_links_file}",logger)
        execute_command(f"sed '1isource,target,links' -i {filtered_links_file}", f"Failed to add header to {filtered_links_file}",logger)


        # # The knee is used to filter HiC signals
        best_knee = find_best_knee(filtered_links_file, f"{args.output_prefix}.chr{chr_num}.knees")
        cut_value_step = 1

        if best_knee < 5:
            cut_value_max = 5
        else:
            cut_value_max = best_knee


        cluster_output = f"{args.output_prefix}.chr{chr_num}.cluster.txt"
        # 计算 uncollapse 的长度和
        avg_uncollapse_num = get_avg_uncollapse_num(utg_rescue_file, args.collapse_num_file, args.hap_number)
        logger.info(f"chr_num:{chr_num}\tbest_knee:{best_knee}")

        # 记录每个满足要求的 cluster，在遍历至best_knee时则找最优的clusters
        # tuple ： (check_cluster_dict, max_group_allele_value）
        satisfy_clusters_list = list()
        check_cluster_dict = defaultdict(list)

        while cut_value < cut_value_max:

            cut_links_file = f"{args.output_prefix}.chr{chr_num}.links.nor.filterAllele.c{float(cut_value)}.csv"
            command = (
                f"awk -F ',' '($3> {float(cut_value)})' "
                f"{filtered_links_file} > {cut_links_file}"
            )
            execute_command(command, f"Failed to filter hic links for {cut_links_file}",logger)

            try:
                # # Run louvain_nei.py
                louvain_nei_result = louvain_nei(args.collapse_num_file, utg_file, cut_links_file, partig_file)
                if not louvain_nei_result:
                    raise
                check_cluster_dict, max_group_allele_value = multilevel_cluster("louvain_nei.csv", cluster_output, float(1), "check", utg_file, no_expand_partig_file, int(args.hap_number))

                optimal_r = multiple_adjust_r_and_cluster(
                    initial_r=1.0,
                    min_r=0.01,
                    max_r=3,
                    step=0.03,
                    cluster_output=cluster_output,
                    csv_file="louvain_nei.csv", 
                    utg_file=utg_file,
                    partig_file=partig_file,
                    hap_number=args.hap_number, logger=logger
                )
                check_cluster_dict, max_group_allele_value = multilevel_cluster("louvain_nei.csv", cluster_output, float(optimal_r), "check", utg_file, no_expand_partig_file, int(args.hap_number))
                
                # Check the number of clusters
                num_clusters = len(check_cluster_dict)
                if num_clusters != args.hap_number:
                    raise
                else:
                    satisfy_clusters_list.append((check_cluster_dict, max_group_allele_value))

                # if max_group_allele_value > min([int(avg_uncollapse_num/20), 1000000]):
                if max_group_allele_value > 1000000:
                    logger.error(f"Local Clustering -> Chr:{chr_num} -> cut_value:{cut_value}, Optimal_r:{optimal_r}, max_group_allele_value:{max_group_allele_value}")
                    raise
                else:
                    logger.info(f"Local Clustering -> Chr:{chr_num} -> cut_value:{cut_value}, Optimal_r:{optimal_r}, max_group_allele_value:{max_group_allele_value}")
                    break

            except:

                logger.error(f"Chr:{chr_num} louvain_nei clustering error!")
                chr_num_collapse_num_file = f"{args.output_prefix}.chr{chr_num}.utgs.uncollapse.txt"
                command = (
                    "awk -F '[, \\t]' 'NR==FNR{lines[$1]=$2;next}{if(lines[$1]<=1){print $0}}' "
                    f"{args.collapse_num_file} {utg_file} > {chr_num_collapse_num_file}"
                )
                execute_command(command, f"Failed to filter collapse Contig for {chr_num_collapse_num_file}",logger)

                chr_num_uncollapse_hic_file = f"{args.output_prefix}.chr{chr_num}.links.uncollapse.csv"
                command = (
                    "awk -F '[, \\t]' 'NR==FNR{lines[$1];next}($1 in lines && $2 in lines)' "
                    f"{chr_num_collapse_num_file} {cut_links_file} > {chr_num_uncollapse_hic_file}"
                )
                execute_command(command, f"Failed to filter hic for {chr_num_uncollapse_hic_file}",logger)
                execute_command(f"sed '1isource,target,links' -i {chr_num_uncollapse_hic_file}", f"Failed to add header to {chr_num_uncollapse_hic_file}",logger)
                
                chr_num_uncollapse_hic_cut_file = f"{args.output_prefix}.chr{chr_num}.links.uncollapse.c{float(cut_value)}.csv"
                command = (
                    f"awk -F ',' '($3> {float(cut_value)})' "
                    f"{chr_num_uncollapse_hic_file} > {chr_num_uncollapse_hic_cut_file}"
                )
                execute_command(command, f"Failed to filter hic links for {chr_num_uncollapse_hic_cut_file}",logger)

                # Adjust r and run multilevel_cluster.py
                try:
                    check_cluster_dict, max_group_allele_value = multilevel_cluster(chr_num_uncollapse_hic_cut_file, cluster_output, 1, "check", utg_file, no_expand_partig_file, args.hap_number)
                    optimal_r = multiple_adjust_r_and_cluster(
                        initial_r=1.0,
                        min_r=0.01,
                        max_r=3,
                        step=0.03,
                        cluster_output=cluster_output,
                        csv_file=chr_num_uncollapse_hic_file, 
                        utg_file=utg_file,
                        partig_file=partig_file,
                        hap_number=args.hap_number, logger=logger
                    )
                    check_cluster_dict, max_group_allele_value = multilevel_cluster(chr_num_uncollapse_hic_cut_file, cluster_output, optimal_r, "check", utg_file, no_expand_partig_file, args.hap_number)

                    # Check the number of clusters
                    num_clusters = len(check_cluster_dict)

                    if num_clusters != args.hap_number:
                        raise
                    else:
                        satisfy_clusters_list.append((check_cluster_dict, max_group_allele_value))

                    # if max_group_allele_value > min([avg_uncollapse_num/20, 1000000]):
                    if max_group_allele_value > 1000000:
                        logger.error(f"Direct Clustering -> Chr:{chr_num} -> cut_value:{cut_value}, Optimal_r:{optimal_r}, max_group_allele_value:{max_group_allele_value}")
                        raise
                    else:
                        logger.info(f"Direct Clustering -> Chr:{chr_num} -> cut_value:{cut_value}, Optimal_r:{optimal_r}, max_group_allele_value:{max_group_allele_value}")
                        break
                except:
                    cut_value += cut_value_step
                    continue

        # 未得到满足的聚类结果
        if cut_value >= cut_value_max:
            if satisfy_clusters_list:   
                sorted_satisfy_clusters_list = sorted(satisfy_clusters_list, key=lambda x: x[1])
                logger.info(f"Chr:{chr_num}: the detection is not met, take the result with the smallest detection value ({sorted_satisfy_clusters_list[0][1]}) as the clustering result.")
                final_clusters = sorted_satisfy_clusters_list[0][0]

                with open(cluster_output, 'w') as file:
                    for idx, group in enumerate(final_clusters):
                        file.write(f"group{idx+1}\t{len(final_clusters[group])}\t{' '.join(final_clusters[group])}\n")

                    
        # Check the number of clusters
        with open(cluster_output, 'r') as f:
            clusters = set(line.split(',')[0].strip() for line in f if line.strip())
        num_clusters = len(clusters)
        if num_clusters != args.hap_number:
            logger.error(f"Chr:{chr_num} clustering error!")
            return False


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
        find_best_isolated, isolated_threshold = True, 0
        min_variance_idx = louvain_reassign_allele(args.collapse_num_file, utg_rescue_file, filtered_links_file, cluster_file , args.RE_file, partig_file, args.output_prefix, find_best_isolated, no_expand_partig_file, isolated_threshold)
        logger.info(f"Chr:{chr_num} louvain_reassign_allele find best isolated : {min_variance_idx}")

        for i in range(int(args.reassign_number)-1):
            min_variance_idx = louvain_reassign_allele(args.collapse_num_file, utg_rescue_file, filtered_links_file, f"{args.output_prefix}.reassign.cluster.txt", args.RE_file, partig_file, args.output_prefix, find_best_isolated, no_expand_partig_file, isolated_threshold)

        
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
    args = parse_arguments()
    pwd = os.getcwd()
    logger = setup_logging('cluster_hap.log')
    log_start(logger, "cluster_hap.py", "1.0.0", args)

    # # # Step 1: Run partig
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

    if args.rescue:
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
            logger.info(result)



if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")