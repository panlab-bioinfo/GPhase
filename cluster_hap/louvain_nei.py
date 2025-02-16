
from collections import defaultdict
import matplotlib.pyplot as plt
import networkx as nx
import igraph as ig
import pandas as pd
import numpy as np
import subprocess
import argparse
import copy
import os
import sys



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

def read_chr_utgs(chr_file):
    utgs_list = list()
    with open(chr_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0].startswith("utg") or line[0].startswith("utig"):
                utgs_list.append(line[0])
    return utgs_list


def read_c(c):
    cluster_dict = defaultdict(list)
    utg_group_dict = defaultdict(list)
    with open(c, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            for utg in line[2].split():
                cluster_dict[line[0]].append(utg)
                utg_group_dict[utg].append(line[0])
    return cluster_dict, utg_group_dict


def read_l(l):
    hic_nei_dict = defaultdict(set)
    hic_links_dict = defaultdict()
    with open(l, 'r') as file:
        for line in file:
            if line.startswith("utg") or line.startswith("utig") :
                line = line.strip().split(',')
                hic_links_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
                hic_nei_dict[line[0]].add(line[1])
                hic_nei_dict[line[1]].add(line[0])
    return hic_links_dict, hic_nei_dict

def read_allele(allele_file):

    allele_utg_dict = defaultdict(set)
    allele_key_dict = defaultdict()
    with open(allele_file, 'r') as file:
        for line in file:
            if line.startswith("utg") or line.startswith("utig"):
                line = line.strip().split(',')
                allele_utg_dict[line[0]].add(line[1])
                allele_key_dict[tuple(sorted([line[0], line[1]]))] = 1
    return allele_utg_dict, allele_key_dict

def cal_hic_links(hic_links_dict, utg, utgs_list):

    hic_num = 0
    for utg2 in utgs_list:
        tmp = tuple(sorted([utg, utg2]))
        if tmp in hic_links_dict:
            hic_num += hic_links_dict[tmp]
    return  hic_num



def run(collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, allele_utg_dict, allele_key_dict):

    # collapse_num_file = "06.genes.round.cn"
    # chr_file = "chr06.txt"
    # l = "rice4.links.nor.csv"
    # allele_file =  "chr06.allel.csv"

    # collapse_num_dict = read_collapse_num(collapse_num_file)
    # utgs_list = read_chr_utgs(chr_file)
    # hic_links_dict, hic_nei_dict = read_l(l)
    # allele_utg_dict, allele_key_dict = read_allele(allele_file)

    for utg in utgs_list:
        if utg not in collapse_num_dict:
            collapse_num_dict[utg] = int(1)

    louvin_nei_log_file = open('louvain_nei.log', 'w')
    uncollapse_list = [ utg for utg in utgs_list \
        if collapse_num_dict[utg] \
        and collapse_num_dict[utg] < 2]

    chr_utg_links = { pair:value  for pair, value in hic_links_dict.items() if pair[0] in utgs_list and pair[1] in utgs_list }

    utg_group_links_dict = defaultdict(lambda : defaultdict(int))

    re_dict = defaultdict(lambda: defaultdict(int))

    G = nx.Graph()

    os.makedirs("utgs_cluster", exist_ok=True)

    script_path = os.path.abspath(sys.path[0])
    script_path_add = os.path.join(script_path,"multilevel_cluster.py")

    for utg in utgs_list:
        # print(utg)

        nei_utg = list(hic_nei_dict[utg])
        nei_uncollapse_utg = list(set(nei_utg) & set(uncollapse_list))

        # 对坍缩节点进行邻居非坍缩节点聚类
        if collapse_num_dict[utg] > 1 and nei_utg and len(nei_uncollapse_utg) > 3 :

            # 获取 nei_uncollapse_utg 之间信号
            nei_uncollapse_utg_links = { pair:value  for pair, value in hic_links_dict.items() if pair[0] in nei_uncollapse_utg and pair[1] in nei_uncollapse_utg }
            with open(f"utgs_cluster/{utg}.links.nor.csv", 'w') as file:
                file.write("source,target,links\n")
                for pair,value in nei_uncollapse_utg_links.items():
                    file.write(f"{pair[0]},{pair[1]},{value}\n")

            try:
                script = f"""python {script_path_add} -c utgs_cluster/{utg}.links.nor.csv -o utgs_cluster/{utg}.cluster.txt"""
                process = subprocess.run(script, shell=True, executable='/bin/bash', capture_output=True, text=True)

                cluster_file = f"utgs_cluster/{utg}.cluster.txt"
                cluster_dict, utg_group_dict = read_c(cluster_file)
            except:
                louvin_nei_log_file.write(f"error: {utg}\t{','.join(nei_utg)}\t{','.join(nei_uncollapse_utg)}\n")
                cluster_dict = defaultdict()
                # print(f"error: {utg}\t{','.join(nei_utg)}\t{','.join(nei_uncollapse_utg)}")

            
            # cluster filter : 对聚类结果进行过滤，过滤掉等位的基因组（和不等位的基因links小的过滤掉）
            if not nei_uncollapse_utg_links or not cluster_dict:
                continue
            for group_ in cluster_dict:
                Allele_utg_pair = set( tuple(sorted([utg1, utg2])) 
                                        for utg1 in cluster_dict[group_] 
                                        for utg2 in cluster_dict[group_] 
                                        if utg1 != utg2 and tuple(sorted([utg1, utg2])) in allele_key_dict 
                                    )
                Allele_utgs = set( utg for pair in Allele_utg_pair for utg in pair )
                unAllele_utgs = set(utg for utg in cluster_dict[group_] if utg not in Allele_utgs)

                for (utg1, utg2) in Allele_utg_pair:
                    if utg1 in cluster_dict[group_] and utg2 in cluster_dict[group_]:

                        utg1_num = cal_hic_links(hic_links_dict, utg1, unAllele_utgs)
                        utg2_num = cal_hic_links(hic_links_dict, utg2, unAllele_utgs)

                        min_utg = utg1 if utg1_num < utg2_num else utg2
                        try:
                            cluster_dict[group_].remove(min_utg)
                        except:
                            louvin_nei_log_file.write(f"error:{min_utg} not in {utg}.cluster.txt {group_}\n")
                            # print(f"error:{min_utg} not in {utg}.cluster.txt {group_}")
                louvin_nei_log_file.write(f"{utg}\t{group_}\t{','.join(cluster_dict[group_])}\n")
                # print(f"{utg}\t{group_}\t{','.join(cluster_dict[group_])}")

            for group_ in cluster_dict:
                for idx1, utg_1 in enumerate(cluster_dict[group_]):
                    for idx2, utg_2 in enumerate(cluster_dict[group_]):
                        if idx1 < idx2 :
                            if G.has_edge(utg_1, utg_2):
                                G[utg_1][utg_2]['links'] += 1  
                            else:
                                G.add_edge(utg_1, utg_2, links=1)

    edges = list(G.edges(data=True))

    # 转换为 DataFrame
    df_edges = pd.DataFrame(edges, columns=['source', 'target', 'links'])
    if not df_edges['links'].empty:
        if df_edges['links'].iloc[0] is not None:
            attributes_df = pd.json_normalize(df_edges['links'])
            df_edges = pd.concat([df_edges[['source', 'target']], attributes_df], axis=1)

    # 保存为 CSV
    df_edges.to_csv('louvain_nei.csv', index=False)
    louvin_nei_log_file.close()





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="louvain nei cluster")
    parser.add_argument('-c', '--collase_num', required=True,
                        help='<filepath> collapse num for utgs')
    parser.add_argument('-chr', '--chromosome', required=True,
                        help='<filepath>all utg for chrompsome')
    parser.add_argument('-l', '--links', required=True,
                        help='<filepath>hic links for utgs')
    parser.add_argument('-a', '--allele', required=True,
                        help='<filepath>allele utgs pair for chromosome')

    # argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # collapse_num_file = "06.genes.round.cn"
    # chr_file = "chr06.txt"
    # l = "rice4.links.nor.csv"
    # allele_file =  "chr06.allel.csv"

    collapse_num_file = args.collase_num
    chr_file = args.chromosome
    l = args.links
    allele_file = args.allele

    collapse_num_dict = read_collapse_num(collapse_num_file)
    utgs_list = read_chr_utgs(chr_file)
    hic_links_dict, hic_nei_dict = read_l(l)
    allele_utg_dict, allele_key_dict = read_allele(allele_file)

    run(collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, allele_utg_dict, allele_key_dict)
        



            






        






