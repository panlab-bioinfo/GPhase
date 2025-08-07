import igraph as ig
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import networkx as nx
import community as louvain
from collections import defaultdict

def read_c(c):
    cluster_dict = defaultdict(list)
    subgraph_group_dict = defaultdict(list)
    with open(c, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            for utg in line[2].split():
                cluster_dict[line[0]].append(utg)
                subgraph_group_dict[utg].append(line[0])
    return cluster_dict, subgraph_group_dict

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len


def Multilevel_cluster(csv_file, output_file, resolution, check=False, RE_file=None, Allele_cluster=None, n_chr=None):

    if RE_file:
        ctg_RE_len = read_REs(RE_file)
    
    if Allele_cluster:
        cluster_dict, subgraph_group_dict = read_c(Allele_cluster)

    df = pd.read_csv(csv_file)
    df['source'] = df['source'].astype(str)
    df['target'] = df['target'].astype(str)
    df['links'] = df['links'].astype(float)
    g = ig.Graph()  
    nodes = list(set(df['source']).union(set(df['target'])))
    g.add_vertices(nodes)

    for i, row in df.iterrows():
        source = row['source']
        target = row['target']
        weight = row['links']
        g.add_edge(source, target, weight=weight)

    communities = g.community_multilevel(weights='weight',resolution=float(resolution)) 

    if check:
        chr_len_dict = defaultdict(int)
        chr_cluster_dict = defaultdict(set)
        for idx, community in enumerate(communities):
            sum_length = 0
            groups = [g.vs[i]['name'] for i in communities[idx]]
            for group in groups:
                for utg in cluster_dict[group]:
                    chr_cluster_dict[idx].add(utg)
                    if utg in ctg_RE_len:
                        sum_length += ctg_RE_len[utg][1]
            chr_len_dict[idx] = sum_length

        # 检查有效聚类簇数目
        # 1:阈值设置为平均聚类簇长度的 1/10
        threshold = sum(chr_len_dict.values()) / int(n_chr) / 2
        filtered_chr_list = [ key for key, value in chr_len_dict.items() if value > threshold ]
        
        group = 0
        with open(output_file, 'w') as file:
            for idx in chr_cluster_dict:
                if idx in filtered_chr_list:
                    group += 1
                    utgs = list(chr_cluster_dict[idx])
                    file.write(f"group{group}\t{len(utgs)}\t{' '.join(utgs)}\n")

        return len(filtered_chr_list)

    else:
        with open(output_file, 'w') as file:
            for idx, community in enumerate(communities):
                utgs = [g.vs[i]['name'] for i in communities[idx]]
                file.write(f"group{idx+1}\t{len(communities[idx])}\t{' '.join(utgs)}\n")
        return len(communities)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="multilevel clustering algorithm is used and cluster file is generated")
    parser.add_argument('-c', '--csv_file', required=True,
                        help='<filepath> csv file of the hic signal')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath> output file')
    parser.add_argument('-r', '--resolution', default=1,
                    help='<float> The resolution parameter of multilevel cluster algorithm, the larger the value, the more clusters with fewer nodes')
    parser.add_argument('--check', action='store_true', help='get the number of vaild clusters')
    parser.add_argument('--RE_file',help='File required when --check is enabled')
    parser.add_argument('--Allele_cluster',help='File required when --check is enabled')
    parser.add_argument("--n_chr", type=int, metavar='\b', help="Desired number of clusters corresponding to chromosomes, required when --check is enabled.")

    # argcomplete.autocomplete(parser)
    args = parser.parse_args()
    if args.check and args.RE_file is None:
        parser.error("--RE_file is required when --check is enabled")
    
    if args.check and args.Allele_cluster is None:
        parser.error("--Allele_cluster is required when --check is enabled")
    
    if args.check and args.n_chr is None:
        parser.error("--n_chr is required when --check is enabled")
    Multilevel_cluster(args.csv_file, args.output_file, args.resolution, args.check, args.RE_file, args.Allele_cluster, args.n_chr)
