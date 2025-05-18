import igraph as ig
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import networkx as nx
import community as louvain
from collections import defaultdict
import statistics
import networkit as nk

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len

def read_Allele(Allele_file):
    allele_dict = defaultdict()
    with open(Allele_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if line[0].startswith("u") and line[1].startswith("u"):
                allele_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
    return allele_dict

def Louvain_cluster(g, id_to_node, resolution):
    partitioner = nk.community.PLM(g, refine=True, gamma=float(resolution))
    partitioner.run()
    partition = partitioner.getPartition()

    communities = defaultdict(list)
    for node in range(g.numberOfNodes()):
        cluster_id = partition.subsetOf(node)
        label = id_to_node[node]
        communities[cluster_id].append(label)

    return communities

def multilevel_cluster(csv_file, output_file, resolution=1, check=None, RE_file=None, Allele_file=None, n_hap=None):

    df = pd.read_csv(csv_file)
    df['source'] = df['source'].astype(str)
    df['target'] = df['target'].astype(str)
    df['links'] = df['links'].astype(float)

    # 建立节点映射
    g = nk.graph.Graph(weighted=True, directed=False)
    nodes = list(set(df['source']).union(set(df['target'])))
    node_to_id = {node: idx for idx, node in enumerate(nodes)}
    id_to_node = {idx: node for node, idx in node_to_id.items()}
    g.addNodes(len(nodes))

    # 添加边（避免重复添加）
    for _, row in df.iterrows():
        u = node_to_id[row['source']]
        v = node_to_id[row['target']]
        weight = row['links']
        if not g.hasEdge(u, v):
            g.addEdge(u, v, weight)

    # 社区划分
    communities = Louvain_cluster(g, id_to_node, resolution)

    if check:
        # 检查有效聚类簇数目
        # 1:阈值设置为平均聚类簇长度的 1/10
        # 2: 簇中平均contigs长度小于平均contig长度中位数的1/7 （防止聚类到核糖体）
        
        ctg_RE_len = read_REs(RE_file)
        allele_dict = read_Allele(Allele_file)
        cluster_dict = defaultdict(list)
        group_len_dict = defaultdict()
        cluster_num = len(communities)
        avg_contig_len= defaultdict()

        for idx, community in enumerate(communities.values()):
            utgs = list(community)
            utgs_len = sum([ int(ctg_RE_len[utg][1]) for utg in utgs])
            group_len_dict[idx] = utgs_len
            avg_contig_len[idx] = utgs_len/len(list(community))
        
        avg_len = sum(group_len_dict.values()) / n_hap
        # print(avg_len)
        median_avg_contig_len = statistics.median(avg_contig_len.values())

        save_group = [ k for k,v in group_len_dict.items() if v > float(avg_len)/5 and avg_contig_len[k] > median_avg_contig_len/7]

        with open(output_file, 'w') as file:
            flag = 0
            for idx, community in enumerate(communities.values()):
                if idx in save_group:
                    flag += 1
                    utgs = list(community)
                    cluster_dict[flag] = utgs
                    file.write(f"group{flag}\t{len(utgs)}\t{' '.join(utgs)}\n")
        
        # 计算聚类簇中最大同源contigs长度
        group_allele_list = list()
        for group, ctgs in cluster_dict.items():
            allele_sum = 0
            for idx1, ctg_1 in enumerate(ctgs):
                for idx2, ctg_2 in enumerate(ctgs):
                    if idx1 < idx2:
                        if tuple(sorted([ctg_1, ctg_2])) in allele_dict:
                            allele_sum += min(ctg_RE_len[ctg_1][1], ctg_RE_len[ctg_2][1])
            group_allele_list.append(allele_sum)
        # print(max(group_allele_list))
        
        return cluster_dict, max(group_allele_list)




    else:
        with open(output_file, 'w') as file:
            for idx, community in enumerate(communities.values()):
                utgs = list(community)
                file.write(f"group{idx+1}\t{len(communities[idx])}\t{' '.join(utgs)}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="multilevel clustering algorithm is used and cluster file is generated")
    parser.add_argument('-c', '--csv_file', required=True,
                        help='<filepath> csv file of the hic signal')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath> output file')
    parser.add_argument('-r', '--resolution', default=1,type=float,
                    help='<float> The resolution parameter of multilevel cluster algorithm, the larger the value, the more clusters with fewer nodes')
    parser.add_argument('--check', action='store_true', help='get the number of vaild clusters')
    parser.add_argument('--RE_file',help='File required when --check is enabled')

    parser.add_argument('--Allele_file',help='File required when --check is enabled')
    parser.add_argument('--n_hap', type=int, metavar='\b', help='Expected number of hap, required when --check is enabled')


    args = parser.parse_args()

    if args.check and args.RE_file is None:
        parser.error("--RE_file is required when --check is enabled")
    
    if args.check and args.Allele_file is None:
        parser.error("--Allele_file is required when --check is enabled")
    
    if args.check and args.n_hap is None:
        parser.error("--n_hap is required when --check is enabled")

    multilevel_cluster(args.csv_file, args.output_file, args.resolution, args.check, args.RE_file, args.Allele_file, args.n_hap)
