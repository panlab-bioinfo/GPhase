#!/usr/bin/env python3

import igraph as ig
import pandas as pd
import numpy as np
import copy as cp
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import networkx as nx
import community as louvain
from collections import defaultdict
import statistics
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import numpy as np


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
            if line[0] == "source":
                continue
            allele_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
    return allele_dict


def hac(communities, g, n_hap):

    contig_to_cluster, cluster_dict = {}, defaultdict(list) 
    for i, comm in enumerate(communities):
        for node in comm:
            contig_to_cluster[g.vs[node]['name']] = i

    cluster_ids = sorted(set(contig_to_cluster.values()))
    cluster_index = {cid: idx for idx, cid in enumerate(cluster_ids)}
    num_clusters = len(cluster_ids)
    cluster_contact = defaultdict(lambda: defaultdict(float))

    for edge in g.es:
        v1 = g.vs[edge.source]['name']
        v2 = g.vs[edge.target]['name']
        cl1 = contig_to_cluster[v1]
        cl2 = contig_to_cluster[v2]
        if cl1 != cl2:
            cluster_contact[cl1][cl2] += edge['weight']

    all_contacts = []
    for c1 in cluster_contact:
        for c2 in cluster_contact[c1]:
            all_contacts.append(cluster_contact[c1][c2])
    max_contact = max(all_contacts) if all_contacts else 1.0

    dist_mat = np.zeros((num_clusters, num_clusters))
    for i, c1 in enumerate(cluster_ids):
        for j, c2 in enumerate(cluster_ids):
            if i >= j:
                continue
            contact = cluster_contact[c1].get(c2, 1e-6)
            dist = max_contact - contact
            dist_mat[i, j] = dist
            dist_mat[j, i] = dist

    condensed_dist = squareform(dist_mat)
    Z = linkage(condensed_dist, method='average')
    target_cluster_num = n_hap
    new_labels = fcluster(Z, t=target_cluster_num, criterion='maxclust')

    cluster_to_newlabel = {
        old_cluster: new_label
        for old_cluster, new_label in zip(cluster_ids, new_labels)
    }
    new_contig_to_cluster = {
        ctg: cluster_to_newlabel[old_cl]
        for ctg, old_cl in contig_to_cluster.items()
    }

    for contig, cluster_id in new_contig_to_cluster.items():
        cluster_dict[cluster_id].append(contig)

    return cluster_dict




def multilevel_cluster(csv_file, output_file, resolution, check=None, RE_file=None, Allele_file=None, n_hap=None):

    df = pd.read_csv(csv_file, sep=',', header=0)
    if list(df.columns)[0] != "source":
        df = pd.read_csv(csv_file, sep=',')
        df.columns = ['source', 'target', 'links']
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
        # Check the number of valid clusters
        # 1: Set the threshold to 1/3 of the average cluster length
        # 2: The average contig length in the cluster is less than 1/7 of the median average contig length (to prevent clustering into ribosomes)
        
        ctg_RE_len = read_REs(RE_file)
        allele_dict = read_Allele(Allele_file)
        cluster_dict = defaultdict(list)
        group_len_dict = defaultdict()
        cluster_num = len(communities)
        avg_contig_len= defaultdict()


        # when r<1 & save_group>n_hap, try to merge clusters using HiC links
        if float(resolution) < 1 and len(communities) > n_hap:
            cluster_dict = hac(communities, g, n_hap)
            for group, utgs in cluster_dict.items():
                utgs_len = sum([ int(ctg_RE_len[utg][1]) for utg in utgs])
                group_len_dict[group] = utgs_len
                avg_contig_len[group] = utgs_len/n_hap
        else:
            for idx, community in enumerate(communities):
                utgs = [g.vs[i]['name'] for i in communities[idx]]
                utgs_len = sum([ int(ctg_RE_len[utg][1]) for utg in utgs])
                group_len_dict[idx] = utgs_len
                avg_contig_len[idx] = utgs_len/len(list(community))

        
        avg_len = sum(group_len_dict.values()) / n_hap
        median_avg_contig_len = statistics.median(avg_contig_len.values())

        save_group = [ k for k,v in group_len_dict.items() if v > float(avg_len)/5 and avg_contig_len[k] > median_avg_contig_len/7]

        if float(resolution) < 1 and len(communities) > n_hap:
            cluster_dict_cp = cp.deepcopy(cluster_dict)
            with open(output_file, 'w') as file:
                flag = 0
                for group, utgs in cluster_dict.items():
                    if group in save_group:
                        flag += 1
                        cluster_dict_cp[flag] = utgs
                        file.write(f"group{flag}\t{len(utgs)}\t{' '.join(utgs)}\n")
            cluster_dict = cluster_dict_cp
        
        else:
            with open(output_file, 'w') as file:
                flag = 0
                for idx, community in enumerate(communities):
                    if idx in save_group:
                        flag += 1
                        utgs = [g.vs[i]['name'] for i in communities[idx]]
                        cluster_dict[flag] = utgs
                        file.write(f"group{flag}\t{len(utgs)}\t{' '.join(utgs)}\n")
        
        group_allele_list = list()
        for group, ctgs in cluster_dict.items():
            allele_sum = 0
            for idx1, ctg_1 in enumerate(ctgs):
                for idx2, ctg_2 in enumerate(ctgs):
                    if idx1 < idx2:
                        if tuple(sorted([ctg_1, ctg_2])) in allele_dict:
                            allele_sum += min(ctg_RE_len[ctg_1][1], ctg_RE_len[ctg_2][1])
            group_allele_list.append(allele_sum)

        return cluster_dict, max(group_allele_list)

    else:
        with open(output_file, 'w') as file:
            for idx, community in enumerate(communities):
                utgs = [g.vs[i]['name'] for i in communities[idx]]
                file.write(f"group{idx+1}\t{len(communities[idx])}\t{' '.join(utgs)}\n")


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
