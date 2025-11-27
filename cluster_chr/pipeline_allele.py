#!/usr/bin/env python3

from collections import defaultdict,Counter
from multilevel_cluster_v2 import Multilevel_cluster
import csv
import networkx as nx
import subprocess
import argparse
import os 
import sys

def read_gfa(gfa_filePath):

    utgs_list = list()
    # graph = nx.MultiDiGraph()
    graph = nx.Graph()
    with open(gfa_filePath, 'r') as file:
        for line in file:
            line = line.strip().split()
            if(line[0] == "S"):
                utgs_list.append(line[1])
                graph.add_node(line[1], length = len(line[2]), seq = line[2], coverage = int(line[4][5:]), visit_count = 1, out=set(), enter=set(), type=None)
            elif(line[0] == "L"): 
                graph.add_edge(line[1],line[3],strand1=line[2], strand2=line[4], match=int(line[5][:-1]))
                if line[2] == '+' :
                    graph.nodes[line[1]]['out'].add(line[3]) 
                elif line[2] == '-':
                    graph.nodes[line[1]]['enter'].add(line[3])


    return graph, utgs_list

def read_RE(REFile):
    ctg_RE_dict = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (line[1], line[2])
    return ctg_RE_dict

def read_l(l):
    hic_nei_dict = defaultdict(set)
    hic_links_dict = {}

    with open(l, newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            if not row:
                continue
            if row[0].startswith(('utg', 'utig')):
                a, b, w = row[0], row[1], float(row[2])

                key = (a, b) if a < b else (b, a)

                hic_links_dict[key] = w
                hic_nei_dict[a].add(b)
                hic_nei_dict[b].add(a)

    return hic_links_dict, hic_nei_dict

def read_allele(allele_file):
    allele_dict = defaultdict()
    with open(allele_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if line[0].startswith("u") and line[1].startswith("u"):
                allele_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])

    return allele_dict

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

def get_N80(len_list, threshold=0.8):

    len_list.sort(reverse=True)
    total_length = sum(len_list)

    cumulative_length = 0
    N80_length = 0
    for length in len_list:
        cumulative_length += length
        if cumulative_length >= threshold * total_length:
            N80_length = length
            break

    return int(N80_length)

def get_All_N80(utgs_list, ctg_RE_dict):

    len_list = [int(ctg_RE_dict[utg][1]) for utg in utgs_list if utg in ctg_RE_dict]
    return get_N80(len_list, threshold=0.9)
    


def trans_net(utgs_list, graph, ctg_RE_dict):

    all_N80 = get_All_N80(utgs_list, ctg_RE_dict)
    ctg_group_dict = defaultdict(int)
    group_ctg_dict = defaultdict(lambda: dict)
    group_ctg_N80_dict = defaultdict(lambda: dict)
    group_ctg_filter_dict = defaultdict(lambda: dict)

    for idx,component in enumerate(nx.connected_components(graph), 1):

        group_list = list(component)
        # all node < N80
        if all(int(ctg_RE_dict[ctg][1]) < all_N80 for ctg in group_list if ctg in ctg_RE_dict):
            group_ctg_filter_dict[idx] = {'ctgs':group_list,'length':1}
            continue
        

        for ctg in group_list:
            ctg_group_dict[ctg] = idx

        group_len = sum([int(ctg_RE_dict[ctg][1]) for ctg in component if ctg in ctg_RE_dict])
        group_ctg_dict[idx] = {'ctgs':group_list, 'length':group_len}

    len_list = [group_ctg_dict[idx]['length'] for idx in group_ctg_dict]

    # N80_length = get_N80(len_list, threshold=1)
    # print(f"N80: {N80_length}")
    # group_ctg_N80_dict = {group:dict_ for group, dict_ in group_ctg_dict.items() if dict_['length'] > N80_length}

    # group_ctg_N80_dict = {group:dict_ for group, dict_ in group_ctg_dict.items() if dict_['length'] > 0}

    # with open("group_ctgs_N80.txt", 'w') as file:
    #     for group in group_ctg_N80_dict:
    #         file.write(f"{group}\t{group_ctg_dict[group]['length']}\t{','.join(list(group_ctg_N80_dict[group]['ctgs']))}\n")

    with open("group_ctgs_save.txt", 'w') as file:
        for group in group_ctg_dict:
            file.write(f"{group}\t{group_ctg_dict[group]['length']}\t{','.join(list(group_ctg_dict[group]['ctgs']))}\n")

    with open("group_ctgs_filter.txt", 'w') as file:
        for group in group_ctg_filter_dict:
            file.write(f"{group}\t{group_ctg_filter_dict[group]['length']}\t{','.join(list(group_ctg_filter_dict[group]['ctgs']))}\n")

    with open("group_ctgs_All.txt", 'w') as f:
        subprocess.run(['cat', 'group_ctgs_save.txt', 'group_ctgs_filter.txt'], stdout=f, text=True)
    
    return group_ctg_dict, ctg_group_dict


def cluster(group_ctg_dict, ctg_group_dict, hic_links_dict, allele_dict):
    
    group_links_dict = defaultdict(float)
    group_allele_dict = defaultdict(float)

    for pair, value in hic_links_dict.items():
        (ctg1, ctg2) = pair
        group_links_dict[ tuple(sorted([ctg_group_dict[ctg1], ctg_group_dict[ctg2]]))] += value if value > 0 else 0

    for (ctg1, ctg2) in allele_dict:
        group_allele_dict[ tuple(sorted([ctg_group_dict[ctg1], ctg_group_dict[ctg2]]))] += float(allele_dict[tuple(sorted([ctg1,ctg2]))])
        


    with open('group.links.save.csv', 'w') as file:

        file.write("source,target,links\n")

        for ctg1, ctg2 in group_links_dict:

            if  ctg1 != ctg2 and ctg1 in group_ctg_dict and ctg2 in group_ctg_dict and group_links_dict[tuple(sorted([ctg1, ctg2]))] > 0:
                file.write(f"{ctg1},{ctg2},{group_links_dict[tuple(sorted([ctg1, ctg2]))]}\n")

    
    with open('group.allele.csv', 'w') as file:
        file.write("source,target,links\n")
        for ctg1, ctg2 in group_allele_dict:
            if  ctg1 != ctg2 and ctg1 in group_ctg_dict and ctg2 in group_ctg_dict:
                file.write(f"{ctg1},{ctg2},{group_allele_dict[tuple(sorted([ctg1, ctg2]))]}\n")
    

def expand_cluster(csv_file, group_ctg_dict, output_prefix, n_chr):

    cluster_file = f"{output_prefix}.allele.cluster.txt"
    with open(csv_file, 'r') as f:
        row_count = sum(1 for _ in f) - 1 

    if row_count > n_chr:
        Multilevel_cluster(csv_file, cluster_file, 1, False, None, None, None)
        cluster_dict, utg_group_dict = read_c(cluster_file)
    else:
        cluster_dict = defaultdict(list)
        utg_group_dict = defaultdict(int)
        for idx, group in enumerate(group_ctg_dict, 1):
            cluster_dict["group"+str(idx)].append(str(group))
            utg_group_dict[group] = "group"+str(idx)

    # expand
    len_iter = len(cluster_dict)
    for group in group_ctg_dict:
        if group not in utg_group_dict:
            cluster_dict["group"+str(len_iter+1)].append(str(group))
            len_iter += 1
    

    with open(f"{output_prefix}.allele.cluster.expand.txt", 'w') as file:
        for group in cluster_dict:
            file.write(f"{group}\t{len(cluster_dict[group])}\t{' '.join(cluster_dict[group])}\n")


def Pipeline_allele(split_gfa_file, HiC_file, RE_file, partig_file, output_prefix, n_chr):

    graph, utgs_list = read_gfa(split_gfa_file)
    ctg_RE_dict = read_RE(RE_file)
    hic_links_dict, hic_nei_dict = read_l(HiC_file)
    allele_dict = read_allele(partig_file)

    group_ctg_dict, ctg_group_dict = trans_net(utgs_list, graph, ctg_RE_dict)
    cluster(group_ctg_dict, ctg_group_dict, hic_links_dict, allele_dict)

    csv_file = 'group.allele.csv'
    expand_cluster(csv_file, group_ctg_dict,output_prefix, n_chr)


    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="expand and filter allele pair use subgraph")
    parser.add_argument('-g', '--gfa', required=True,
                        help='<filepath> gfa file')
    parser.add_argument('-l', '--links', required=True,
                        help='<filepath>hic links for ctgs')
    parser.add_argument('-r', '--REs', required=True,
                        help='<filepath>REs and length for utgs')
    parser.add_argument('-a', '--allele', required=True,
                        help='<filepath>allele utgs pair for chromosome')
    parser.add_argument('--n_chr', required=True,
                        help='<int>Desired number of clusters corresponding to chromosomes.')
    parser.add_argument('-n', '--name', required=True,
                        help='<str> output prefix')



    args = parser.parse_args()

    REFile = args.REs
    gfa_filePath = args.gfa
    l = args.links
    allele_file = args.allele
    output_prefix = args.name
    n_chr = args.n_chr

    Pipeline_allele(gfa_filePath, l, REFile, allele_file, output_prefix, n_chr)