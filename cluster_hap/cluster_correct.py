from matplotlib.colors import to_hex
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



def read_collapse_num(collapse_num_file):

    collapse_num_dict = defaultdict()
    with open(collapse_num_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0].startswith("utg") or line[0].startswith("utig") :
                try:
                    collapse_num_dict[line[0]] = int(line[1])
                except:
                    print(f"error : {line[0]}\t{line[1]}")
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


def find_other_group(hic_links_dict, cluster_modify_dict, min_utg):

    links_dict = defaultdict(float)

    for group in cluster_modify_dict:

        if min_utg in cluster_modify_dict[group]:
            continue

        for utg2 in cluster_modify_dict[group]:
            if tuple(sorted([min_utg, utg2])) in hic_links_dict:
                links_dict[group] += hic_links_dict[ tuple(sorted([min_utg, utg2]))]

        links_dict[group] /= len(cluster_modify_dict[group])

    links_sorted_dict = dict(sorted(links_dict.items(), key=lambda item: item[1], reverse=True))

    return links_sorted_dict

def run_cor_cluster(collapse_num_dict, hic_links_dict, cluster_dict, utg_group_dict, allele_key_dict, utgs_list, n_hap):

    cluster_modify_dict = copy.deepcopy(cluster_dict)
    cluster_copy_dict = defaultdict()


    cluster_modify_flag_dict = { utg:[ "group"+str(i) for i in range(1,int(n_hap)+1)] for utg in utgs_list if collapse_num_dict[utg] <2 }


    cluster_copy_dict = copy.deepcopy(cluster_modify_dict)

    
    for group_ in list(cluster_modify_dict.keys()):  
        Allele_utg_pair = set(tuple(sorted([utg1, utg2]))
                            for utg1 in cluster_modify_dict[group_]
                            for utg2 in cluster_modify_dict[group_]
                            if utg1 != utg2 and tuple(sorted([utg1, utg2])) in allele_key_dict)
        Allele_utgs = set(utg for pair in Allele_utg_pair for utg in pair)
        unAllele_utgs = set(utg for utg in cluster_modify_dict[group_] if utg not in Allele_utgs)

        for (utg1, utg2) in Allele_utg_pair:
            if utg1 in utg_group_dict and utg2 in utg_group_dict and utg1 in cluster_modify_dict[group_] and utg2 in cluster_modify_dict[group_]:
                utg1_num = cal_hic_links(hic_links_dict, utg1, unAllele_utgs)
                utg2_num = cal_hic_links(hic_links_dict, utg2, unAllele_utgs)
                min_utg = utg1 if utg1_num < utg2_num else utg2

                try:
                    select_group_dict = find_other_group(hic_links_dict, cluster_modify_dict, min_utg)
                    if select_group_dict:
                        select_group_list = list(select_group_dict)
                        select_group = next((x for x in select_group_list if x in cluster_modify_flag_dict[min_utg]), None)
                        if select_group:
                            cluster_modify_dict[group_].remove(min_utg)
                            cluster_modify_dict[select_group].append(min_utg)
                            cluster_modify_flag_dict[min_utg].remove(group_)
                except:
                    print(f"error:{min_utg} not in cluster.txt {group_}\n")

    with open("cor.cluster.txt", "w") as file:
        for group in cluster_modify_dict:
            file.write(f"{group}\t{len(cluster_modify_dict[group])}\t")
            for utg in cluster_modify_dict[group]:
                file.write(f"{utg} ")
            file.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="correct the cluster base Allele file")
    parser.add_argument('-cn', '--collapse_num', required=True,
                        help='<filepath> collapse num file ')
    parser.add_argument('-chr', '--chr_file', required=True,
                        help='<filepath> chr file')
    parser.add_argument('-l', '--hic_liks', required=True,
                        help='<filepath> hic links file')
    parser.add_argument('-a', '--allele', required=True,
                        help='<filepath> allele file')
    parser.add_argument('-c', '--clusters', required=True,
                        help='<filepath> clusters file')
    parser.add_argument('-n_hap', '--n_hap', required=True,
                        help='<filepath> hap number')
            
    args = parser.parse_args()

    collapse_num_dict = read_collapse_num(args.collapse_num)
    utgs_list = read_chr_utgs(args.chr_file)
    hic_links_dict, hic_nei_dict = read_l(args.hic_liks)
    allele_utg_dict, allele_key_dict = read_allele(args.allele)
    cluster_dict, utg_group_dict = read_c(args.clusters)

    run_cor_cluster(collapse_num_dict, hic_links_dict, cluster_dict, utg_group_dict, allele_key_dict, utgs_list, args.n_hap)