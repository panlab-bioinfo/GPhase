from collections import defaultdict
import argparse
import logging
import os
import os.path as op
import sys
import re
import numpy as np
import pandas as pd
import itertools

def read_clm(clm_file):
    utg_utg_pos_dict = defaultdict(list)
    with open(clm_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            utg_utg_pos_dict[tuple(sorted([line[0],line[1]]))] = line[3:]
    return utg_utg_pos_dict

def read_cluster(cluster_file):
    group_utgs_dict = defaultdict(list)
    with open(cluster_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            group_utgs_dict[line[0]] = line[2:]
    return group_utgs_dict


def main(args):
    parser = argparse.ArgumentParser(description='Classify clm file based on the input cluster files')
    parser.add_argument('clm', help='clm file')
    parser.add_argument('cluster', help='cluster file')

    args = parser.parse_args(args)

    clm_file = args.clm 
    cluster_file = args.cluster

    utg_utg_pos_dict = read_clm(clm_file)
    group_utgs_dict = read_cluster(cluster_file)

    idr = ['+', '-']
    prefix = cluster_file[:-7]
    os.makedirs("split_clms", exist_ok=True)
    for group, utgs in group_utgs_dict.items():
        group_file = open(f"split_clms/{group}.clm", 'w')
        utg_pairs = [tuple([sorted([utg1,utg2])]) for utg1 in utgs for utg2 in utgs if utg1 != utg2]
        for pair in utg_pairs:
            pair_idr_4 = [ str(utg)+i for utg in list(pair)[0] for i in idr]
            pair_idr_4_list = [tuple(sorted([utg1, utg2])) for utg1 in pair_idr_4 for utg2 in pair_idr_4 if utg1[:-1] != utg2[:-1]]
            # print(pair_idr_4_list)
            for pair_idr in pair_idr_4_list:
                if pair_idr in utg_utg_pos_dict:
                    group_file.write(f"{pair_idr[0]} {pair_idr[1]}\t{len(utg_utg_pos_dict[pair_idr])}\t{' '.join(utg_utg_pos_dict[pair_idr])}\n")
        group_file.close()


main(sys.argv[1:])
