from collections import defaultdict, deque
import csv
import networkx as nx
import pandas as pd
import argparse


#--------------------------------------#
# 过滤 等位contig对
# 非同一子图：两个contig长度大于N80
# 同一子图：相同前驱和相同后继
#--------------------------------------#

def read_digraph(digraph_file):

    df = pd.read_csv(digraph_file)

    # Create a directed graph and load the edge list with weights
    digraph = nx.from_pandas_edgelist(df, source='source', target='target', create_using=nx.DiGraph())
    
    return digraph


def get_N80(len_list, threshold = 0.8):

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

def read_RE(REFile):
    ctg_RE_dict = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_dict








def read_subgraph(subgraph_file):

    subgraph_ctgs_dict = defaultdict(list)
    ctg_subgraph_dict = defaultdict()
    with open(subgraph_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            subgraph_num = line[0]
            ctgs = list(line[2].split(','))
            subgraph_ctgs_dict[subgraph_num] = ctgs
            for ctg in ctgs:
                ctg_subgraph_dict[ctg] = subgraph_num

    return subgraph_ctgs_dict, ctg_subgraph_dict

def read_partig(partig_file):
    partig_dict = defaultdict()
    with open(partig_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if line[0].startswith("u") and line[1].startswith("u"):
                partig_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
    return partig_dict


def filter_allele(digraph, partig_dict,subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_len, threshold = 0.8):

    filted_dict = defaultdict()
    len_list = [ float( ctg_RE_len[utg][1]) for utg in ctg_RE_len]
    N80 = get_N80(len_list, threshold=threshold)

    for (utg1, utg2), value in partig_dict.items():

        if not (utg1 in digraph and utg2 in digraph):
            # filted_dict[tuple(sorted([utg1, utg2]))] = value
            continue

        if digraph.has_edge(utg1, utg2) or digraph.has_edge(utg2, utg1):
            filted_dict[tuple(sorted([utg1, utg2]))] = value

        if utg1 in ctg_subgraph_dict and utg2 in ctg_subgraph_dict:

            subgraph1 = ctg_subgraph_dict[utg1]
            subgraph2 = ctg_subgraph_dict[utg2]

        else:

            filted_dict[tuple(sorted([utg1, utg2]))] = value
            continue

        if subgraph1 == subgraph2:
            predecessors_utg1 = set(digraph.predecessors(utg1))
            predecessors_utg2 = set(digraph.predecessors(utg2))

            successors_utg1 = set(digraph.successors(utg1))
            successors_utg2 = set(digraph.successors(utg2))

            if not (predecessors_utg1 & predecessors_utg2 ) and \
                not (successors_utg1 & successors_utg2 ):
                filted_dict[tuple(sorted([utg1, utg2]))] = value
        else:
            if utg1 not in ctg_RE_len or ctg_RE_len[utg1][1] < N80 or \
                utg2 not in ctg_RE_len or ctg_RE_len[utg2][1] < N80:

                filted_dict[tuple(sorted([utg1, utg2]))] = value

        with open("filter_partig.csv", 'w') as file:
            for (ctg1, ctg2), value in partig_dict.items():
                if (ctg1, ctg2) not in filted_dict:
                    file.write(f"{ctg1},{ctg2},{value}\n")

    return filted_dict


if __name__ == '__main__':

    def validate_n_threshold(value=0.8):
        value = float(value)
        if value <= 0 or value >= 1:
            raise argparse.ArgumentTypeError("num_N must be between 0 and 1 (exclusive).")
        return value

    parser = argparse.ArgumentParser(description="filter partig")
    parser.add_argument('-d', '--digraph', required=True,
                        help='<filepath> GFA digraph')
    parser.add_argument('-r', '--REs', required=True,
                        help='<filepath> REs file for utgs')
    parser.add_argument('-s', '--subgraph_file', required=True,
                        help='<filepath> subgraph ctgs')
    parser.add_argument('-p', '--partig_file', required=True,
                        help='<filepath> partig file')
    parser.add_argument('-n', '--n_threshold',  default=0.8,type=validate_n_threshold,
                        help='<filepath>  filter utg(default: N80)')



    # digraph_file = "digraph.csv"
    # REFile = "HiC.filtered.sort.counts_GATC.txt"
    # subgraph_file = "group_ctgs_All.txt"
    # partig_file = "rice4.partig.csv"

    args = parser.parse_args()
    digraph_file = args.digraph
    REFile = args.REs
    subgraph_file = args.subgraph_file
    partig_file = args.partig_file
    n_threshold = args.n_threshold


    digraph = read_digraph(digraph_file)
    ctg_RE_dict = read_RE(REFile)
    partig_dict = read_partig(partig_file)


    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)

    filter_allele(digraph, partig_dict, subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_dict, n_threshold)






