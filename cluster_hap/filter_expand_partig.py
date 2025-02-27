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


def get_N80(len_list):

    len_list.sort(reverse=True)
    total_length = sum(len_list)

    cumulative_length = 0
    N80_length = 0
    for length in len_list:
        cumulative_length += length
        if cumulative_length >= 0.95 * total_length:
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


def filter_allele(digraph, partig_dict,subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_len):

    filted_dict = defaultdict()
    len_list = [ float( ctg_RE_len[utg][1]) for utg in ctg_RE_len]
    N80 = get_N80(len_list)
    for (utg1, utg2), value in partig_dict.items():

        if utg1 not in ctg_RE_len or ctg_RE_len[utg1][1] < N80 or \
            utg2 not in ctg_RE_len or ctg_RE_len[utg2][1] < N80:
            filted_dict[tuple(sorted([utg1, utg2]))] = value

        # if not (utg1 in digraph and utg2 in digraph):
        #     continue

        if digraph.has_edge(utg1, utg2) or digraph.has_edge(utg2, utg1):
            filted_dict[tuple(sorted([utg1, utg2]))] = value

        # if utg1 in ctg_subgraph_dict and utg2 in ctg_subgraph_dict:

        #     subgraph1 = ctg_subgraph_dict[utg1]
        #     subgraph2 = ctg_subgraph_dict[utg2]

        # else:

        #     filted_dict[tuple(sorted([utg1, utg2]))] = value
        #     continue

        # if subgraph1 == subgraph2:
        #     predecessors_utg1 = set(digraph.predecessors(utg1))
        #     predecessors_utg2 = set(digraph.predecessors(utg2))

        #     successors_utg1 = set(digraph.successors(utg1))
        #     successors_utg2 = set(digraph.successors(utg2))

        #     # if not (predecessors_utg1 & predecessors_utg2 ) and \
        #     #     not (successors_utg1 & successors_utg2 ):
        #     #     filted_dict[tuple(sorted([utg1, utg2]))] = value
        # else:
        #     if utg1 not in ctg_RE_len or ctg_RE_len[utg1][1] < N80 or \
        #         utg2 not in ctg_RE_len or ctg_RE_len[utg2][1] < N80:

        #         filted_dict[tuple(sorted([utg1, utg2]))] = value
    
    save_dict = {k:v for k,v in partig_dict.items() if k not in filted_dict}

    with open("filter_partig.csv", 'w') as file:
        for (ctg1, ctg2), value in save_dict.items():
            file.write(f"{ctg1},{ctg2},{value}\n")


    return save_dict

def expand_allele(save_dict, subgraph_ctgs_dict, ctg_subgraph_dict):

    expand_dict = defaultdict()

    for (utg1, utg2) in save_dict:

        expand_dict[tuple(sorted([utg1, utg2]))] = 1

        if utg1 in ctg_subgraph_dict and utg2 in ctg_subgraph_dict and \
            ctg_subgraph_dict[utg1] != ctg_subgraph_dict[utg2]:

            subgraph1 = subgraph_ctgs_dict[ctg_subgraph_dict[utg1]]
            subgraph2 = subgraph_ctgs_dict[ctg_subgraph_dict[utg2]]

            if subgraph1 != subgraph2:
                for sub_utg1 in subgraph1:
                    for sub_utg2 in subgraph2:
                        expand_dict[tuple(sorted([sub_utg1, sub_utg2]))] = 1


    with open("merge.partig.csv", 'w') as file:
        file.write("source,target,links\n")
        for (utg1, utg2) in expand_dict:
            file.write(f"{utg1},{utg2},1\n")

    return expand_dict

# def merge_allele(allele_dict, filted_dict, subgraph_ctgs_dict, ctg_subgraph_dict):

#     re_allele_dict = defaultdict()
#     allele_filter_dict = defaultdict()
#     # print(f"{len(allele_dict)}\t{len(filted_dict)}\t{len(expand_dict)}")

#     for (utg1, utg2) in allele_dict:

#         if (utg1, utg2)  not in filted_dict:
#             allele_filter_dict[(utg1, utg2)] = 1


#     with open("filter.partig.csv", 'w') as file:
#         file.write("source,target,links\n")
#         for (utg1, utg2) in allele_filter_dict:
#             file.write(f"{utg1},{utg2},1\n")

#     expand_dict = expand_allele(allele_filter_dict, subgraph_ctgs_dict, ctg_subgraph_dict)

#     merge_dict = {**allele_filter_dict, **expand_dict}

#     for (utg1, utg2) in merge_dict:

#         if (utg1, utg2)  not in filted_dict:
#             re_allele_dict[(utg1, utg2)] = 1


#     with open("merge.partig.csv", 'w') as file:
#         file.write("source,target,links\n")
        
#         for (utg1, utg2) in re_allele_dict:
#             file.write(f"{utg1},{utg2},1\n")


#     return re_allele_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="filter partig")
    parser.add_argument('-d', '--digraph', required=True,
                        help='<filepath> GFA digraph')
    parser.add_argument('-r', '--REs', required=True,
                        help='<filepath> REs file for utgs')
    parser.add_argument('-s', '--subgraph_file', required=True,
                        help='<filepath> subgraph ctgs')
    parser.add_argument('-p', '--partig_file', required=True,
                        help='<filepath> partig file')



    # digraph_file = "digraph.csv"
    # REFile = "HiC.filtered.sort.counts_GATC.txt"
    # subgraph_file = "group_ctgs_All.txt"
    # partig_file = "rice4.partig.csv"

    args = parser.parse_args()
    digraph_file = args.digraph
    REFile = args.REs
    subgraph_file = args.subgraph_file
    partig_file = args.partig_file


    digraph = read_digraph(digraph_file)
    ctg_RE_dict = read_RE(REFile)
    partig_dict = read_partig(partig_file)


    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)

    save_dict = filter_allele(digraph, partig_dict, subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_dict)
    expand_allele(save_dict, subgraph_ctgs_dict, ctg_subgraph_dict)







