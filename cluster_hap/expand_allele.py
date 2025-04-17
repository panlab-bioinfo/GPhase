from collections import defaultdict, deque
import csv
import networkx as nx
import argparse

#--------------------------------------#
# 过滤 等位contig对
# 非同一子图：两个contig长度大于N80
# 同一子图：相同前驱和相同后继
#--------------------------------------#



def get_N80(len_list):

        len_list.sort(reverse=True)
        total_length = sum(len_list)

        cumulative_length = 0
        N80_length = 0
        for length in len_list:
            cumulative_length += length
            if cumulative_length >= 0.8 * total_length:
                N80_length = length
                break
            
        return int(N80_length)





def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len


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

def read_allele(allele_file):
    allele_dict = defaultdict()
    with open(allele_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if line[0].startswith("u") and line[1].startswith("u"):
                allele_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
    return allele_dict


def filter_allele(digraph, utgs_list, allele_dict,subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_len):

    filted_dict = defaultdict()
    len_list = [  ctg_RE_len[utg][1] for utg in utgs_list if utg in ctg_RE_len]
    N80 = get_N80(len_list)

    for (utg1, utg2) in allele_dict:

        if digraph.has_edge(utg1, utg2) or digraph.has_edge(utg2, utg1):
            filted_dict[tuple(sorted([utg1, utg2]))] = 1

        if utg1 in ctg_subgraph_dict and utg2 in ctg_subgraph_dict:

            subgraph1 = ctg_subgraph_dict[utg1]
            subgraph2 = ctg_subgraph_dict[utg2]

        else:

            filted_dict[tuple(sorted([utg1, utg2]))] = 1
            continue

        if subgraph1 == subgraph2:
            predecessors_utg1 = set(digraph.predecessors(utg1))
            predecessors_utg2 = set(digraph.predecessors(utg2))

            successors_utg1 = set(digraph.successors(utg1))
            successors_utg2 = set(digraph.successors(utg2))

            if not (predecessors_utg1 & predecessors_utg2 ) and \
                not (successors_utg1 & successors_utg2 ):
                filted_dict[tuple(sorted([utg1, utg2]))] = 1
        else:
            if utg1 not in ctg_RE_len or ctg_RE_len[utg1][1] < N80 or \
                utg2 not in ctg_RE_len or ctg_RE_len[utg2][1] < N80:

                filted_dict[tuple(sorted([utg1, utg2]))] = 1

    return filted_dict


def expand_allele(allele_dict, subgraph_ctgs_dict, ctg_subgraph_dict):

    expand_dict = defaultdict()

    for (utg1, utg2) in allele_dict:

        
        if utg1 in ctg_subgraph_dict and utg2 in ctg_subgraph_dict and \
            ctg_subgraph_dict[utg1] != ctg_subgraph_dict[utg2]:

            subgraph1 = subgraph_ctgs_dict[ctg_subgraph_dict[utg1]]
            subgraph2 = subgraph_ctgs_dict[ctg_subgraph_dict[utg2]]

            for sub_utg1 in subgraph1:
                for sub_utg2 in subgraph2:
                    if sub_utg1 == utg1 and sub_utg2 == utg2:
                        continue
                    if sub_utg2 == utg1 and sub_utg1 == utg2:
                        continue
                    expand_dict[tuple(sorted([sub_utg1, sub_utg2]))] = 1

    return expand_dict

def merge_allele(allele_dict, filted_dict, subgraph_ctgs_dict, ctg_subgraph_dict):

    re_allele_dict = defaultdict()
    allele_filter_dict = defaultdict()

    for (utg1, utg2) in allele_dict:

        if (utg1, utg2)  not in filted_dict:
            allele_filter_dict[(utg1, utg2)] = 1


    with open("filter.partig.csv", 'w') as file:
        file.write("source,target,links\n")
        for (utg1, utg2) in allele_filter_dict:
            file.write(f"{utg1},{utg2},1\n")

    expand_dict = expand_allele(allele_filter_dict, subgraph_ctgs_dict, ctg_subgraph_dict)

    merge_dict = {**allele_filter_dict, **expand_dict}

    for (utg1, utg2) in merge_dict:

        if (utg1, utg2)  not in filted_dict:
            re_allele_dict[(utg1, utg2)] = 1


    with open("merge.partig.csv", 'w') as file:
        file.write("source,target,links\n")
        
        for (utg1, utg2) in re_allele_dict:
            file.write(f"{utg1},{utg2},1\n")


    return re_allele_dict

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="expand and filter allele pair use subgraph")
    parser.add_argument('-g', '--gfa', required=True,
                        help='<filepath> gfa file')
    parser.add_argument('-s', '--subgraph', required=True,
                        help='<filepath>group_ctgs_All.txt')
    parser.add_argument('-r', '--REs', required=True,
                        help='<filepath>REs and length for utgs')
    parser.add_argument('-a', '--allele', required=True,
                        help='<filepath>allele utgs pair for chromosome')


    args = parser.parse_args()

    REFile = args.REs
    gfa_filePath = args.gfa
    subgraph_file = args.subgraph
    allele_file = args.allele


    # gfa_filePath = "waxapple.asm.bp.p_utg.gfa"
    # subgraph_file = "group_ctgs_merge.txt"
    # # allele_file = "wa.partig.19.95.trans.txt"
    # # allele_file = "wa.partig.19.trans.txt"
    # allele_file = "wa.partig.21_0.5.trans.txt"
    # REFile = "alignment.counts_GATC.txt"

    ctg_RE_len = read_REs(REFile)
    graph, utgs_list = read_gfa(gfa_filePath)
    digraph = run_trans_digraph(graph)

    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)
    allele_dict = read_allele(allele_file)

    filted_dict = filter_allele(digraph, utgs_list, allele_dict,subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_len)
    expand_dict = expand_allele(allele_dict, subgraph_ctgs_dict, ctg_subgraph_dict)
    re_allele_dict = merge_allele(allele_dict, filted_dict, subgraph_ctgs_dict, ctg_subgraph_dict)







                    









