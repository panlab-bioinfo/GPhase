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
    G = nx.from_pandas_edgelist(df, source='source', target='target', create_using=nx.DiGraph())
    
    return G

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
            elif(line[0] == "L"): # no self loop
                # gfa.add_edge(line[1],line[3])
                graph.add_edge(line[1],line[3],strand1=line[2], strand2=line[4], match=int(line[5][:-1]))
                if line[2] == '+' :
                    graph.nodes[line[1]]['out'].add(line[3]) 
                elif line[2] == '-':
                    graph.nodes[line[1]]['enter'].add(line[3])


    return graph, utgs_list

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

def read_RE(REFile):
    ctg_RE_dict = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_dict

def get_All_N80(utgs_list, ctg_RE_dict):

    len_list = [int(ctg_RE_dict[utg][1]) for utg in utgs_list if utg in ctg_RE_dict]
    return get_N80(len_list)



def trans_net(utgs_list, graph, ctg_RE_dict):

    all_N80 = get_All_N80(utgs_list, ctg_RE_dict)
    ctg_group_dict = defaultdict(int)
    group_ctg_dict = defaultdict(lambda: dict)
    group_ctg_N80_dict = defaultdict(lambda: dict)
    group_ctg_rm_dict = defaultdict(lambda: dict)

    for idx,component in enumerate(nx.connected_components(graph), 1):

        group_list = list(component)

        if all(int(ctg_RE_dict[ctg][1]) < all_N80 for ctg in group_list if ctg in ctg_RE_dict):
            group_ctg_rm_dict[idx] = {'ctgs':group_list,'length':1000000}
            continue
        

        for ctg in group_list:
            ctg_group_dict[ctg] = idx

        group_len = sum([int(ctg_RE_dict[ctg][1]) for ctg in component if ctg in ctg_RE_dict])
        group_ctg_dict[idx] = {'ctgs':group_list, 'length':group_len}

    len_list = [group_ctg_dict[idx]['length'] for idx in group_ctg_dict]
    N80_length = get_N80(len_list)
    group_ctg_N80_dict = {group:dict_ for group, dict_ in group_ctg_dict.items() if dict_['length'] > N80_length}

    with open("group_ctgs_N80.txt", 'w') as file:
        for group in group_ctg_N80_dict:
            file.write(f"{group}\t{group_ctg_dict[group]['length']}\t{','.join(list(group_ctg_N80_dict[group]['ctgs']))}\n")

    with open("group_ctgs_All.txt", 'w') as file:
        for group in group_ctg_dict:
            file.write(f"{group}\t{group_ctg_dict[group]['length']}\t{','.join(list(group_ctg_dict[group]['ctgs']))}\n")

    with open("group_ctgs_Rm.txt", 'w') as file:
        for group in group_ctg_rm_dict:
            file.write(f"{group}\t{group_ctg_rm_dict[idx]['length']}\t{','.join(list(group_ctg_rm_dict[group]['ctgs']))}\n")
    
    return group_ctg_dict, group_ctg_N80_dict, ctg_group_dict


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
    len_list = [ float( ctg_RE_len[utg][1]) for utg in utgs_list if utg in ctg_RE_len]
    N80 = get_N80(len_list)

    for (utg1, utg2), value in allele_dict.items():

        if not (utg1 in digraph and utg2 in digraph):
            filted_dict[tuple(sorted([utg1, utg2]))] = value
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

        with open("filter_allele.csv", 'w') as file:
            for (ctg1, ctg2), value in filted_dict.items():
                file.write(f"{ctg1},{ctg2},{value}\n")

    return filted_dict



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="filter partig")
    parser.add_argument('-d', '--digraph_file', required=True,
                        help='<filepath>clusters file')
    parser.add_argument('-g', '--gfa_filePath', required=True,
                        help='<filepath>subgraphs file')
    parser.add_argument('-r', '--REFile', required=True,
                        help='<filepath>subgraphs file')
    parser.add_argument('-p', '--partig_file', required=True,
                        help='<filepath>partig file')

    # digraph_file = "c88.split.digraph.csv"
    # gfa_filePath = "c88.bp.p_utg.noseq.gfa"
    # REFile = "c88.counts_RE.txt"
    subgraph_file = "group_ctgs_All.txt"
    # allele_file = "c88.partig.csv"

    args = parser.parse_args()
    digraph_file = args.digraph_file
    gfa_filePath = args.gfa_filePath
    REFile = args.REFile
    allele_file = args.partig_file



    digraph = read_digraph(digraph_file)
    graph, utgs_list = read_gfa(gfa_filePath)
    ctg_RE_dict = read_RE(REFile)
    allele_dict = read_allele(allele_file)

    trans_net(utgs_list, graph, ctg_RE_dict)

    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)

    filter_allele(digraph, utgs_list, allele_dict,subgraph_ctgs_dict, ctg_subgraph_dict, ctg_RE_dict)







