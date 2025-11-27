#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import networkx as nx
from math import log
import argparse

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len

def read_fai(fai_file):
    fai_dict = defaultdict()
    fai_reverse_dict = defaultdict(str)

    with open(fai_file, 'r') as file:
        idx = 1
        for line in file:
            line = line.strip().split()
            fai_dict[str(idx)] = line[0]
            fai_reverse_dict[line[0]] = idx
            idx += 1

    return fai_dict, fai_reverse_dict

def add_graph_allele(digraph_file):

    graph_allele_dict = defaultdict(tuple)

    df = pd.read_csv(digraph_file)
    G = nx.DiGraph() 
    for _, row in df.iterrows():
        G.add_edge(row['source'], row['target'])

    target_nodes = [n for n in G.nodes if G.out_degree(n) == 1 and G.in_degree(n) == 1]

    for utg1 in target_nodes:
        for utg2 in target_nodes:
            
            if utg1 >= utg2:
                continue
            utg1_predecessors = set(G.predecessors(utg1))
            utg1_successors = set(G.successors(utg1))
            utg2_predecessors = set(G.predecessors(utg2))
            utg2_successors = set(G.successors(utg2))
            
            if len(utg1_predecessors & utg2_predecessors) == 1 and len(utg1_successors & utg2_successors) == 1:
                graph_allele_dict[tuple(sorted([utg1, utg2]))] = 1
                # print(f"{utg1}\t{utg2}")

    return graph_allele_dict


def nor_partig(weight, ctg1_len, ctg2_len, method="no"):

    if method ==  "no":
        sim = weight
    elif method == "min":
        sim = weight * min([ctg1_len, ctg2_len])
    elif method == "len-penalty":
        L_mean = (ctg1_len + ctg2_len) / 2
        length_ratio = min([ctg1_len, ctg2_len]) / max([ctg1_len, ctg2_len])
        beta = 0.1    
        gamma = 500000 
        penalty = (1 - length_ratio) * beta
        sim = weight * (1-penalty) * L_mean
    else:
        sim = 0
    return sim





def read_partig(partig_file, fai_dict, fai_reverse_dict, output_file, ctg_RE_len, graph_allele_dict, method="no"):
    
    partig_re = defaultdict()
    flag_list = list()
    with open(partig_file, 'r') as file, open(output_file, 'w') as output_f:
        for line in file:
            line = line.strip().split()
            if line[0] == "S":
                scontig1 = line[1][1:]
                scontig2 = line[2][1:]
                contig1 = fai_dict[scontig1]
                contig2 = fai_dict[scontig2]
                if tuple(sorted([contig1, contig2])) in graph_allele_dict:
                    continue
                if contig1 in ctg_RE_len and contig2 in ctg_RE_len and ctg_RE_len[contig1][1] > 50000 and ctg_RE_len[contig2][1] > 50000:
                # if contig1 in ctg_RE_len and contig2 in ctg_RE_len:
                    sim = nor_partig(float(line[7]), ctg_RE_len[contig1][1], ctg_RE_len[contig2][1], method=method)
                else:
                    continue
                _tmp = tuple(sorted([contig1, contig2]))
                if _tmp not in flag_list:
                    flag_list.append(tuple(sorted([contig1, contig2])))
                    output_f.write(f"{contig1},{contig2},{sim}\n")
    
    with open(output_file, 'a') as output_f2:
        for (contig1, contig2) in graph_allele_dict:
            if contig1 in ctg_RE_len and contig2 in ctg_RE_len and ctg_RE_len[contig1][1] > 50000 and ctg_RE_len[contig2][1] > 50000:
            # if contig1 in ctg_RE_len and contig2 in ctg_RE_len:
                sim = nor_partig(1, ctg_RE_len[contig1][1], ctg_RE_len[contig2][1], method=method)
            else:
                continue
            _tmp = tuple(sorted([contig1, contig2]))
            if _tmp not in flag_list:
                flag_list.append(tuple(sorted([contig1, contig2])))
                output_f2.write(f"{contig1},{contig2},{sim}\n")


def main():
    parser = argparse.ArgumentParser(description="trans partig file", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-fai', '--fai', required=True, help='fai for asm.fa')
    parser.add_argument('-r', '--RE_file', required=True, help='length for contig')
    parser.add_argument('-p', '--partig_file', required=True, help='partig file')
    parser.add_argument('-d', '--digraph', required=True, help='digraph file')
    parser.add_argument('-m', '--method', default="no", help='standardized method: no, min, len-penalty')
    parser.add_argument('-o', '--output', required=True, help='output file')
    

    args = parser.parse_args()  
    fai_file = args.fai
    partig_file = args.partig_file
    output_file = args.output
    digraph_file = args.digraph
    RE_file = args.RE_file
    method = args.method

    valid_methods = {"no", "min", "len-penalty"}
    if method not in valid_methods:
        parser.error(f"Error: Invalid method '{method}'. Must be one of: {', '.join(valid_methods)}")

    fai_dict, fai_reverse_dict = read_fai(fai_file)
    ctg_RE_len = read_REs(RE_file)
    graph_allele_dict = add_graph_allele(digraph_file)
    read_partig(partig_file, fai_dict, fai_reverse_dict, output_file, ctg_RE_len, graph_allele_dict, method)

if __name__ == "__main__":
    main()
