from collections import defaultdict
import pandas as pd
import networkx as nx
import argparse

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

def read_partig(partig_file, fai_dict, fai_reverse_dict, output_file, graph_allele_dict):
    
    partig_re = defaultdict()

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

                output_f.write(f"{contig1},{contig2},{line[7]}\n")

    with open(output_file, 'a') as output_f2:
        for (utg1, utg2) in graph_allele_dict:
            sim = 1 
            output_f2.write(f"{utg1},{utg2},{sim}\n")


def main():
    parser = argparse.ArgumentParser(description="trans partig file", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-fai', '--fai', required=True, help='fai for asm.fa')
    parser.add_argument('-p', '--partig_file', required=True, help='partig file')
    parser.add_argument('-d', '--digraph', required=True, help='digraph file')
    parser.add_argument('-o', '--output', required=True, help='output file')

    args = parser.parse_args()  
    fai_file = args.fai
    partig_file = args.partig_file
    output_file = args.output
    digraph_file = args.digraph

    fai_dict, fai_reverse_dict = read_fai(fai_file)
    graph_allele_dict = add_graph_allele(digraph_file)
    read_partig(partig_file, fai_dict, fai_reverse_dict, output_file, graph_allele_dict)

if __name__ == "__main__":
    main()
