from collections import defaultdict, deque
import csv
import networkx as nx
from transDiGraph import run_trans_digraph
import pandas as pd
import os
import argparse


def find_n_nei(digraph, utg, n, dir):
    neighbors = set()
    current_level = {utg}
    
    for _ in range(n):
        next_level = set()
        for current_node in current_level:
            if str(dir) == "predecessors":
                next_level.update(digraph.predecessors(current_node))
            elif str(dir) == "successors":
                next_level.update(digraph.successors(current_node))
        neighbors.update(next_level)
        current_level = next_level

    neighbors.add(utg)
    return neighbors
    
def check_end_point(digraph, utg):

    predecessors = set(digraph.predecessors(utg))
    successors = set(digraph.successors(utg)) 

    if (not predecessors or not successors) and len(predecessors | successors) == 1:
            return True

    return False




def split_GFA(gfa_filePath, digraph, nei_level:2, output_prefix):

    # unsafety_set = set(utg for utg in digraph.nodes() if utg not in safety_set)
    # print(unsafety_set)
    predecessors_rmEdges_set = set()
    successors_rmEdges_set = set()
    split_utg_set = set()

    for utg in digraph.nodes():


        predecessors = set(digraph.predecessors(utg))
        successors = set(digraph.successors(utg)) 


        if len(successors) > 1:


            successors_nei_dict = {utg:find_n_nei(digraph, utg, nei_level, "successors") for utg in successors}
            successors_re_dict = { tuple(sorted([utg1, utg2])):set(set1 & set2)
                                    for utg1, set1 in successors_nei_dict.items() 
                                    for utg2, set2 in successors_nei_dict.items()
                                    if utg1 != utg2}
            
            has_empty_set = any(len(value) == 0 for value in successors_re_dict.values())
            if has_empty_set:
                successors_rmEdges_set.update([(utg, successor) for successor in successors])
                split_utg_set.update([successor for successor in successors])
                split_utg_set.add(utg)

        if len(predecessors) > 1:

            predecessors_nei_dict = {utg:find_n_nei(digraph, utg, nei_level, "predecessors") for utg in predecessors}
            
            predecessors_re_dict = { tuple(sorted([utg1, utg2])):set(set1 & set2)
                                    for utg1, set1 in predecessors_nei_dict.items() 
                                    for utg2, set2 in predecessors_nei_dict.items()
                                    if utg1 != utg2}

            has_empty_set = any(len(value) == 0 for value in predecessors_re_dict.values())
            if has_empty_set:
                predecessors_rmEdges_set.update([(predecessor, utg) for predecessor in predecessors])
                split_utg_set.update([predecessor for predecessor in predecessors])
                split_utg_set.add(utg)


    all_rmEdges_set = successors_rmEdges_set | predecessors_rmEdges_set

    with open(gfa_filePath, 'r') as file, open(f"{output_prefix}.rmTip.split.gfa",  'w') as output_file:
        for line in file:
            line = line.strip().split()
            if(line[0] == "S"):
                output_file.write('\t'.join(line) + '\n')
            elif(line[0] == "L"):
                if (line[1], line[3]) in all_rmEdges_set or (line[3], line[1]) in all_rmEdges_set:
                    continue
                else:
                    output_file.write('\t'.join(line) + '\n')             




if __name__ == '__main__':

    gfa_filePath = "tetra.asm.bp.p_utg.noseq.gfa"


    parser = argparse.ArgumentParser(description="split gfa", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gfa_file', required=True, help='gfa file')
    parser.add_argument('-n', '--n', required=False,default=2, type=int, help='Nth-order neighbor')
    parser.add_argument('-iter', '--iter_number', required=False, default=3,type=int, help='Number of iterations')
    parser.add_argument('-o', '--output_prefix', required=True, help='output prefix')


    args = parser.parse_args() 

    gfa_file = args.gfa_file
    n = args.n
    iter_number = args.iter_number
    output_prefix = args.output_prefix


    digraph = run_trans_digraph(gfa_file, output_prefix, flag_output_csv=True)
    split_GFA(gfa_file, digraph, n, output_prefix)

    for i in range(iter_number):
        old_name = f"{output_prefix}.rmTip.split.gfa"
        new_name = f"{output_prefix}.rmTip.split.{i+1}.gfa"
        os.rename(old_name, new_name)

        digraph = run_trans_digraph(new_name, output_prefix)
        split_GFA(new_name, digraph, n, output_prefix)









