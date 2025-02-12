from collections import defaultdict, Counter
import csv
import networkx as nx
import argparse

def read_c(c):
    cluster_dict = defaultdict(list)
    utg_group_dict = defaultdict(list)
    with open(c, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            for utg in line[2].split():
                cluster_dict[line[0]].append(str(utg))
                utg_group_dict[utg].append(line[0])
    return cluster_dict, utg_group_dict


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


def subgraph_chr(subgraph_ctgs_dict, ctg_chr_dict):

    subgraph_chr_dict = defaultdict()
    for subgraph_num, ctgs in subgraph_ctgs_dict.items():

        ctg_chr_subgraph_list = [ctg_chr_dict.get(ctg, None) for ctg in ctgs]
        counter = Counter(ctg_chr_subgraph_list)
        most_common_element, most_common_count = counter.most_common(1)[0]
        most_common_two = counter.most_common(2)

        if (most_common_element and most_common_count > len(ctgs) -3) \
            or (len(counter) == 2 and most_common_two == None):
            subgraph_chr_dict[subgraph_num] = most_common_element
        else:
            print(f"ERROR...{subgraph_num} \tlen:{len(ctgs)}\tdifferent!")
            for idx, ctg in enumerate(ctgs, 1):
                print(f"{ctg}\t{idx}\t{ctg_chr_dict.get(ctg, None)}")
    
    return subgraph_chr_dict

def trans_cluster(cluster_dict, subgraph_ctgs_dict, ctg_subgraph_dict, output_prefix):

    ctg_cluster = defaultdict(list)
    ctg_chr_clustr = defaultdict(list)
    for group in cluster_dict:
        for subgraph in cluster_dict[group]:
            ctg_cluster[group].extend(subgraph_ctgs_dict[subgraph])
            # ctg_chr_clustr[group].extend(subgraph_chr_dict[subgraph])

    with open(f"{output_prefix}.allele.cluster.ctg.txt", 'w') as file:
        for group in ctg_cluster:
            file.write(f"{group}\t{len(ctg_cluster[group])}\t{' '.join(ctg_cluster[group])}\n")


            

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="trans cluster to ctgs")
    parser.add_argument('-c', '--clusters', required=True,
                        help='<filepath>clusters file')
    parser.add_argument('-s', '--subgraphs', required=True,
                        help='<filepath>subgraphs file')
    parser.add_argument('-n', '--name', required=True,
                        help='<str> output prefix')

    args = parser.parse_args()

    c = args.clusters
    subgraph_file = args.subgraphs
    output_prefix = args.name



    cluster_dict, utg_group_dict = read_c(c)
    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)
    trans_cluster(cluster_dict, subgraph_ctgs_dict, ctg_subgraph_dict, output_prefix)