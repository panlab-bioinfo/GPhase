from collections import defaultdict
import csv
import networkx as nx
import argparse

def read_c(c):
    cluster_dict = defaultdict(list)
    subgraph_group_dict = defaultdict(list)
    with open(c, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            for utg in line[2].split():
                cluster_dict[line[0]].append(utg)
                subgraph_group_dict[utg].append(line[0])
    return cluster_dict, subgraph_group_dict

def run_trans_allele(hic_groups_cluster_dict, group_hic_dict, allele_ctgs_cluster_dict, ctgs_allele_dict, output_prefix):

    hic_cluster_re_dict = defaultdict(set)

    for hic_group in hic_groups_cluster_dict:
        groups = hic_groups_cluster_dict[hic_group]

        for group in groups:
            hic_cluster_re_dict[hic_group].update(allele_ctgs_cluster_dict[group])

    with open(f"{output_prefix}.chr.cluster.ctg.txt", 'w') as file:
        for group, ctgs in hic_cluster_re_dict.items():
            file.write(f"{group}\t{len(ctgs)}\t{' '.join(list(ctgs))}\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="cluster chr base hic")
    parser.add_argument('-c1', '--clusters1', required=True,
                        help='<filepath>clusters file1')
    parser.add_argument('-c2', '--clusters2', required=True,
                        help='<filepath>clusters file2')
    parser.add_argument('-n', '--name', required=True,
                        help='<str> output prefix')

    args = parser.parse_args()

    hic_groups_cluster = args.clusters1
    allele_ctgs_cluster = args.clusters2
    output_prefix = args.name

    # hic_groups_cluster = "PT4.chr.cluster.txt"
    # allele_ctgs_cluster = "PT4.allele.cluster.ctg.txt"

    # allele ，hic 聚类过后的cluster
    hic_groups_cluster_dict, group_hic_dict = read_c(hic_groups_cluster)

    # allele 的聚类cluster ,value 为ctgs
    allele_ctgs_cluster_dict, ctgs_allele_dict = read_c(allele_ctgs_cluster)

    run_trans_allele(hic_groups_cluster_dict, group_hic_dict, allele_ctgs_cluster_dict, ctgs_allele_dict, output_prefix)





