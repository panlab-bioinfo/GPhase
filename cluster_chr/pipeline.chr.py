#-----------------------#
# 对allele 聚类的簇再次根据hic进行聚类
#-----------------------#


from collections import defaultdict,Counter
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

# def read_l(l):
#     hic_nei_dict = defaultdict(set)
#     hic_links_dict = defaultdict()
#     with open(l, 'r') as file:
#         for line in file:
#             if not line.startswith("source"):
#                 line = line.strip().split(',')
#                 hic_links_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
#                 hic_nei_dict[line[0]].add(line[1])
#                 hic_nei_dict[line[1]].add(line[0])
#     return hic_links_dict, hic_nei_dict

def read_l(l):
    hic_nei_dict = defaultdict(set)
    hic_links_dict = {}

    with open(l, newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            if not row:
                continue
            if row[0].startswith(('utg', 'utig')):
                a, b, w = row[0], row[1], float(row[2])

                # 用 tuple 排序缓存减少运算
                key = (a, b) if a < b else (b, a)

                hic_links_dict[key] = w
                hic_nei_dict[a].add(b)
                hic_nei_dict[b].add(a)

    return hic_links_dict, hic_nei_dict


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

def cluster(cluster_dict, subgraph_group_dict, subgraph_ctgs_dict, ctg_subgraph_dict, hic_links_dict, output_prefix):


    allele_hic_dict = defaultdict(float)
    for (ctg1, ctg2), value  in hic_links_dict.items():

        subgraph1 =  ctg_subgraph_dict.get(ctg1, None)
        subgraph2 =  ctg_subgraph_dict.get(ctg2, None)
        
        allele_group_1 = subgraph_group_dict.get(subgraph1, [])
        allele_group_2 = subgraph_group_dict.get(subgraph2, [])

        print(f"{ctg1}\t{ctg2}\t{allele_group_1}\t{allele_group_2}")

        if allele_group_1 and allele_group_2:
            allele_hic_dict[tuple(sorted([allele_group_1[0], allele_group_2[0]]))] += value

    with open(f'{output_prefix}.allele.hic.csv', 'w') as file:
        file.write("source,target,links\n")

        for pair, value in allele_hic_dict.items():
            group1 = pair[0]
            group2 = pair[1]
            if  group1 != group2:
                file.write(f"{group1},{group2},{value}\n")
    

    return allele_hic_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="cluster chr base hic")
    parser.add_argument('-c', '--clusters', required=True,
                        help='<filepath>clusters file')
    parser.add_argument('-l', '--links', required=True,
                        help='<filepath>hic links for ctgs')
    parser.add_argument('-s', '--subgraphs', required=True,
                        help='<filepath>subgraphs file')
    parser.add_argument('-n', '--name', required=True,
                        help='<str> output prefix')

    args = parser.parse_args()

    c = args.clusters
    l = args.links
    subgraph_file = args.subgraphs
    output_prefix = args.name


    cluster_dict, subgraph_group_dict = read_c(c)
    hic_links_dict, hic_nei_dict = read_l(l)
    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)
    cluster(cluster_dict, subgraph_group_dict, subgraph_ctgs_dict, ctg_subgraph_dict, hic_links_dict, output_prefix)




    
    






