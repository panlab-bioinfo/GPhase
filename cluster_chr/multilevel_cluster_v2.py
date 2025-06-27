import pandas as pd
import networkit as nk
from collections import defaultdict
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

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len

def Louvain_cluster(g, id_to_node, resolution):

    nk.setSeed(42, True)
    partitioner = nk.community.PLM(g, refine=True, gamma=float(resolution))
    partitioner.run()
    partition = partitioner.getPartition()

    communities = defaultdict(list)
    for node in range(g.numberOfNodes()):
        cluster_id = partition.subsetOf(node)
        label = id_to_node[node]
        communities[cluster_id].append(label)

    return communities

def Multilevel_cluster(csv_file, output_file, resolution, check=None, RE_file=None, Allele_cluster=None, n_chr=None):
    if RE_file:
        ctg_RE_len = read_REs(RE_file)
    if Allele_cluster:
        cluster_dict, subgraph_group_dict = read_c(Allele_cluster)

    df = pd.read_csv(csv_file)
    df['source'] = df['source'].astype(str)
    df['target'] = df['target'].astype(str)
    df['links'] = df['links'].astype(float)

    g = nk.graph.Graph(weighted=True, directed=False)
    nodes = list(set(df['source']).union(set(df['target'])))

    node_to_id, id_to_node = dict(), dict()
    node_to_id = {node: idx for idx, node in enumerate(nodes)}
    id_to_node = {idx: node for node, idx in node_to_id.items()}
    g.addNodes(len(nodes))


    for _, row in df.iterrows():
        u = node_to_id[row['source']]
        v = node_to_id[row['target']]
        weight = row['links']
        if not g.hasEdge(u, v):
            g.addEdge(u, v, weight)

    communities = Louvain_cluster(g, id_to_node, resolution)

    if check:
        filtered_chr_list = list()
        chr_len_dict = defaultdict(int)
        chr_cluster_dict = defaultdict(set)
        for idx, group in communities.items():
            sum_length = 0
            for group_member in group:
                for utg in cluster_dict.get(group_member, []):
                    chr_cluster_dict[idx].add(utg)
                    if utg in ctg_RE_len:
                        sum_length += ctg_RE_len[utg][1]
            chr_len_dict[idx] = sum_length

        # threshold = ctg_len_All / int(n_chr) / 3
        threshold = sum(chr_len_dict.values()) / int(n_chr) / 2
        filtered_chr_list = [key for key, value in chr_len_dict.items() if value > threshold]

        group = 0
        with open(output_file, 'w') as file:
            for idx in chr_cluster_dict:
                if idx in filtered_chr_list:
                    group += 1
                    utgs = list(chr_cluster_dict[idx])
                    file.write(f"group{group}\t{len(utgs)}\t{' '.join(utgs)}\n")

        return len(filtered_chr_list)
    else:
        with open(output_file, 'w') as file:
            for idx, group in communities.items():
                file.write(f"group{idx+1}\t{len(group)}\t{' '.join(group)}\n")
        return len(communities)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Networkit Louvain clustering algorithm based tool.")
    parser.add_argument('-c', '--csv_file', required=True, help='CSV file path containing Hi-C links')
    parser.add_argument('-o', '--output_file', required=True, help='Output file path')
    parser.add_argument('-r', '--resolution', default=1, type=float, help='Resolution parameter for Louvain (gamma)')
    parser.add_argument('--check', action='store_true', help='Enable check mode to filter clusters')
    parser.add_argument('--RE_file', help='File path required when --check is enabled')
    parser.add_argument('--Allele_cluster', help='File path required when --check is enabled')
    parser.add_argument('--n_chr', type=int, metavar='\b', help='Expected number of chromosomes, required when --check is enabled')

    args = parser.parse_args()

    if args.check:
        if not args.RE_file:
            parser.error("--RE_file is required when --check is enabled")
        if not args.Allele_cluster:
            parser.error("--Allele_cluster is required when --check is enabled")
        if args.n_chr is None:
            parser.error("--n_chr is required when --check is enabled")

    if not args.check:
        Multilevel_cluster(args.csv_file, args.output_file, args.resolution)
    else:
        Multilevel_cluster(args.csv_file, args.output_file, args.resolution, args.check, args.RE_file, args.Allele_cluster, int(args.n_chr))
