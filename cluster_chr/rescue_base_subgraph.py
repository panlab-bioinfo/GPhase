from collections import defaultdict,Counter
import copy as cp
import csv
import networkx as nx
import copy
import argcomplete
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



def get_chr_cluster(subgraph_unCluster_ctgs_dict, ctg_subgraph_unCluster_dict, rmSubgraph_ctgs_dict, ctg_rmSubgraph_dict,hicGroup_cluster_dict, ctg_hicGroups_dict, hic_links_dict):

    subgraph_unClusterRm_hicGroup_dict = defaultdict(lambda: defaultdict(int))
    hicGroup_cluster_re_dict = copy.deepcopy(hicGroup_cluster_dict)

    # merge subgraph
    merge_ctg_subgraph_dict = {**ctg_subgraph_unCluster_dict, **ctg_rmSubgraph_dict}
    merge_subgraph_ctg_dict = {**subgraph_unCluster_ctgs_dict, **rmSubgraph_ctgs_dict}


    for (ctg1,ctg2), value in hic_links_dict.items():

        subgraph1 = merge_ctg_subgraph_dict.get(ctg1, None)
        subgraph2 = merge_ctg_subgraph_dict.get(ctg2, None)

        if (not subgraph1 and not subgraph2) or (subgraph1 and subgraph2):
            continue
        else:
            if not subgraph1:

                subgraph1_list = ctg_hicGroups_dict.get(ctg1, None)
                if not subgraph1_list:
                    # print(f"error... {ctg2} in subgraph{subgraph2}\t{ctg1} not in hic clusters.")
                    pass
                else:
                    subgraph1 = subgraph1_list[0]

                subgraph_unClusterRm_hicGroup_dict[subgraph2][subgraph1] += value

            elif not subgraph2:

                subgraph2_list = ctg_hicGroups_dict.get(ctg2, None)
                if not subgraph2_list:
                    # print(f"error... {ctg1} in subgraph{subgraph1}\t{ctg2} not in hic clusters.")
                    pass
                else:
                    subgraph2 = subgraph2_list[0]

                subgraph_unClusterRm_hicGroup_dict[subgraph1][subgraph2] += value
            else:
                print("error...")

    for group, dict_ in subgraph_unClusterRm_hicGroup_dict.items():
        sorted_dict_ = dict(sorted(dict_.items(), key=lambda item: item[1], reverse=True))
        if group == 166 or group == "166":
            print(sorted_dict_)

        if not sorted_dict_:
            continue

        target_group = list(sorted_dict_)[0]

        ctgs = merge_subgraph_ctg_dict[group]
        hicGroup_cluster_re_dict[target_group].extend(ctgs)
    
    
    with open("rescue.cluster.ctg.txt", "w") as file:

        for group, ctgs in hicGroup_cluster_re_dict.items():
            file.write(f"{group}\t{len(ctgs)}\t{' '.join(ctgs)}\n")




def run_rescue(cluster_dict, subgraph_group_dict, subgraph_ctgs_dict, ctg_subgraph_dict, hicGroup_cluster_dict, ctg_hicGroups_dict):

    subgraph_all_list = list(subgraph_ctgs_dict)

    # check 子图或ctg是否因hic link过滤丢失
    hicGroup_cluster_dict_cp = cp.deepcopy(hicGroup_cluster_dict)

    for group in hicGroup_cluster_dict_cp:
        for ctg in hicGroup_cluster_dict_cp[group]:
            subgraph = ctg_subgraph_dict[ctg]
            for same_subgraph_ctg in subgraph_ctgs_dict[subgraph]:
                if same_subgraph_ctg not in ctg_hicGroups_dict:
                    print(same_subgraph_ctg)
                    hicGroup_cluster_dict[group].append(same_subgraph_ctg)
    
    uncluster_ctg_list = [ ctg for ctg in ctg_subgraph_dict if ctg not in ctg_hicGroups_dict ]
    subgraph_unCluster_list = [ ctg_subgraph_dict[ctg] for ctg in uncluster_ctg_list ]

    # subgraph_unCluster_list = set(subgraph_all_list) - set(subgraph_in_clusuter_list)


    subgraph_unCluster_ctgs_dict = { subgraph_num:ctgs for subgraph_num,ctgs in subgraph_ctgs_dict.items() \
                                                        if subgraph_num in subgraph_unCluster_list}

    ctg_subgraph_unCluster_dict = { ctg:subgraph_num for subgraph_num,ctgs in subgraph_unCluster_ctgs_dict.items() \
                                                         for ctg in ctgs}

    return subgraph_unCluster_ctgs_dict, ctg_subgraph_unCluster_dict

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="rescue base subgraph")
    parser.add_argument("--allele_cluster", required=True, help="allele cluster file")
    parser.add_argument("--chr_cluster_ctg", required=True, help="chr cluster ctgs file")
    parser.add_argument("-s", "--subgraph_file", required=True, help="Subgraph file")
    parser.add_argument("-rs", "--rm_subgraph_file", required=True, help="rm subgraph")
    parser.add_argument("-l", "--hic_link", required=True, help="hic links")


    # c = "wax.allele.cluster.epxand.txt"
    # subgraph_file = "group_ctgs_All.txt"
    # hicGroup_clusters_file = "wax.chr.cluster.ctg.txt"
    # l = "wax.links.nor.csv"
    # rm_subgraph_file = "group_ctgs_filter.txt"

    args = parser.parse_args()
    argcomplete.autocomplete(parser)

    c = args.allele_cluster
    subgraph_file = args.subgraph_file
    hicGroup_clusters_file = args.chr_cluster_ctg
    l = args.hic_link
    rm_subgraph_file = args.rm_subgraph_file



    hic_links_dict, hic_nei_dict = read_l(l)

    cluster_dict, subgraph_group_dict = read_c(c)
    hicGroup_cluster_dict, ctg_hicGroups_dict = read_c(hicGroup_clusters_file)

    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)
    rmSubgraph_ctgs_dict, ctg_rmSubgraph_dict = read_subgraph(rm_subgraph_file)

    subgraph_unCluster_ctgs_dict, ctg_subgraph_unCluster_dict = run_rescue(cluster_dict, subgraph_group_dict, subgraph_ctgs_dict, ctg_subgraph_dict, hicGroup_cluster_dict,ctg_hicGroups_dict)

    get_chr_cluster(subgraph_unCluster_ctgs_dict, ctg_subgraph_unCluster_dict, rmSubgraph_ctgs_dict, ctg_rmSubgraph_dict,hicGroup_cluster_dict, ctg_hicGroups_dict, hic_links_dict)


    
