from collections import defaultdict
import argparse
import statistics
import copy
import numpy as np

def read_collapse_num(collapse_num_file):

    collapse_num_dict = defaultdict()
    with open(collapse_num_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0].startswith("utg") or line[0].startswith("utig"):
                try:
                    collapse_num_dict[line[0]] = int(line[1])
                except:
                    collapse_num_dict[line[0]] = 1
    return collapse_num_dict



def read_chr_utgs(chr_file):
    utgs_list = list()
    with open(chr_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0].startswith("utg") or line[0].startswith("utig"):
                utgs_list.append(line[0])
    return utgs_list

def read_c(c):
    cluster_dict = defaultdict(list)
    utg_group_dict = defaultdict(list)
    with open(c, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            for utg in line[2].split():
                cluster_dict[line[0]].append(utg)
                utg_group_dict[utg].append(line[0])
    return cluster_dict, utg_group_dict

def correct_collapse_num(utgs_list, collapse_num_dict, cluster_dict, utg_group_dict):

    correct_collapse_num_dict = copy.deepcopy(collapse_num_dict)

    for utg in utgs_list:
        if utg not in utg_group_dict and collapse_num_dict[utg] == 1:
            correct_collapse_num_dict[utg] = -1

    for group in cluster_dict:
        for ctg in cluster_dict[group]:
            if ctg in correct_collapse_num_dict:
                correct_collapse_num_dict[ctg] -= 1
    return correct_collapse_num_dict

def read_l(l):
    hic_nei_dict = defaultdict(set)
    hic_links_dict = defaultdict()
    with open(l, 'r') as file:
        for line in file:
            if line.startswith("utg") or line.startswith("utig"):
                line = line.strip().split(',')
                hic_links_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
                hic_nei_dict[line[0]].add(line[1])
                hic_nei_dict[line[1]].add(line[0])
    return hic_links_dict, hic_nei_dict


def read_REs(REFile):

    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))

    return ctg_RE_len

def read_allele(allele_file):

    allele_dict = defaultdict()
    ctg_allele_dict = defaultdict(set)
    with open(allele_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if line[0].startswith("u") and line[1].startswith("u"):

                allele_dict[tuple(sorted([line[0], line[1]]))] = float(line[2]) if len(line) > 2 else 1
                ctg_allele_dict[line[0]].add(line[1])
                ctg_allele_dict[line[1]].add(line[0])

    return allele_dict, ctg_allele_dict

def allele_sort(unreassign_groups_hic):
    def sort_key(item):
        key, value = item
        if value == 0:
            # 如果值为0，使用dict2中的值进行逆序排序
            return (-unreassign_groups_hic.get(key, float('inf')), key)
        else:
            # 如果值不为0，按照升序排序
            return (value, key)
    return sort_key



def run(correct_collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, cluster_dict, utg_group_dict, ctg_RE_len, allele_dict, ctg_allele_dict, isolated_threshold, output_prefix):

    log_file = open("reassign_collapse.log", 'a')
    collapse_utgs_list = [ utg for utg in correct_collapse_num_dict if (correct_collapse_num_dict[utg] > 1 or correct_collapse_num_dict[utg] == -1)and utg in utgs_list]

    # 计算collapse utg 对每个cluster的hic信号总数和allele的片段总长
    collapse_utg_group_links_dict = defaultdict(lambda: defaultdict(float))
    collapse_utg_group_allele_dict = defaultdict(lambda: defaultdict(float))
    cluster_uncollapse_dict = copy.deepcopy(cluster_dict)
    n_hap = len(cluster_dict)

    for collapse_utg in collapse_utgs_list:
        for group in cluster_dict:
            group_links, allele_links = 0, 0
            for utg in cluster_dict[group]:
                if utg != collapse_utg:

                    tmp = tuple(sorted([utg, collapse_utg]))
                    group_links += float(hic_links_dict[tmp]) if tmp in hic_links_dict else 0

                    allele_links += ctg_RE_len[utg][1] \
                                    if utg in ctg_allele_dict[collapse_utg] and utg in ctg_RE_len \
                                    else 0



            collapse_utg_group_links_dict[collapse_utg][group] = float(group_links)
            collapse_utg_group_allele_dict[collapse_utg][group] = float(allele_links)


        log_file.write(f"{collapse_utg}\thic links: {collapse_utg_group_links_dict[collapse_utg]}\n")
        log_file.write(f"{collapse_utg}\talleles: {collapse_utg_group_allele_dict[collapse_utg]}\n")



    # 检测hic信号中最大值是否属于离群值，若不为离群值则使用按照allele选择；若为离群值，则按照hic分配此值
    for collapse_utg in collapse_utg_group_links_dict:


        reassign_list = list()
        for i in range(abs(correct_collapse_num_dict[collapse_utg])):

            if i >= n_hap:
                break

            # 未分配的group
            unreassign_groups_hic = { group:value for group, value in collapse_utg_group_links_dict[collapse_utg].items() \
                                            if group not in reassign_list }
            unreassign_groups_allele = { group:value for group, value in collapse_utg_group_allele_dict[collapse_utg].items() \
                                            if group not in reassign_list }

            unreassign_groups_hic_sorted = dict(sorted(unreassign_groups_hic.items(), \
                                            key=lambda item: item[1], reverse=True))
            # 在存在多个0的情况下，则按照hic大小排序
            unreassign_groups_allele_sorted = dict(sorted(unreassign_groups_allele.items(), \
                                            key=allele_sort(unreassign_groups_hic)))

            # print(f"{collapse_utg}\t{unreassign_groups_hic_sorted}\n{unreassign_groups_allele}\t{unreassign_groups_allele_sorted}")

            
            max_hic_group = list(unreassign_groups_hic_sorted)[0]
            min_allele_group = list(unreassign_groups_allele_sorted)[0]
            min_allele_group, max_allele_group = list(unreassign_groups_allele_sorted)[0], list(unreassign_groups_allele_sorted)[n_hap - len(reassign_list) -1]

            if unreassign_groups_hic_sorted[max_hic_group] == 0 and unreassign_groups_allele_sorted[max_allele_group] ==0:
                continue

            # 记录 可以分配的 group
            hic_list = list(unreassign_groups_hic_sorted.keys())
            allele_list = list(unreassign_groups_allele_sorted.values())

            if len(hic_list) == 1:
                cluster_dict[max_hic_group].append(collapse_utg)
            else:  
                list_ = list(unreassign_groups_hic_sorted.values())
                if len(list_) == 2:
                    list_.remove(float(list(unreassign_groups_hic_sorted.values())[0]))
                mean = statistics.mean(list_)
                # print(f"{collapse_utg}\t{max_hic_group}\t{mean*isolated_threshold}")

                if unreassign_groups_hic[max_hic_group] > mean*isolated_threshold:
                    reassign_list.append(max_hic_group)
                    cluster_dict[max_hic_group].append(collapse_utg)
                    if i==0:
                        cluster_uncollapse_dict[max_hic_group].append(collapse_utg)
                
                # 不显著离群，使用 allele 信息 
                else:

                    reassign_list.append(min_allele_group)
                    cluster_dict[min_allele_group].append(collapse_utg)
                    if i==0:
                        cluster_uncollapse_dict[max_hic_group].append(collapse_utg)


    
    with open(f"{output_prefix}.reassign.cluster.txt", 'w') as file:
        for group in cluster_dict:
            file.write(f"{group}\t{len(cluster_dict[group])}\t")
            for utg in cluster_dict[group]:
                file.write(f"{utg} ")
            file.write("\n")

    with open(f"{output_prefix}.reassign.uncopy.cluster.txt", 'w') as file:
        for group in cluster_uncollapse_dict:
            file.write(f"{group}\t{len(cluster_uncollapse_dict[group])}\t")
            for utg in cluster_uncollapse_dict[group]:
                file.write(f"{utg} ")
            file.write("\n")


    # 计算 cluster 中每个簇的长度方差
    group_len_sum_list = list()
    for group, ctgs in cluster_dict.items():
        group_len_sum = sum([ ctg_RE_len[ctg][1] for ctg in ctgs if ctg in ctg_RE_len ])
        group_len_sum_list.append(group_len_sum)
    variance_stats  = statistics.variance(group_len_sum_list)

    return float(variance_stats)




def louvain_reassign_allele(collapse_num_file, chr_file, l, c ,r, a, output_prefix, find_best_isolated,isolated_threshold=5):

    collapse_num_dict = read_collapse_num(collapse_num_file)
    utgs_list = read_chr_utgs(chr_file)
    hic_links_dict, hic_nei_dict = read_l(l)
    cluster_dict, utg_group_dict = read_c(c)
    ctg_RE_len = read_REs(r)
    allele_dict, ctg_allele_dict = read_allele(a)

    correct_collapse_num_dict = correct_collapse_num(utgs_list, collapse_num_dict, cluster_dict, utg_group_dict)
    if find_best_isolated:
        variance_list = list()
        isolated_list = np.arange(0, 10.5, 0.5).tolist()
        for isolated in isolated_list:
            cluster_dict, utg_group_dict = read_c(c)
            variance = run(correct_collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, cluster_dict, utg_group_dict, ctg_RE_len, allele_dict, ctg_allele_dict, isolated, output_prefix)
            variance_list.append(variance)
        
        min_variance_idx = variance_list.index(min(variance_list))
        cluster_dict, utg_group_dict = read_c(c)
        run(correct_collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, cluster_dict, utg_group_dict, ctg_RE_len, allele_dict, ctg_allele_dict, isolated_list[min_variance_idx], output_prefix)
        return isolated_list[min_variance_idx]
    else:
        cluster_dict, utg_group_dict = read_c(c)
        run(correct_collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, cluster_dict, utg_group_dict, ctg_RE_len, allele_dict, ctg_allele_dict, isolated_threshold, output_prefix)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="reassign collapse utgs for clusters")
    parser.add_argument('-c', '--collase_num', required=True,
                        help='<filepath> collapse num for utgs')
    parser.add_argument('-chr', '--chromosome', required=True,
                        help='<filepath>all utg for chrompsome')
    parser.add_argument('-l', '--links', required=True,
                        help='<filepath>hic links for utgs')
    parser.add_argument('-r', '--REs', required=True,
                        help='<filepath>REs and Length for utgs')
    parser.add_argument('-a', '--allele', required=True,
                        help='<filepath>alleles links for utgs')
    parser.add_argument('--clusters', required=True,
                        help='<filepath>clusters for utgs')
    parser.add_argument('-op','--output_prefix',required=True,
                        help='output file prefix')
    parser.add_argument('--find_best_isolated', action='store_true', help='Iterate to find the best isolated_threshold')
    parser.add_argument('--isolated_threshold',type=float,default=5,
                        help='<int>Detect whether the intensity of the hic signal is an outlier(required when --iter is disabled)')


    # argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # collapse_num_file = "06.genes.round.cn"
    # chr_file = "chr06.txt"
    # l = "rice4.links.nor.csv"

    collapse_num_file = args.collase_num
    chr_file = args.chromosome
    l = args.links
    c = args.clusters
    r = args.REs
    a = args.allele
    output_prefix = args.output_prefix
    isolated_threshold = args.isolated_threshold

    # if args.find_best_isolated and args.isolated_threshold:
    #     parser.error("--isolated_threshold is disable when --find_best_isolated is enabled")

    louvain_reassign_allele(collapse_num_file, chr_file, l, c ,r, a, output_prefix, args.find_best_isolated,isolated_threshold)


