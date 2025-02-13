from collections import defaultdict
import argparse
import statistics
import copy

def read_collapse_num(collapse_num_file):

    collapse_num_dict = defaultdict()
    with open(collapse_num_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[0].startswith("utg") or line[0].startswith("utig"):
                try:
                    collapse_num_dict[line[0]] = int(line[1])
                except:
                    print(f"error : {line[0]}\t{line[1]}")
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



def run(collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, cluster_dict, utg_group_dict, ctg_RE_len, allele_dict, ctg_allele_dict):

    log_file = open("reassign_collapse.log", 'w')
    collapse_utgs_list = [ utg for utg in collapse_num_dict if collapse_num_dict[utg] > 1 and utg in utgs_list]

    # 计算collapse utg 对每个cluster的hic信号总数和allele的片段总长
    collapse_utg_group_links_dict = defaultdict(lambda: defaultdict(float))
    collapse_utg_group_allele_dict = defaultdict(lambda: defaultdict(float))
    cluster_uncollapse_dict = copy.deepcopy(cluster_dict)

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
        for i in range(collapse_num_dict[collapse_utg]):

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

            hic_list = list(unreassign_groups_hic_sorted.values())
            allele_list = list(unreassign_groups_allele_sorted.values())

            if len(hic_list) == 1:
                cluster_dict[max_hic_group].append(collapse_utg)
            else:
                # 检测最大值是否离群，是否中值的三倍    
                median = statistics.median(hic_list[1:])

                if unreassign_groups_hic[max_hic_group] > median * 2:
                    reassign_list.append(max_hic_group)
                    cluster_dict[max_hic_group].append(collapse_utg)
                    if i==0:
                        cluster_uncollapse_dict[max_hic_group].append(collapse_utg)
                
                else:
                    reassign_list.append(min_allele_group)
                    cluster_dict[min_allele_group].append(collapse_utg)
                    if i==0:
                        cluster_uncollapse_dict[max_hic_group].append(collapse_utg)


    
    with open("reassign_collapse.cluster.txt", 'w') as file:
        for group in cluster_dict:
            file.write(f"{group}\t{len(cluster_dict[group])}\t")
            for utg in cluster_dict[group]:
                file.write(f"{utg} ")
            file.write("\n")

    with open("reassign_collapse.uncollapse.cluster.txt", 'w') as file:
        for group in cluster_uncollapse_dict:
            file.write(f"{group}\t{len(cluster_uncollapse_dict[group])}\t")
            for utg in cluster_uncollapse_dict[group]:
                file.write(f"{utg} ")
            file.write("\n")

            








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
    parser.add_argument('-clusters', '--clusters', required=True,
                        help='<filepath>clusters for utgs')


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

    collapse_num_dict = read_collapse_num(collapse_num_file)
    utgs_list = read_chr_utgs(chr_file)
    hic_links_dict, hic_nei_dict = read_l(l)
    cluster_dict, utg_group_dict = read_c(c)
    ctg_RE_len = read_REs(r)
    allele_dict, ctg_allele_dict = read_allele(a)


    run(collapse_num_dict, utgs_list, hic_links_dict, hic_nei_dict, cluster_dict, utg_group_dict, ctg_RE_len, allele_dict, ctg_allele_dict)


