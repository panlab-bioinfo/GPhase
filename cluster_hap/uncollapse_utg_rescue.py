from collections import defaultdict
import copy
import argparse


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


def run_uncollapse_rescue(c, l, collapse_num_file):
    cluster_dict, utg_group_dict = read_c(c)
    hic_links_dict, hic_nei_dict = read_l(l)
    collapse_num_dict = read_collapse_num(collapse_num_file)

    uncollapse_all_list = [ utg for utg in collapse_num_dict if collapse_num_dict[utg] < 2 ]
    rescue_uncollapse_list = set(uncollapse_all_list) - set(utg_group_dict.keys())

    # print(f"uncollapse_all_list: {rescue_uncollapse_list}")
    re_cluster_dict = copy.deepcopy(cluster_dict)
    rescue_utg_links_dict = defaultdict(lambda: defaultdict(float))
    for utg in rescue_uncollapse_list:
        group_sum = 0
        for group in cluster_dict:
            for utg2 in cluster_dict[group]:
                if tuple(sorted([utg, utg2])) in hic_links_dict:
                    group_sum += hic_links_dict[tuple(sorted([utg, utg2]))]
            rescue_utg_links_dict[utg][group] = group_sum

        sorted_dict = dict(sorted(rescue_utg_links_dict[utg].items(), key=lambda item: item[1], reverse=True))
        # print(f"{utg}\t{sorted_dict}")
        
        if sorted_dict[list(sorted_dict)[0]] > 0:
            re_cluster_dict[list(sorted_dict)[0]].append(utg)
    
    with open("uncollapse_rescue.cluster.txt", 'w') as file:
        for group in re_cluster_dict:
            file.write(f"{group}\t{len(re_cluster_dict[group])}\t")
            for utg in re_cluster_dict[group]:
                file.write(f"{utg} ")
            file.write("\n")
            

if __name__ == "__main__":
    # c = "test.cluster.txt"
    # l = "c88.chr10.links.nor.csv"
    # collapse_num_file = "06.genes.round.cn"
    # run_uncollapse_rescue(c, l, collapse_num_file)

    parser = argparse.ArgumentParser(description="rescue uncollapse utg")
    parser.add_argument('-c', '--cluster', required=True,
                        help='<filepath> clusters file')
    parser.add_argument('-l', '--links', required=True,
                        help='<filepath> hic links file')
    parser.add_argument('-n', '--collapse_num', required=True,
                        help='<filepath> collapse num for utgs')

    args = parser.parse_args()

    run_uncollapse_rescue(args.cluster, args.links, args.collapse_num)
