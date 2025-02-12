from collections import defaultdict
import copy as cp
import argparse
import os
import re

def parse_size(size_str):
    """
    解析带单位的字符串，如 500k, 500M, 1G 等，并将其转换为整数。
    支持的单位有 'k' (千), 'm' (百万), 'g' (十亿) 等，不区分大小写。
    """
    size_str = size_str.strip().lower()  # 转换为小写，确保大小写不敏感
    
    # 正则匹配：数字 + 单位
    match = re.match(r"(\d+)([kmg]?)", size_str)
    
    if match:
        num = int(match.group(1))
        unit = match.group(2)
        
        if unit == "k":
            return num * 1000
        elif unit == "m":
            return num * 1000000
        elif unit == "g":
            return num * 1000000000
        else:
            return num  # 如果没有单位，直接返回数字
    
    raise ValueError(f"Invalid size string: {size_str}")

def parse_agp(file_path):
    """解析 AGP 文件为 scaffolds 数据结构。"""

    scaffolds = defaultdict(list)
    scaffolds_length_dict = defaultdict(int)
    with open(file_path, 'r') as infile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue
            columns = line.strip().split('\t')
            scaffold_name = columns[0]
            start = int(columns[1])
            end = int(columns[2])
            scaffolds_length_dict[scaffold_name] += int(end - start)
            idx = int(columns[3])
            object_type = columns[4]
            if object_type == 'W':  
                contig_id = columns[5]
                strand = columns[8]
                scaffolds[scaffold_name].append((start, end, idx ,strand, contig_id, columns))

    return scaffolds, scaffolds_length_dict


def trans_scaffolds(scaffolds, scaffolds_length_dict, len_threshold:"1m", output_file):

    len_threshold_trans =  parse_size("1m")

    if len_threshold != "1m":
        try:
            len_threshold_trans = parse_size(len_threshold)
        except:
            print(f"ERROR... Scaffold length threshold {len_threshold} parsing error.")

    scaffolds_length_sorted_dict = dict(sorted(scaffolds_length_dict.items(), key=lambda item: item[1], reverse=True))




    scaffold_offset, idx_offset = 0,0
    not_gap = False
    for idx, scaffold_name in enumerate(list(scaffolds_length_sorted_dict)):

        if idx == len(list(scaffolds_length_sorted_dict))-1:
            not_gap = True
        else:
            next_scaffold_name = list(scaffolds_length_sorted_dict)[idx+1] 
            if scaffolds_length_dict[next_scaffold_name] < len_threshold_trans:
                not_gap = True

        if scaffolds_length_dict[scaffold_name] >= len_threshold_trans:

            with open(output_file, 'a') as file:
                for idx, ( _, _, _,_, _, columns) in enumerate(scaffolds[scaffold_name]):
                    new_columns = cp.deepcopy(columns)

                    new_columns[0] = "scaffold_1"
                    new_columns[1] = str(int(columns[1]) + scaffold_offset)
                    new_columns[2] = str(int(columns[2]) + scaffold_offset)
                    new_columns[3] = str(int(columns[3]) + idx_offset)

                    file.write('\t'.join(new_columns)+'\n')

                    if idx == len(scaffolds[scaffold_name]) -1:
                        scaffold_offset += int(columns[2]) + 100
                        idx_offset += int(columns[3]) + 1

                    if not not_gap or idx != len(scaffolds[scaffold_name])-1:
                        file.write(f"scaffold_{'1'}\t{int(new_columns[2])+1}\t{int(new_columns[2])+100}\t{int(new_columns[3])+1}\t{'U'}\t{100}\t{'scaffold'}\t{'yes'}\t{'proximity_ligation'}\n")


def main():
    parser = argparse.ArgumentParser(description="The scaffold is merged according to the relationship between the AGPs")
    parser.add_argument('-agp', '--agp_file', help='target AGP file', required=True)
    parser.add_argument('-cut_off', '--cut_off', help='length threshold')
    parser.add_argument('-o', '--output_file', help='Output file name', required=True)

    args = parser.parse_args()
    agp_file = args.agp_file 
    cut_off = args.cut_off 
    output_file = args.output_file

                    

    # agp1_file = "yahs_rescaffold_scaffolds_final.agp"  
    # agp2_file = "yahs_scaffold_scaffolds_final.agp" 
    # output_file = "merged.agp" 

    agp_scaffolds, agp_scaffolds_length_dict = parse_agp(agp_file)

    if cut_off:
        trans_scaffolds(agp_scaffolds, agp_scaffolds_length_dict, cut_off, output_file)
    else:
        trans_scaffolds(agp_scaffolds, agp_scaffolds_length_dict, "1m",output_file)


if __name__=="__main__":
    main()