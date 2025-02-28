from collections import defaultdict
import pickle
import argcomplete
from argcomplete.completers import FilesCompleter
import argparse
import pysam

def read_agp(agp_file):

    scaffold_ctgs_dict = defaultdict(lambda: defaultdict(tuple))
    ctg_scaffold_dict = defaultdict(list)
    scaffold_length_dict = defaultdict(int)
    with open(agp_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[5][0]=="u":

                # scaffold_idx_1, scaffold_idx_2, contig_idx_1, contig_idx_1
                scaffold_ctgs_dict[line[0]][line[5]] = (int(line[1]), int(line[2]), int(line[6]), int(line[7]))
                ctg_scaffold_dict[line[5]].append((int(line[6]), int(line[7]), line[0]))

    for scaffold, ctgs in scaffold_ctgs_dict.items():
        scaffold_length_dict[scaffold] = ctgs[list(ctgs.keys())[-1]][1]

    return scaffold_ctgs_dict, ctg_scaffold_dict, scaffold_length_dict

def get_RE(RE_file, scaffold_ctgs_dict, scaffold_length_dict, output_prefix):

    ctg_RE_dict = defaultdict(tuple)
    scaffold_RE_dict = defaultdict(int)
    with open(RE_file, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (int(line[1]), int(line[2]))

    for scaffold in scaffold_ctgs_dict:
        for contig, pair in scaffold_ctgs_dict[scaffold].items():
            if contig in ctg_RE_dict:
                scaffold_RE_dict[scaffold] += ((pair[3]-pair[2]+1) / ctg_RE_dict[contig][1]) * ctg_RE_dict[contig][0]
    
    with open(f"{output_prefix}.scaffold.txt", 'w') as file:
        file.write("#Contig\tRECounts\tLength\n")
        for scaffold, value in scaffold_RE_dict.items():
            file.write(f"{scaffold}\t{int(value)}\t{scaffold_length_dict[scaffold]}\n")


    return scaffold_RE_dict


def read_map_file(map_file, map_file_type, scaffold_ctgs_dict, ctg_scaffold_dict, scaffold_length_dict, output_prefix):

    scaffold_HT_dict = defaultdict(int)
    scaffold_clm_dict = defaultdict(lambda: defaultdict(list))

    if map_file_type == "bam":
        map_file_iter = pysam.AlignmentFile(map_file, "rb")

        for line in map_file_iter:

            ctg1, ctg1_map_idx = line.reference_name, line.reference_start
            ctg2, ctg2_map_idx = line.next_reference_name, line.next_reference_start

            if ctg1 == ctg2 or ctg2 == "=":
                continue
            
            scaffold_1_None_flag , scaffold_2_None_flag= True, True
            for (idx_1, idx_2, scaffold) in ctg_scaffold_dict[ctg1]:
                if int(ctg1_map_idx) > idx_1 and int(ctg1_map_idx) < idx_2:
                    scaffold_1 = scaffold
                    scaffold_1_None_flag = False
                    break
            for (idx_1, idx_2, scaffold) in ctg_scaffold_dict[ctg2]:
                if int(ctg2_map_idx) > idx_1 and int(ctg2_map_idx) < idx_2:
                    scaffold_2 = scaffold
                    scaffold_2_None_flag = False
                    break

            if scaffold_1_None_flag or scaffold_2_None_flag or scaffold_1 == scaffold_2:
                continue

            ctg1_in_scaffold_idx = scaffold_ctgs_dict[scaffold_1][ctg1][0] + int(ctg1_map_idx)
            ctg2_in_scaffold_idx = scaffold_ctgs_dict[scaffold_2][ctg2][0] + int(ctg2_map_idx)

            # HT
            if ctg1_in_scaffold_idx * 2 > scaffold_length_dict[scaffold_1]:
                scaffold_1_pos = "T"
            else:
                scaffold_1_pos = "H"
            
            if ctg2_in_scaffold_idx * 2 > scaffold_length_dict[scaffold_2]:
                scaffold_2_pos = "T"
            else:
                scaffold_2_pos = "H"

            scaffold_HT_dict[tuple(sorted([scaffold_1+"_"+scaffold_1_pos, scaffold_2+"_"+scaffold_2_pos]))] += 1

            # clm
            dir_0 = scaffold_length_dict[scaffold_1] - ctg1_in_scaffold_idx + ctg2_in_scaffold_idx
            dir_1 = scaffold_length_dict[scaffold_1] - ctg1_in_scaffold_idx + scaffold_length_dict[scaffold_2] - ctg2_in_scaffold_idx
            dir_2 = ctg1_in_scaffold_idx + ctg2_in_scaffold_idx
            dir_3 = ctg1_in_scaffold_idx + scaffold_length_dict[scaffold_2] - ctg2_in_scaffold_idx

            scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][0].append(dir_0)
            scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][1].append(dir_1)
            scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][2].append(dir_2)
            scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][3].append(dir_3)

    elif map_file_type == "pa5":

        with open(map_file, "r") as file:

            for line in file:
                line = line.strip().split()
                if line[0][0] == '#':
                    continue
                ctg1, ctg1_map_idx = line[1], line[2]
                ctg2, ctg2_map_idx = line[3], line[4]

                if ctg1 == ctg2 or ctg2 == "=":
                    continue

                scaffold_1_flag, scaffold_2_flag = False, False
                for (idx_1, idx_2, scaffold) in ctg_scaffold_dict[ctg1]:
                    if int(ctg1_map_idx) > idx_1 and int(ctg1_map_idx) < idx_2:
                        scaffold_1 = scaffold
                        scaffold_1_flag = True
                        break
                for (idx_1, idx_2, scaffold) in ctg_scaffold_dict[ctg2]:
                    if int(ctg2_map_idx) > idx_1 and int(ctg2_map_idx) < idx_2:
                        scaffold_2 = scaffold
                        scaffold_2_flag = True
                        break

                if not scaffold_1_flag or not scaffold_2_flag or scaffold_1 == scaffold_2:
                    continue

                ctg1_in_scaffold_idx = scaffold_ctgs_dict[scaffold_1][ctg1][0] + int(ctg1_map_idx)
                ctg2_in_scaffold_idx = scaffold_ctgs_dict[scaffold_2][ctg2][0] + int(ctg2_map_idx)

                # HT
                if ctg1_in_scaffold_idx * 2 > scaffold_length_dict[scaffold_1]:
                    scaffold_1_pos = "T"
                else:
                    scaffold_1_pos = "H"
                
                if ctg2_in_scaffold_idx * 2 > scaffold_length_dict[scaffold_2]:
                    scaffold_2_pos = "T"
                else:
                    scaffold_2_pos = "H"

                scaffold_HT_dict[tuple(sorted([scaffold_1+"_"+scaffold_1_pos, scaffold_2+"_"+scaffold_2_pos]))] += 1

                # clm
                dir_0 = scaffold_length_dict[scaffold_1] - ctg1_in_scaffold_idx + ctg2_in_scaffold_idx
                dir_1 = scaffold_length_dict[scaffold_1] - ctg1_in_scaffold_idx + scaffold_length_dict[scaffold_2] - ctg2_in_scaffold_idx
                dir_2 = ctg1_in_scaffold_idx + ctg2_in_scaffold_idx
                dir_3 = ctg1_in_scaffold_idx + scaffold_length_dict[scaffold_2] - ctg2_in_scaffold_idx

                scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][0].append(dir_0)
                scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][1].append(dir_1)
                scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][2].append(dir_2)
                scaffold_clm_dict[tuple(sorted([scaffold_1, scaffold_2]))][3].append(dir_3)



    with open(f"{output_prefix}.scaffold.HT.pkl", 'wb') as file:
        pickle.dump(scaffold_HT_dict, file)
    
    with open(f"{output_prefix}.scaffold.clm", 'w') as file:
        for (scaffold_1, scaffold_2) in scaffold_clm_dict:
            for i in range(4):

                if i == 0:
                    scaffold_1_dir = "+"
                    scaffold_2_dir = "+"
                elif i == 1:
                    scaffold_1_dir = "+"
                    scaffold_2_dir = "-"
                elif i == 2:
                    scaffold_1_dir = "-"
                    scaffold_2_dir = "+"
                elif i == 3:
                    scaffold_1_dir = "-"
                    scaffold_2_dir = "-"

                idx_list = [ str(idx) for idx in list(scaffold_clm_dict[(scaffold_1, scaffold_2)][i])]

                file.write(f"{scaffold_1}{scaffold_1_dir} {scaffold_2}{scaffold_2_dir}\t{len(idx_list)}\t{' '.join(idx_list)}\n")


    return scaffold_HT_dict, scaffold_clm_dict

def main():
    parser = argparse.ArgumentParser(description="Gets the data for running HapHiC sort.")

    # Required arguments
    parser.add_argument("-m", "--map_file", required=True, help="mapping file.")
    parser.add_argument("--map_file_type", default="pa5", choices=['bam', 'pa5'], type=str.lower, help="mapping file type.")
    parser.add_argument("-a", "--agp", required=True, help="agp files.")
    parser.add_argument("-r", "--RE_file", required=True, help="Path to the restriction enzyme file.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files.")


    # pairs_file = "chr1g1.pairs"
    # agp_file = "subgraphGroup.agp"
    # RE_file = "rice4.RE_counts.txt"
    
    args = parser.parse_args()
    map_file = args.map_file
    map_file_type = args.map_file_type
    agp_file = args.agp
    RE_file = args.RE_file
    output_prefix = args.output_prefix


    scaffold_ctgs_dict, ctg_scaffold_dict, scaffold_length_dict = read_agp(agp_file)
    read_map_file(map_file, map_file_type, scaffold_ctgs_dict, ctg_scaffold_dict, scaffold_length_dict, output_prefix)
    get_RE(RE_file, scaffold_ctgs_dict, scaffold_length_dict, output_prefix)


if __name__ == "__main__":
    main()