from collections import defaultdict
import pandas as pd
import argparse
from pathlib import Path
import sys


def read_noNorFile(noNorFile):
    utg_utg_link_dict = defaultdict()
    with open(noNorFile, "r") as file:
        for line in file:
            if line != "source,target,links\n":
                line = line.strip().split(',')
                utg_utg_link_dict[tuple([line[0],  line[1]])] = int(line[2])
    return utg_utg_link_dict


def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len



def normalize_links(utg_utg_link_dict, ctg_RE_len_dict, output_file):
    with Path(output_file).open("w") as file:
        # if head == True:
        file.write("source,target,links\n")
        for pair, links in utg_utg_link_dict.items():

            # r1 = float(ctg_RE_len_dict[pair[0]][0])
            # r2 = float(ctg_RE_len_dict[pair[1]][0])


            # r1 = float(ctg_RE_len_dict[pair[0]][0]**2 / ctg_RE_len_dict[pair[0]][1]) 
            # r2 = float(ctg_RE_len_dict[pair[1]][0]**2 / ctg_RE_len_dict[pair[1]][1])

            r1 = float(ctg_RE_len_dict[pair[0]][0] / ctg_RE_len_dict[pair[0]][1]) 
            r2 = float(ctg_RE_len_dict[pair[1]][0] / ctg_RE_len_dict[pair[1]][1])

            # r1 = float(ctg_RE_len_dict[pair[0]][0] * ctg_RE_len_dict[pair[0]][1]) 
            # r2 = float(ctg_RE_len_dict[pair[1]][0] * ctg_RE_len_dict[pair[1]][1])

            links /= (r1 * r2)
            links /= 1e6
            line = pair[0] + "," + pair[1] + "," + str(links) + "\n"
            file.write(line)

def main():
    parser = argparse.ArgumentParser(description="Standardize the hic signal file according to the RE of the ctg", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--file',required=True, help='HiC links file')
    parser.add_argument('-r', '--REfile',required=True, help='REs file')
    parser.add_argument('-o', '--output', required=True,help='Output file name')
    args = parser.parse_args()
    csv_file = args.file
    res_file = args.REfile
    output_file = args.output

    utg_utg_link_dict = read_noNorFile(csv_file)
    ctg_RE_len_dict = read_REs(res_file)
    normalize_links(utg_utg_link_dict, ctg_RE_len_dict, output_file)



if __name__=="__main__":
    main()


