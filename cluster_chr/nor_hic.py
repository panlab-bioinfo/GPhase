from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import statsmodels.api as sm
import sys


def read_noNorFile(noNorFile):
    utg_utg_link_dict = defaultdict()
    with open(noNorFile, "r") as file:
        for line in file:
            if line != "source,target,links\n":
                line = line.strip().split(',')
                utg_utg_link_dict[tuple([line[0],  line[1]])] = float(line[2])
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

def normalize_links_v2(utg_utg_link_dict, ctg_RE_len_dict, output_file):
    rows = []

    for (u, v), links in utg_utg_link_dict.items():
        len1, re1 = ctg_RE_len_dict[u]
        len2, re2 = ctg_RE_len_dict[v]
        
        rows.append({
            "source": u,
            "target": v,
            "log_len_prod": np.log(len1 * len2 + 1),
            "log_re_prod": np.log(re1 * re2 + 1),
            "log_links": np.log(links + 1)
        })

    df = pd.DataFrame(rows)

    X = df[["log_len_prod", "log_re_prod"]]
    X = sm.add_constant(X)  
    y = df["log_links"]
    model = sm.OLS(y, X).fit()

    df["log_pred"] = model.predict(X)
    df["residual"] = df["log_links"] - df["log_pred"]
    df["exp_residual"] = np.exp(df["residual"])  

    with Path(output_file).open("w") as file:
        file.write("source,target,links\n")
        for _, row in df.iterrows():
            file.write(f"{row['source']},{row['target']},{row['exp_residual']:.6f}\n")

    # print(model.summary())


def normalize_links(utg_utg_link_dict, ctg_RE_len_dict, output_file):
    with Path(output_file).open("w") as file:
        # if head == True:
        file.write("source,target,links\n")
        for pair, links in utg_utg_link_dict.items():

            # r1 = float(ctg_RE_len_dict[pair[0]][0])
            # r2 = float(ctg_RE_len_dict[pair[1]][0])

            if pair[0] in ctg_RE_len_dict and pair[1] in ctg_RE_len_dict:
                r1 = float(ctg_RE_len_dict[pair[0]][0] / ctg_RE_len_dict[pair[0]][1]) 
                r2 = float(ctg_RE_len_dict[pair[1]][0] / ctg_RE_len_dict[pair[1]][1])
            else:
                continue
            

            # r1 = float(ctg_RE_len_dict[pair[0]][0] * ctg_RE_len_dict[pair[0]][1]) 
            # r2 = float(ctg_RE_len_dict[pair[1]][0] * ctg_RE_len_dict[pair[1]][1])

            links /= (r1 * r2)
            links /= 1e6
            # links *= 1e6

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


