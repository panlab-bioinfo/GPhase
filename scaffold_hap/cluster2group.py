import argparse
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Process cluster and counts files.')
    parser.add_argument('-c','--cluster', required=True, help='Path to the cluster file')
    parser.add_argument('-r', '--re', required=True, help='REs file')
    parser.add_argument('-m', '--mark',default="",help='group_name_[mark].txt')

    args = parser.parse_args()


    counts_dict = defaultdict(list)
    with open(args.re, 'r') as f:
        for line in f:
            line = line.rstrip('\n')  
            if line.strip():
                parts = line.split()
                if parts:
                    id = parts[0]
                    counts_dict[id].append(line)

    with open(args.cluster, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            group_name = parts[0]
            ids = parts[2].split()

            matched_lines = []
            for id in ids:
                matched_lines.extend(counts_dict.get(id, []))

            if matched_lines:
                if not args.mark:
                    with open(f"{group_name}.txt", 'w') as out_f:
                        matched_lines = list(set(matched_lines))
                        out_f.write('\n'.join(matched_lines) + '\n')
                else:
                    with open(f"{group_name}_{args.mark}.txt", 'w') as out_f:
                        matched_lines = list(set(matched_lines))
                        out_f.write('\n'.join(matched_lines) + '\n')


if __name__ == "__main__":
    main()