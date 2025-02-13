from collections import defaultdict
import argparse

def read_fai(fai_file):
    fai_dict = defaultdict()
    fai_reverse_dict = defaultdict(str)

    with open(fai_file, 'r') as file:
        idx = 1
        for line in file:
            line = line.strip().split()
            # if line[0].startswith('u') or line[0].startswith('c'):
            fai_dict[str(idx)] = line[0]
            fai_reverse_dict[line[0]] = idx
            idx += 1

    return fai_dict, fai_reverse_dict

def read_partig(partig_file, fai_dict, fai_reverse_dict, output_file):
    
    partig_re = defaultdict()

    with open(partig_file, 'r') as file, open(output_file, 'w') as output_file:
        for line in file:
            line = line.strip().split()
            if line[0] == "S":
                scontig1 = line[1][1:]
                scontig2 = line[2][1:]
                contig1 = fai_dict[scontig1]
                contig2 = fai_dict[scontig2]

                output_file.write(f"{contig1},{contig2},{line[7]}\n")


def main():
    parser = argparse.ArgumentParser(description="trans partig file", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-fai', '--fai', required=True, help='fai for asm.fa')
    parser.add_argument('-p', '--partig_file', required=True, help='partig file')
    parser.add_argument('-o', '--output', required=True, help='output file')

    args = parser.parse_args()  
    fai_file = args.fai
    partig_file = args.partig_file
    output_file = args.output

    fai_dict, fai_reverse_dict = read_fai(fai_file)
    read_partig(partig_file, fai_dict, fai_reverse_dict, output_file)

if __name__ == "__main__":
    main()
