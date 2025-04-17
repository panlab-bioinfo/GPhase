from Bio import SeqIO
import argcomplete
from argcomplete.completers import FilesCompleter
import argparse

def count_restriction_sites(fasta_file, enzyme_site):

    result = []

    enzyme_site = enzyme_site.lower()

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).lower()   
        count = int(sequence.count(enzyme_site)) + 1     
        seq_length = len(sequence)
        result.append((record.id, count, seq_length))

    return result

def write_output_to_file(results, output_prefix):


    with open(f"{output_prefix}.RE_counts.txt", "w") as f:
        f.write("#Contig\tRECounts\tLength\n")

        for seq_id, count, seq_length in results:
            f.write(f"{seq_id}\t{count}\t{seq_length}\n")

def main():
    parser = argparse.ArgumentParser(description="Gets the RE file form fasta.")
    parser.add_argument("-f", "--fasta", required=True, help="fasta file.")
    parser.add_argument("-e", "--enzyme_site", default="GATC", help="restriction enzyme file.")
    parser.add_argument("-op", "--output_prefix", required=True, help="Prefix for output files.")



    args = parser.parse_args()
    argcomplete.autocomplete(parser)
    fasta_file = args.fasta
    enzyme_site = args.enzyme_site
    output_prefix = args.output_prefix


    results = count_restriction_sites(fasta_file, enzyme_site)


    write_output_to_file(results, output_prefix)

if __name__ == "__main__":
    main()
