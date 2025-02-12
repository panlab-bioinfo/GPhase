from Bio import SeqIO
import argcomplete
from argcomplete.completers import FilesCompleter
import argparse

def count_restriction_sites(fasta_file, enzyme_site):
    """
    统计每个序列中的内切酶位点数目，并返回包含序列名称、内切酶位点数目和序列长度的列表。
    
    参数:
    fasta_file (str): FASTA 文件路径
    enzyme_site (str): 内切酶识别位点的序列
    
    返回:
    list: 包含每条序列名称、内切酶位点数目和序列长度的元组列表
    """
    result = []

    enzyme_site = enzyme_site.lower()

    for record in SeqIO.parse(fasta_file, "fasta"):
        # print(record)
        sequence = str(record.seq).lower()   
        count = int(sequence.count(enzyme_site)) + 1     
        seq_length = len(sequence)
        result.append((record.id, count, seq_length))

    return result

def write_output_to_file(results, output_prefix):
    """
    将结果写入CSV文件
    
    参数:
    results (list): 包含每条序列名称、内切酶位点数目和序列长度的元组列表
    output_file (str): 输出文件路径
    """
    with open(f"{output_prefix}.RE_counts.txt", "w") as f:
        f.write("#Contig\tRECounts\tLength\n")

        for seq_id, count, seq_length in results:
            f.write(f"{seq_id}\t{count}\t{seq_length}\n")

def main():
    parser = argparse.ArgumentParser(description="Gets the RE file form fasta.")
    parser.add_argument("-f", "--fasta", required=True, help="fasta file.")
    parser.add_argument("-e", "--enzyme_site", default="GATC", help="restriction enzyme file.")
    parser.add_argument("-op", "--output_prefix", required=True, help="Prefix for output files.")

    # fasta_file = "asm.fa"  # 替换为你的FASTA文件路径
    # enzyme_site = "GATC"  # 替换为你感兴趣的内切酶位点（如 EcoRI 的 GAATTC）
    # output_file = "wax.test.counts.txt"  # 输出文件路径

    args = parser.parse_args()
    argcomplete.autocomplete(parser)
    fasta_file = args.fasta
    enzyme_site = args.enzyme_site
    output_prefix = args.output_prefix


    # 获取每个序列中的内切酶位点数目
    results = count_restriction_sites(fasta_file, enzyme_site)

    # 将结果写入CSV文件
    write_output_to_file(results, output_prefix)

if __name__ == "__main__":
    main()
