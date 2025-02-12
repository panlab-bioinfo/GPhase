import argparse

def agp_to_fasta_with_gaps(agp_file, sequence_dict, output_file):
    scaffold_sequences = {}

    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith("#"):  # 跳过注释行
                continue
            
            # 解析 AGP 文件中的每一行
            columns = line.strip().split('\t')
            scaffold_name = columns[0]
            fragment_start = int(columns[1])
            fragment_end = int(columns[2])
            fragment_type = columns[4]
            fragment_sequence_name = columns[5] if fragment_type == "W" else None
            gap_length = fragment_end - fragment_start
            
            # 如果是已知序列（W 类型），则从 sequence_dict 获取序列
            if fragment_type == "W" and fragment_sequence_name:
                fragment_sequence = sequence_dict.get(fragment_sequence_name, "")
            # 如果是 gap 区域（U 类型），则插入 N
            elif fragment_type == "U":
                fragment_sequence = "N" * (gap_length+1)
            else:
                continue  # 跳过不相关的行

            # 如果 scaffold 不在字典中，初始化一个空序列
            if scaffold_name not in scaffold_sequences:
                scaffold_sequences[scaffold_name] = ""

            # 拼接片段序列到对应 scaffold
            scaffold_sequences[scaffold_name] += fragment_sequence

    # 将拼接的 scaffold 序列写入 FASTA 文件
    with open(output_file, 'w') as out_f:
        for scaffold, sequence in scaffold_sequences.items():
            out_f.write(f">{scaffold}\n{sequence}\n")


def load_sequence_dict(sequence_file):
    """加载序列字典文件，假设每行是一个 FASTA 格式的条目"""
    sequence_dict = {}
    with open(sequence_file, 'r') as seq_file:
        seq_name = None
        seq_sequence = []
        for line in seq_file:
            line = line.strip()
            if line.startswith(">"):  # 新的序列
                if seq_name:
                    sequence_dict[seq_name] = "".join(seq_sequence)
                seq_name = line[1:]
                seq_sequence = []
            else:
                seq_sequence.append(line)
        if seq_name:
            sequence_dict[seq_name] = "".join(seq_sequence)  # 最后一条序列
    return sequence_dict




def main():
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="Convert AGP file to FASTA format.")
    parser.add_argument("agp_file", help="Path to the AGP file")
    parser.add_argument("fasta_file", help="Path to the reference FASTA file")
    parser.add_argument("output_fasta", help="Path to the output FASTA file")

    # 解析命令行参数
    args = parser.parse_args()

    # 加载序列字典
    sequence_dict = load_sequence_dict(args.fasta_file)

    # 调用转换函数
    agp_to_fasta_with_gaps(args.agp_file, sequence_dict, args.output_fasta)


if __name__ == "__main__":
    main()
