import argparse
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import igraph as ig
import re



def read_RE(REFile):
    """读取RE文件，返回 {ctg: (RE1, RE2)}"""
    ctg_RE_dict = {}
    with open(REFile, 'r') as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                ctg_RE_dict[parts[0]] = (int(parts[1]), int(parts[2]))
    return ctg_RE_dict

def read_gfa(file_path):
    """读取GFA文件，返回 dict，记录 forward/reverse list"""
    gfa_graph = defaultdict(lambda: {"forward_list": [], "reverse_list": [], "edges": {}})
    utgs_set = set()
    with open(file_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[0] == 'S':
                _ = gfa_graph[fields[1]] 
                utgs_set.add(fields[1])
            elif fields[0] == 'L':
                from_node, from_orient, to_node, to_orient, overlap = (
                    fields[1], fields[2], fields[3], fields[4], fields[5][:-1]
                )
                if from_orient == '+':
                    gfa_graph[from_node]['forward_list'].append((to_node, to_orient))
                else:
                    gfa_graph[from_node]['reverse_list'].append((to_node, to_orient))

                overlap_length = int(''.join(filter(str.isdigit, overlap)))

                gfa_graph[from_node]['edges'][(to_node, to_orient)] = overlap_length

    return gfa_graph, utgs_set

def read_graph_igraph(edge_file):
    """读取边文件并构建 igraph Graph"""
    edges = []
    nodes = set()
    with open(edge_file, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 2:
                u, v = parts[0], parts[1]
                edges.append((u, v))
                nodes.add(u)
                nodes.add(v)
    node_list = sorted(nodes)
    name_to_idx = {name: idx for idx, name in enumerate(node_list)}
    edges_idx = [(name_to_idx[u], name_to_idx[v]) for u, v in edges]
    g = ig.Graph(edges=edges_idx, directed=True)
    g.vs["name"] = node_list
    return g, name_to_idx

def read_agp(agp_file):
    """读取AGP，返回 {scaffold: [(utg, dir), ...]}"""
    scaffolds = defaultdict(list)
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[4] == 'W':
                scaffolds[parts[0]].append((parts[5], parts[8]))
    return scaffolds

def read_agp_pd(agp_file):
    return pd.read_csv(
        agp_file,
        sep='\t',
        header=None,
        comment='#',
        names=['scaffold', 'start', 'end', 'part_num', 'type',
               'object', 'object_beg', 'object_end', 'orientation']
    )


def has_single_shortest_path_igraph(g, name_to_idx, source, target, shortest_len_T=10):
    src_idx = name_to_idx.get(source)
    tgt_idx = name_to_idx.get(target)
    if src_idx is None or tgt_idx is None:
        return False
    dist = g.shortest_paths_dijkstra(src_idx, tgt_idx)[0][0]
    if dist == float("inf") or dist > shortest_len_T:
        return False
    paths = g.get_all_shortest_paths(src_idx, tgt_idx)
    return len(paths) == 1


def rescue(edge_file, agp_file, gfa_file, REFile):
    utg_rescue_dict = {}

    g, name_to_idx = read_graph_igraph(edge_file)
    scaffolds = read_agp(agp_file)
    gfa_graph, utgs_set = read_gfa(gfa_file)
    ctg_RE_dict = read_RE(REFile)

    for scaffold_id, utgs in scaffolds.items():
        for i in range(len(utgs) - 1):
            (utg1, dir1), (utg2, dir2) = utgs[i], utgs[i + 1]
            if utg1 not in name_to_idx or utg2 not in name_to_idx:
                continue

            if has_single_shortest_path_igraph(g, name_to_idx, utg1, utg2):
                try:
                    vpaths = g.get_shortest_paths(name_to_idx[utg1], to=name_to_idx[utg2], output="vpath")[0]
                    path_names = [g.vs[idx]["name"] for idx in vpaths]

                    path_utg_dir_list = [(utg1, dir1)]
                    before_utg, before_dir = utg1, dir1
                    for _utg in path_names[1:]:
                        if _utg == utg2:
                            path_utg_dir_list.append((utg2, dir2))
                        else:
                            if any(x[0] == _utg for x in gfa_graph[before_utg]['forward_list']):
                                before_dir = "+"
                            elif any(x[0] == _utg for x in gfa_graph[before_utg]['reverse_list']):
                                before_dir = "-"
                            path_utg_dir_list.append((_utg, before_dir))
                            before_utg = _utg

                    utg_rescue_dict[(utg1, utg2)] = path_utg_dir_list
                except Exception:
                    continue
    return utg_rescue_dict


def update_agp_with_insert_lists(agp_df, insert_dict, ctg_RE_dict, utgs_set):
    updated_rows = []
    grouped = agp_df.groupby("scaffold", sort=False)

    for scaffold, group in grouped:
        before_utg, bp_set_off, idx_set_off = None, 0, 0
        group = group.sort_values(by="start").reset_index(drop=True)
        new_group = []
        i = 0

        while i < len(group):
            row = group.iloc[i]

            if row['type'] == 'W':
                if not before_utg:
                    new_group.append(row.to_dict())
                    before_utg = row['object']
                    i += 1
                    continue
                else:
                    now_utg = row['object']
                    key = (before_utg, now_utg)

                    if (now_utg not in ctg_RE_dict) or (key not in insert_dict):
                        new_row = row.copy()
                        new_row['start'] = int(row['start']) + bp_set_off
                        new_row['end'] = int(row['end']) + bp_set_off
                        new_row['part_num'] = int(row['part_num']) + idx_set_off
                        new_group.append(new_row.to_dict())
                        before_utg = now_utg
                        i += 1
                        continue

                    # insert insert_path
                    path = insert_dict[key]
                    insert_path = path[1:-1]

                    for idx, _utg in enumerate(insert_path):

                        if _utg not in utgs_set:
                            continue
                        insert_begin = int(row["start"]) + bp_set_off

                        new_group.append({
                            "scaffold": scaffold,
                            "start": insert_begin,
                            "end": insert_begin + ctg_RE_dict[_utg[0]][1],
                            "part_num": int(row['part_num']) + idx_set_off,
                            "type": 'W',
                            "object": _utg[0],
                            "object_beg": 1,
                            "object_end": ctg_RE_dict[_utg[0]][1],
                            "orientation": _utg[1]
                        })
                        bp_set_off += ctg_RE_dict[_utg[0]][1]
                        idx_set_off += 1

                        new_group.append({
                            "scaffold": scaffold,
                            "start": insert_begin + ctg_RE_dict[_utg[0]][1] + 1,
                            "end": insert_begin + ctg_RE_dict[_utg[0]][1] + 100,
                            "part_num": int(row['part_num']) + idx_set_off,
                            "type": 'U',
                            "object": 100,
                            "object_beg": 'scaffold',
                            "object_end": 'yes',
                            "orientation": 'proximity_ligation'
                        })
                        bp_set_off += 100
                        idx_set_off += 1

                    new_row = row.copy()
                    new_row['start'] = int(row['start']) + bp_set_off + 1
                    new_row['end'] = int(row['end']) + bp_set_off
                    new_row['part_num'] = int(row['part_num']) + idx_set_off
                    new_group.append(new_row.to_dict())

                    before_utg = now_utg
                    i += 1
                    continue

            elif row['type'] in ['U', 'N']:
                new_group.append({
                    "scaffold": scaffold,
                    "start": int(row['start']) + bp_set_off,
                    "end": int(row['end']) + bp_set_off,
                    "part_num": int(row['part_num']) + idx_set_off,
                    "type": row['type'],
                    "object": row['object'],
                    "object_beg": row['object_beg'],
                    "object_end": row['object_end'],
                    "orientation": row['orientation']
                })
                i += 1
                continue
            else:
                i += 1

        updated_rows.extend(new_group)

    updated_df = pd.DataFrame(updated_rows, columns=agp_df.columns)
    updated_df.to_csv("gphase_final_rescue.agp", sep='\t', index=False, header=False)
    return updated_df


def scaffold_sequences_from_agp(updated_df, gfa_graph, fasta_file):
    def reverse_complement(seq):
        return str(Seq(seq).reverse_complement())

    utg_seq_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}
    scaffold_seq_dict = {}
    grouped = updated_df.groupby("scaffold", sort=False)

    for scaffold, group in grouped:
        group = group.reset_index(drop=True)
        i = 0
        while i < len(group):
            row_current = group.iloc[i]
            if row_current['type'] != 'W':
                i += 1
                continue

            utg_id = row_current['object']
            orientation = row_current['orientation']
            seq = utg_seq_dict[utg_id]

            start_pos = int(row_current.iloc[6])  
            end_pos = int(row_current.iloc[7])   

            if start_pos < 1 or end_pos < start_pos or end_pos > len(seq):
                raise ValueError(f"error: start_pos={start_pos}, end_pos={end_pos} for {utg_id}")
            actual_length = end_pos - start_pos + 1
            seq = seq[start_pos-1:end_pos]  
            
            if orientation == '-':
                seq = reverse_complement(seq)
            current_seq = seq
            current_components = [f"{utg_id}_{orientation}"]

            j = i + 1
            while j < len(group):
                row_next = group.iloc[j]
                if row_next['type'] != 'W':
                    j += 1
                    continue

                next_utg_id = row_next['object']
                next_orientation = row_next['orientation']
                connected = False
                if orientation == '+':
                    if (next_utg_id, next_orientation) in gfa_graph[utg_id]['forward_list']:
                        connected = True
                        overlap_length = gfa_graph[utg_id]['edges'].get((next_utg_id, next_orientation), 0)
                elif orientation == '-':
                    if (next_utg_id, next_orientation) in gfa_graph[utg_id]['reverse_list']:
                        connected = True
                        overlap_length = gfa_graph[utg_id]['edges'].get((next_utg_id, next_orientation), 0)

                if connected:
                    next_seq = utg_seq_dict[next_utg_id]
                    # 使用 updated_df 第七列和第八列裁剪实际序列
                    try:
                        start_pos = int(row_next.iloc[6])  
                        end_pos = int(row_next.iloc[7])   
                        if start_pos < 1 or end_pos < start_pos or end_pos > len(next_seq):
                            raise ValueError(f"error: start_pos={start_pos}, end_pos={end_pos} for {next_utg_id}")
                        actual_length = end_pos - start_pos + 1
                        next_seq = next_seq[start_pos-1:end_pos]  
                    except (ValueError, IndexError) as e:
                        actual_length = len(next_seq)

                    # 根据方向转换
                    if next_orientation == '-':
                        next_seq = reverse_complement(next_seq)

                    overlap_length = gfa_graph[utg_id]['edges'].get((next_utg_id, next_orientation), 0)

                    if overlap_length > 0:
                        if overlap_length >= len(next_seq):
                            overlap_length = 0
                        else:
                            next_seq = next_seq[overlap_length:] 

                    current_seq += next_seq
                    current_components.append(f"{next_utg_id}_{next_orientation}")
                    utg_id = next_utg_id
                    orientation = next_orientation
                    j += 1
                else:
                    break

            scaffold_name = "_".join(current_components)
            if scaffold_name:
                scaffold_seq_dict[scaffold_name] = current_seq

            i = j if j > i else i + 1

    # 添加未在 AGP 的 contig
    agp_utg_ids = set(updated_df[updated_df['type'] == 'W']['object'])
    for utg_id in utg_seq_dict:
        if utg_id not in agp_utg_ids:
            scaffold_seq_dict[f"{utg_id}_+"] = utg_seq_dict[utg_id]

    # 输出 fasta
    with open("gphase_final_contig.fasta", 'w') as f,  open("gphase_final_ctg2utg.txt", 'w') as f2:
        for idx, (scaffold_name, seq) in enumerate(scaffold_seq_dict.items(), 1):

            ctg_ID = f"ctg{idx:0>6}l"
            utgs_list = scaffold_name.replace("+", "_").replace("-", "_").split("_")

            f2.write(f"{ctg_ID}\t{scaffold_name}\n")

            f.write(f">{ctg_ID}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    return scaffold_seq_dict


def utg_to_ctg_agp(updated_df, ctg2utg_file, fasta_file, output_agp="gphase_final_contig.agp", gap_len=100):

    ctg_len = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    utg_to_ctg_list = defaultdict(list)  
    pair_pat = re.compile(r"(utg\d+[a-z]?)_([+-])")
    with open(ctg2utg_file) as f:
        for line in f:
            ctg, utgs_str = line.strip().split("\t")
            pairs = pair_pat.findall(utgs_str)
            for utg_id, orient in pairs:
                utg_to_ctg_list[utg_id].append(ctg)

    utg_counter = defaultdict(int)
    agp_lines = []

    for scaffold, group in updated_df.groupby("scaffold", sort=False):
        group = group.reset_index(drop=True)
        part_number = 1
        current_pos = 1
        i = 0
        prev_ctg = None

        while i < len(group):
            row = group.iloc[i]

            if row['type'] != 'W':
                i += 1
                continue 

            utg = row['object']
            idx = utg_counter[utg]
            if idx >= len(utg_to_ctg_list[utg]):
                idx = len(utg_to_ctg_list[utg]) - 1
            ctg = utg_to_ctg_list[utg][idx]
            utg_counter[utg] += 1


            j = i + 1
            while j < len(group):
                nxt_row = group.iloc[j]
                if nxt_row['type'] == 'W':
                    nxt_utg = nxt_row['object']
                    nxt_idx = utg_counter[nxt_utg]
                    if nxt_idx >= len(utg_to_ctg_list[nxt_utg]):
                        nxt_idx = len(utg_to_ctg_list[nxt_utg]) - 1
                    nxt_ctg = utg_to_ctg_list[nxt_utg][nxt_idx]
                    if nxt_ctg != ctg:
                        break
                    utg_counter[nxt_utg] += 1
                j += 1

            if prev_ctg is not None and prev_ctg != ctg:
                agp_lines.append([
                    scaffold,
                    current_pos,
                    current_pos + gap_len - 1,
                    part_number,
                    "U",
                    gap_len,
                    "scaffold",
                    "yes",
                    "proximity_ligation"
                ])
                current_pos += gap_len
                part_number += 1


            length = ctg_len.get(ctg)
            if length is None:
                raise KeyError(f"CTG '{ctg}' 未在 fasta 文件中找到！")
            agp_lines.append([
                scaffold,
                current_pos,
                current_pos + length - 1,
                part_number,
                "W",
                ctg,
                1,
                length,
                "+"
            ])
            current_pos += length
            part_number += 1

            prev_ctg = ctg
            i = j

    with open(output_agp, "w") as f:
        for line in agp_lines:
            f.write("\t".join(map(str, line)) + "\n")

    return output_agp


def Rescue_base_graph(edge_file, agp_file, gfa_file, REFile, fa_file):
    gfa_graph, utgs_set = read_gfa(gfa_file)
    ctg_RE_dict = read_RE(REFile)
    utg_rescue_dict = rescue(edge_file, agp_file, gfa_file, REFile)
    agp_df = read_agp_pd(agp_file)
    updated_df = update_agp_with_insert_lists(agp_df, utg_rescue_dict, ctg_RE_dict, utgs_set)
    scaffold_seq_dict = scaffold_sequences_from_agp(updated_df, gfa_graph, fa_file)

    utg_to_ctg_agp(updated_df, "gphase_final_ctg2utg.txt", "gphase_final_contig.fasta", output_agp="gphase_final_contig.agp")




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rescue scaffolds and generate sequences based on AGP and GFA files.')
    parser.add_argument('-e', '--edge_file', required=True, help='Input edge file')
    parser.add_argument('-a', '--agp_file', required=True, help='Input AGP file')
    parser.add_argument('-g', '--gfa_file', required=True, help='Input GFA file')
    parser.add_argument('-r', '--re_file', required=True, help='Input RE counts file')
    parser.add_argument('-f', '--fasta_file', required=True, help='Input fasta file of unitigs')
    args = parser.parse_args()

    Rescue_base_graph(
        args.edge_file,
        args.agp_file,
        args.gfa_file,
        args.re_file,
        args.fasta_file
    )
