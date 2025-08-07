import networkx as nx
import argparse
import copy as cp
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq




def has_single_shortest_path(G, source, target, weight=None, shortest_len_T=10):

    try:
        shortest_length = nx.shortest_path_length(G, source, target, weight=weight)

        if shortest_length > shortest_len_T:
            return False

        path_count = [0]  

        def dfs(current, target, length_so_far, visited):
            if path_count[0] > 1:  
                return
            if current == target:
                if length_so_far == shortest_length:  
                    path_count[0] += 1
                return
            if length_so_far > shortest_length:  
                return


            for neighbor in G.successors(current):
                if neighbor not in visited:
                    edge_weight = 1 if weight is None else G[current][neighbor][weight]
                    new_length = length_so_far + edge_weight
                    visited.add(neighbor)
                    dfs(neighbor, target, new_length, visited)
                    visited.remove(neighbor)

        dfs(source, target, 0, {source})
        return path_count[0] == 1

    except nx.NetworkXNoPath:
        return False

def read_RE(REFile):
    ctg_RE_dict = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_dict

def read_gfa(file_path):
    G = nx.DiGraph()
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  

            fields = line.split('\t')
            record_type = fields[0]

            if record_type == 'S':
                node_id = fields[1]
                G.add_node(
                    node_id,
                    forward_list=[],
                    reverse_list=[]
                )

            elif record_type == 'L':
                from_node = fields[1]
                from_orient = fields[2]
                to_node = fields[3]
                to_orient = fields[4]
                overlap = fields[5][:-1]
                G.add_edge(
                    from_node,
                    to_node,
                    from_orient=from_orient,
                    to_orient=to_orient,
                    overlap=overlap
                )
                pair = (to_node, to_orient)
                if from_orient == '+':
                    G.nodes[from_node]['forward_list'].append(pair)
                elif from_orient == '-':
                    G.nodes[from_node]['reverse_list'].append(pair)
    return G

def read_graph(edge_file):
    G = nx.DiGraph()
    with open(edge_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 2:
                utgA, utgB = parts[0], parts[1]
                G.add_edge(utgA, utgB)
    
    G_reverse = nx.DiGraph.reverse(G, copy=True)
    return G, G_reverse


def read_agp(agp_file):
    scaffolds = defaultdict(list)
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6 or parts[4] != 'W':
                continue
            scaffold_id = parts[0]
            utg_id = parts[5]
            direction = parts[8]
            scaffolds[scaffold_id].append((utg_id, direction))
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


def rescue(edge_file, agp_file, gfa_file, REFile):

    utg_rescue_dict = defaultdict(list)

    G, G_reverse = read_graph(edge_file)
    scaffolds = read_agp(agp_file)
    gfa_graph = read_gfa(gfa_file)
    ctg_RE_dict = read_RE(REFile)
    find_path_graph = cp.deepcopy(G)


    used_utgs = set() 
    utg_rescue_dict = defaultdict(list)

    for scaffold_id, utgs in scaffolds.items():

        for i in range(len(utgs) - 1):
            (utg1,dir1), (utg2, dir2) = utgs[i], utgs[i + 1]

            if utg1 in used_utgs or utg2 in used_utgs:
                continue

            if utg1 in list(G.nodes()) and utg2 in list(G.nodes()):

                if dir1 == "+":
                    successors = list(G.successors(utg1))
                    utg_forward_set = set(pair[0] for pair in gfa_graph.nodes()[utg1]['forward_list'])
                    if set(successors) & utg_forward_set:
                        find_path_graph = cp.deepcopy(G)
                    else:
                        find_path_graph = cp.deepcopy(G_reverse)

                else:
                    successors = list(G.successors(utg1))
                    utg_reverse_set = set(pair[0] for pair in gfa_graph.nodes()[utg1]['reverse_list'])
                    if set(successors) & utg_reverse_set:
                        find_path_graph = cp.deepcopy(G_reverse)
                    else:
                        find_path_graph = cp.deepcopy(G)


                # if only one path
                if has_single_shortest_path(find_path_graph, source=utg1, target=utg2, weight=None):

                    try:
                        path = nx.shortest_path(find_path_graph, source=utg1, target=utg2)
                        path_utg_dir_list = [tuple([utg1, dir1])]
                        before_utg, before_dir = utg1, dir1

                        if any(utg in used_utgs for utg in path[1:-1]):
                            continue

                        for _utg in path[1:]:
                            if _utg not in ctg_RE_dict:
                                break
                            if _utg == utg2:
                                pair = (utg2, dir2)
                            else:
                                try:
                                    utg_forward_set = set(pair[0] for pair in gfa_graph.nodes()[before_utg]['forward_list'])
                                    utg_reverse_set = set(pair[0] for pair in gfa_graph.nodes()[before_utg]['reverse_list'])
                                    if gfa_graph.has_edge(before_utg, _utg):
                                        before_dir = gfa_graph.edges[before_utg, _utg]['to_orient']
                                    pair = (_utg, before_dir)
                                    before_utg = _utg
                                except:
                                    continue
                            path_utg_dir_list.append(pair)

                        if path_utg_dir_list:     
                            utg_rescue_dict[(utg1, utg2)] = path_utg_dir_list
                    except:
                        print("error...")

    return utg_rescue_dict

def update_agp_with_insert_lists(agp_df, insert_dict, ctg_RE_dict):
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


def scaffold_sequences_from_agp(updated_df, G, fasta_file):

    def reverse_complement(seq):
        return str(Seq(seq).reverse_complement())

    utg_seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        utg_seq_dict[record.id] = str(record.seq)

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
                    if (next_utg_id, next_orientation) in G.nodes[utg_id]['forward_list']:
                        connected = True
                elif orientation == '-':
                    if (next_utg_id, next_orientation) in G.nodes[utg_id]['reverse_list']:
                        connected = True

                if connected and G.has_edge(utg_id, next_utg_id):
                    overlap_len = int(G.edges[utg_id, next_utg_id]['overlap'])
                    next_seq = utg_seq_dict[next_utg_id]
                    if next_orientation == '-':
                        next_seq = reverse_complement(next_seq)

                    # 验证overlap区域
                    if len(current_seq) >= overlap_len and len(next_seq) >= overlap_len:
                        current_tail = current_seq[-overlap_len:]
                        next_head = next_seq[:overlap_len]
                        if current_tail == next_head:
                            current_seq = current_seq[:-overlap_len] + next_seq
                        else:
                            current_seq = current_seq[:-overlap_len] + next_seq
                    else:
                        current_seq = current_seq[:-overlap_len] + next_seq

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

    agp_utg_ids = set(updated_df[updated_df['type'] == 'W']['object'])
    for utg_id in utg_seq_dict:
        if utg_id not in agp_utg_ids:
            scaffold_name = f"{utg_id}_+"
            scaffold_seq_dict[scaffold_name] = utg_seq_dict[utg_id]

    with open("gphase_final_contig.fasta", 'w') as f:
        for scaffold_name, seq in scaffold_seq_dict.items():
            f.write(f">{scaffold_name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    return scaffold_seq_dict


def Rescue_base_graph(edge_file, agp_file, gfa_file, REFile, fa_file):
    G = read_gfa(gfa_file)
    ctg_RE_dict = read_RE(REFile)
    utg_rescue_dict = rescue(edge_file, agp_file, gfa_file, REFile)
    
    agp_df = read_agp_pd(agp_file)

    updated_df = update_agp_with_insert_lists(agp_df, utg_rescue_dict, ctg_RE_dict)

    scaffold_sequences_from_agp(updated_df, G, fa_file)







if __name__ == '__main__':
    # edge_file = 'waxapple.digraph.csv' 
    # agp_file = 'gphase_final.agp'  
    # gfa_file = "waxapple.ho.asm.bp.p_utg.gfa"
    # REFile = "waxapple.RE_counts.txt"
    # fa_file = "waxapple.ho.asm.bp.p_utg.gfa.fa"

    parser = argparse.ArgumentParser(description='Rescue scaffolds and generate sequences based on AGP and GFA files.')
    parser.add_argument('-e', '--edge_file', required=True, help='Input edge file')
    parser.add_argument('-a', '--agp_file', required=True, help='Input AGP file')
    parser.add_argument('-g', '--gfa_file', required=True, help='Input GFA file')
    parser.add_argument('-r', '--re_file', required=True, help='Input RE counts file')
    parser.add_argument('-f', '--fasta_file', required=True, help='Input fasta file of unitigs')

    args = parser.parse_args()

    edge_file = args.edge_file
    agp_file = args.agp_file
    gfa_file = args.gfa_file
    REFile = args.re_file
    fa_file = args.fasta_file

    G = read_gfa(gfa_file)
    ctg_RE_dict = read_RE(REFile)
    utg_rescue_dict = rescue(edge_file, agp_file, gfa_file, REFile)



    agp_df = read_agp_pd(agp_file)

    updated_df = update_agp_with_insert_lists(agp_df, utg_rescue_dict, ctg_RE_dict)

    scaffold_sequences_from_agp(updated_df, G, fa_file)