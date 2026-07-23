import argparse
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import igraph as ig



def read_RE(REFile):
    """read RE file"""
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
    """read gfa file"""
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
    """build igraph Graph"""
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
    """read agp file, return  {scaffold: [(utg, dir), ...]}"""
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
    scaffolded_utgs = {utg for utgs in scaffolds.values() for utg, _ in utgs}
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
                    if any(u in scaffolded_utgs for u in path_names[1:-1]):
                        continue

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
                    n_inserted = 0

                    for idx, _utg in enumerate(insert_path):

                        if _utg[0] not in utgs_set:
                            continue
                        n_inserted += 1
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
                    new_row['start'] = int(row['start']) + bp_set_off + (1 if n_inserted > 0 else 0)
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


def _reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def _gfa_connected(gfa_graph, a, b):
    """Return (connected, overlap_len) for placement a -> b using orientations."""
    utg, ori = a["utg"], a["ori"]
    nutg, nori = b["utg"], b["ori"]
    if ori == "+":
        ok = (nutg, nori) in gfa_graph[utg]["forward_list"]
    elif ori == "-":
        ok = (nutg, nori) in gfa_graph[utg]["reverse_list"]
    else:
        return False, 0
    if not ok:
        return False, 0
    ov = gfa_graph[utg]["edges"].get((nutg, nori), 0)
    return True, ov


def _interrupted_utgs(placements):
    """
    Utgs split into multiple different component ranges on this scaffold.

    Same (beg, end) repeated does not count as 打断; different break positions do.
    """
    ranges_by_utg = defaultdict(set)
    for p in placements:
        ranges_by_utg[p["utg"]].add((p["beg"], p["end"]))
    return {utg for utg, ranges in ranges_by_utg.items() if len(ranges) > 1}


def build_contigs_from_rescue_agp(
    updated_df,
    gfa_graph,
    fasta_file,
    gap_len=100,
    out_ctg2utg="gphase_final_ctg2utg.txt",
    out_contig_fa="gphase_final_contig.fasta",
    out_contig_agp="gphase_final_contig.agp",
):
    """
    Build contig fasta / ctg2utg / contig AGP from rescue AGP placements.

    Source of truth: rescue AGP (updated_df -> gphase_final.unitig.rescue.agp).

    - Each W placement is one instance; never clamp/reuse contig IDs.
    - Different (utg, comp_beg, comp_end) = different fragments.
    - Same path appearing again still gets a new contig ID (no dict overwrite).
    - Consecutive GFA-connected W placements merge into one contig walk.
    - If an utg is interrupted (split into different component ranges on this
      scaffold), stop merging: do not extend a walk through that utg.
    """
    utg_seq_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    placements_by_scaf = {}
    for scaffold, group in updated_df.groupby("scaffold", sort=False):
        group = group.reset_index(drop=True)
        placements = []
        for _, row in group.iterrows():
            if row["type"] != "W":
                continue
            utg_id = row["object"]
            beg = int(row["object_beg"])
            end = int(row["object_end"])
            ori = row["orientation"]
            if utg_id not in utg_seq_dict:
                raise KeyError(f"utg '{utg_id}' not found in fasta file")
            full = utg_seq_dict[utg_id]
            if beg < 1 or end < beg or end > len(full):
                raise ValueError(
                    f"error: start_pos={beg}, end_pos={end} for {utg_id} (len={len(full)})"
                )
            piece = full[beg - 1 : end]
            if ori == "-":
                piece = _reverse_complement(piece)
            placements.append(
                {
                    "utg": utg_id,
                    "beg": beg,
                    "end": end,
                    "ori": ori,
                    "seq": piece,
                    "token": f"{utg_id}:{beg}:{end}_{ori}",
                }
            )
        placements_by_scaf[scaffold] = placements

    contig_records = []
    ctg_idx = 0

    for scaffold, placements in placements_by_scaf.items():
        interrupted = _interrupted_utgs(placements)
        i = 0
        while i < len(placements):
            # Interrupted utg fragment: always emit alone, do not merge.
            if placements[i]["utg"] in interrupted:
                walk = [placements[i]]
                seq = placements[i]["seq"]
                j = i + 1
            else:
                walk = [placements[i]]
                seq = placements[i]["seq"]
                j = i + 1
                while j < len(placements):
                    # Stop before an interrupted utg; do not merge across it.
                    if placements[j]["utg"] in interrupted:
                        break
                    connected, ov = _gfa_connected(gfa_graph, walk[-1], placements[j])
                    if not connected:
                        break
                    nxt = placements[j]["seq"]
                    if ov > 0:
                        if ov >= len(nxt):
                            ov = 0
                        else:
                            nxt = nxt[ov:]
                    seq += nxt
                    walk.append(placements[j])
                    j += 1

            ctg_idx += 1
            ctg_id = f"ctg{ctg_idx:06d}l"
            path = "_".join(p["token"] for p in walk)
            contig_records.append(
                {
                    "ctg": ctg_id,
                    "path": path,
                    "seq": seq,
                    "scaffold": scaffold,
                }
            )
            for k in range(i, j):
                placements[k]["ctg"] = ctg_id
            i = j

    with open(out_ctg2utg, "w") as f2, open(out_contig_fa, "w") as fa:
        for rec in contig_records:
            f2.write(f"{rec['ctg']}\t{rec['path']}\n")
            fa.write(f">{rec['ctg']}\n")
            s = rec["seq"]
            for x in range(0, len(s), 80):
                fa.write(s[x : x + 80] + "\n")

    ctg_len = {rec["ctg"]: len(rec["seq"]) for rec in contig_records}

    agp_lines = []
    for scaffold, placements in placements_by_scaf.items():
        ordered_ctgs = []
        seen = set()
        for p in placements:
            c = p["ctg"]
            if c not in seen:
                ordered_ctgs.append(c)
                seen.add(c)

        part_number = 1
        current_pos = 1
        prev_ctg = None
        for ctg in ordered_ctgs:
            if prev_ctg is not None:
                agp_lines.append(
                    [
                        scaffold,
                        current_pos,
                        current_pos + gap_len - 1,
                        part_number,
                        "U",
                        gap_len,
                        "scaffold",
                        "yes",
                        "proximity_ligation",
                    ]
                )
                current_pos += gap_len
                part_number += 1

            length = ctg_len[ctg]
            agp_lines.append(
                [
                    scaffold,
                    current_pos,
                    current_pos + length - 1,
                    part_number,
                    "W",
                    ctg,
                    1,
                    length,
                    "+",
                ]
            )
            current_pos += length
            part_number += 1
            prev_ctg = ctg

    with open(out_contig_agp, "w") as f:
        for line in agp_lines:
            f.write("\t".join(map(str, line)) + "\n")

    return contig_records


def Rescue_base_graph(edge_file, agp_file, gfa_file, REFile, fa_file):
    gfa_graph, utgs_set = read_gfa(gfa_file)
    ctg_RE_dict = read_RE(REFile)
    utg_rescue_dict = rescue(edge_file, agp_file, gfa_file, REFile)
    agp_df = read_agp_pd(agp_file)
    updated_df = update_agp_with_insert_lists(agp_df, utg_rescue_dict, ctg_RE_dict, utgs_set)
    # updated_df is written as gphase_final_rescue.agp (-> gphase_final.unitig.rescue.agp)
    build_contigs_from_rescue_agp(
        updated_df,
        gfa_graph,
        fa_file,
        gap_len=100,
        out_ctg2utg="gphase_final_ctg2utg.txt",
        out_contig_fa="gphase_final_contig.fasta",
        out_contig_agp="gphase_final_contig.agp",
    )




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
