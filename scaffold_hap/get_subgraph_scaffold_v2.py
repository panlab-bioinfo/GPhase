from collections import defaultdict, deque
import csv
import networkx as nx
import argparse
import copy as cp
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism



def read_group(group_file):

    ctgs_set = set()
    with open(group_file, 'r') as file:
        for line in file:
            if line.startswith('u'):
                line = line.strip().split()
                ctgs_set.add(line[0])

    return list(ctgs_set)

def read_digraph(digraph_file):

    digraph_dict = defaultdict()
    with open(digraph_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if not line[0].startswith("s"):
                digraph_dict[tuple([line[0], line[1]])] = 1

    digraph = nx.DiGraph()

    # 添加边和权重
    for (node1, node2) in digraph_dict.keys():
        digraph.add_edge(node1, node2)
    
    return digraph, digraph_dict

# 读取子图
def read_subgraph(subgraph_file):

    subgraph_ctgs_dict = defaultdict(list)
    ctg_subgraph_dict = defaultdict()
    with open(subgraph_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            subgraph_num = line[0]
            ctgs = list(line[2].split(','))
            subgraph_ctgs_dict[subgraph_num] = ctgs
            for ctg in ctgs:
                ctg_subgraph_dict[ctg] = subgraph_num

    return subgraph_ctgs_dict, ctg_subgraph_dict

def read_RE(RE_file):

    ctg_RE_dict = defaultdict(tuple)
    with open(RE_file, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (line[1], line[2])
    return ctg_RE_dict

# 全部子图之间的连接
def get_subgraph_link(digraph_dict, subgraph_ctgs_dict, ctg_subgraph_dict):

    subgraph_connect_dict = defaultdict()
    for (ctg1, ctg2) in digraph_dict:

        subgraph1 = ctg_subgraph_dict.get(ctg1, None)
        subgraph2 = ctg_subgraph_dict.get(ctg2, None)

        if subgraph1 and subgraph2 and subgraph1 != subgraph2:
            subgraph_connect_dict[tuple([subgraph1, subgraph2])] = 1

    digraph = nx.DiGraph()

    # 添加边和权重
    for (node1, node2) in subgraph_connect_dict.keys():
        digraph.add_edge(node1, node2)

    
    return digraph, subgraph_connect_dict

def break_cycles_keep_low_weight(G, weight_attr='weight'):
    G = G.copy()
    while True:
        try:
            # 找到一个环，返回节点序列
            cycle = next(nx.simple_cycles(G))
            # 在这个环中找出权重最大的边
            max_weight = float('-inf')
            edge_to_remove = None
            for i in range(len(cycle)):
                u = cycle[i]
                v = cycle[(i + 1) % len(cycle)]
                w = G[u][v].get(weight_attr, 1) 
                if w > max_weight:
                    max_weight = w
                    edge_to_remove = (u, v)
            # 删除这条权重最大的边
            if edge_to_remove:
                G.remove_edge(*edge_to_remove)
        except StopIteration:
            break
    return G

def get_nth_order_successors(G, start_node, N, return_by_level=False, reverse = False):
    if reverse:
        reversed_G = G.reverse()
        path_lengths = nx.single_source_shortest_path_length(reversed_G, start_node, cutoff=N)
    else:
        path_lengths = nx.single_source_shortest_path_length(G, start_node, cutoff=N)

    if return_by_level:
        level_dict = defaultdict(set)
        for node, length in path_lengths.items():
            if length > 0:
                level_dict[length].add(node)
        return dict(level_dict)
    else:
        return {node for node, length in path_lengths.items() if 0 < length <= N}


def iter_node_nth_successors(G, subgraph_digraph, ctgs_set, subgraph_ctgs_dict, ctg_subgraph_dict, N):

    first_digraph = nx.DiGraph()
    for node in ctgs_set:
        first_digraph.add_node(node, weight=0)
    
    for node in ctgs_set:

        if node not in G.nodes():
            continue

        level_dict = get_nth_order_successors(G, node, N, return_by_level=True)

        for level in level_dict:
            for node_nei in [ node for node in level_dict[level] if node in ctgs_set]:
                first_digraph.add_edge(node, node_nei, weight=int(level))
                


    nx.write_edgelist(first_digraph, path='edges.csv', delimiter=',',data=True)

    # 去除虚假边
    # 检测两个点之间最短路径，去除三个以上节点
    first_digraph_cp = cp.deepcopy(first_digraph)
    for u, v in first_digraph_cp.edges():
        try:
            # 找路径使用subgraph链接图
            subgraph1 = ctg_subgraph_dict.get(u, None)
            subgraph2 = ctg_subgraph_dict.get(v, None)
            if not subgraph1 or not subgraph2:
                continue
            paths = nx.all_simple_paths(subgraph_digraph, subgraph1, subgraph2, cutoff=20)
            for path in paths:

                filter_path = [ subgraph for subgraph in path if set(ctgs_set) & set(subgraph_ctgs_dict[subgraph]) ]
                if len(filter_path) > 2:
                    first_digraph.remove_edge(u, v)
                    break
        except:
            pass

    # 去环
    first_digraph_DAG = break_cycles_keep_low_weight(first_digraph)
    nx.write_edgelist(first_digraph_DAG, path='edges_DAG.csv', delimiter=',',data=True)



    first_digraph_reduce = nx.transitive_reduction(first_digraph_DAG)

    for u, v, data in first_digraph_DAG.edges(data=True):
        weight = data.get('weight', None)
        if first_digraph_reduce.has_edge(u, v):
            first_digraph_reduce[u][v]['weight'] = weight

    nx.write_edgelist(first_digraph_reduce, path='edges_reduce.csv', delimiter=',',data=True)

    return first_digraph_reduce

def Nth_nei(G, reverse = False):

    # N-th nei
    while True:

        if reverse:
            nodes_with_many = [n for n in G.nodes if len(list(G.predecessors(n))) > 1]
            print(f"N-th nei: {nodes_with_many}")
        else:
            nodes_with_many = [n for n in G.nodes if len(list(G.successors(n))) > 1]
            print(f"N-th nei: {nodes_with_many}")

        if len(nodes_with_many) == 0:
            break
        for node in nodes_with_many:
            node_nei_dict = defaultdict(list)
            node_list = list(G.predecessors(node)) if reverse else list(G.successors(node))

            if reverse:
                for branch_node in node_list:
                    node_nei_dict[branch_node] = get_nth_order_successors(G, branch_node, 5, return_by_level=False)
            else:
                for branch_node in node_list:
                    node_nei_dict[branch_node] = get_nth_order_successors(G, branch_node, 5, return_by_level=False, reverse=True)

            # 找 分支节点N阶邻居之间的交集
            intersect_list = node_nei_dict[node_list[0]]
            for branch_node, nei_1 in node_nei_dict.items():
                intersect_list_cp = cp.deepcopy(intersect_list)
                intersect_list = set(intersect_list_cp) & set(nei_1)
            
            # 去除节点与分支节点之间的连接
            # 若 分支中存在 tip节点 则只去除tip节点
            if node_list == 2 and (len(node_nei_dict[node_list[0]]) == 0 )  !=  (len(node_nei_dict[node_list[1]]) == 0):
                if len(node_nei_dict[node_list[0]]) == 0:
                    if G.has_edge(node, node_list[0]):
                        G.remove_edge(node, node_list[0])
                    if G.has_edge(node_list[0], node):
                        G.remove_edge(node_list[0], node)
                else:
                    if G.has_edge(node, node_list[0]):
                        G.remove_edge(node, node_list[1])
                    if G.has_edge(node_list[0], node):
                        G.remove_edge(node_list[1], node)
            else:
                for branch_node in node_list:
                    if G.has_edge(node, branch_node):
                        G.remove_edge(node, branch_node)
                    if G.has_edge(branch_node, node):
                        G.remove_edge(branch_node, node)

            # 交集存在
            if intersect_list and not reverse:
                for intersect_node in intersect_list:
                    predecessors = list(G.predecessors(intersect_node))
                    successors = list(G.successors(intersect_node))

                    if len(predecessors) > 1 and len(successors) == 1:
                        for predecessor in predecessors:
                            if G.has_edge(predecessor, intersect_node):
                                G.remove_edge(predecessor, intersect_node)
                        if nx.has_path(G, node, intersect_node):
                            G.add_edge(node, intersect_node)
                        break   
    return G


def simple_scaffold_Graph(G):

    # 模式P1
    P1 = nx.DiGraph()
    P1.add_edges_from([
        (1, 2), (1, 4), (3, 4)
    ])

    matcher = isomorphism.DiGraphMatcher(G, P1)
    if matcher.subgraph_is_isomorphic():
        print("找到P1同构子图")
        for subgraph_mapping in matcher.subgraph_isomorphisms_iter():
            start_node = [ node for node,num in subgraph_mapping.items() if num== 1][0]
            end_node = [ node for node,num in subgraph_mapping.items() if num== 4][0]
            if G.has_edge(start_node, end_node):
                G.remove_edge(start_node, end_node)
            if G.has_edge(end_node, start_node):
                G.remove_edge(end_node, start_node)
            print("匹配P1子图映射：", subgraph_mapping)
    else:
        print("未找到P1同构子图")

    P2 = nx.DiGraph()
    P2.add_edges_from([
        (1, 2), (1, 3), (2, 4), (3,4)
    ])

    while True:
        matcher = isomorphism.DiGraphMatcher(G, P2)
        if matcher.subgraph_is_isomorphic():
            print("找到P2同构子图")
            for subgraph_mapping in matcher.subgraph_isomorphisms_iter():
                node_1 = [ node for node,num in subgraph_mapping.items() if num== 1][0]
                node_2 = [ node for node,num in subgraph_mapping.items() if num== 2][0]
                node_3 = [ node for node,num in subgraph_mapping.items() if num== 3][0]
                node_4 = [ node for node,num in subgraph_mapping.items() if num== 4][0]
                if G.has_edge(node_1, node_2):
                    G.remove_edge(node_1, node_2)
                if G.has_edge(node_1, node_3):
                    G.remove_edge(node_1, node_3)
                if G.has_edge(node_2, node_4):
                    G.remove_edge(node_2, node_4)
                if G.has_edge(node_3, node_4):
                    G.remove_edge(node_3, node_4)
                if nx.has_path(G, node_1, node_4):
                    G.add_edge(node_1, node_4)
                print("匹配子图P2映射：", subgraph_mapping)
        else:
            print("未找到同构子图")
            break

    return G

def get_linear_path(G, start_node):

    path = [start_node]
    current = start_node
    while True:
        successors = list(G.successors(current))
        if not successors:
            break 
        current = successors[0]
        path.append(current)
    return path




def orient(G, digraph, GFA_file):
    # G: scaffold 图
    # digraph ：全部的有向边
    # global_digraph ： 边节点朝向的有向图

    scaffold_list = list()
    global_digraph = nx.DiGraph()

    for (node1, node2) in digraph.edges():
        global_digraph.add_edge(node1, node2)
    
    for node in global_digraph.nodes():
        global_digraph.nodes[node]['direction'] = 1

    edge_dict = defaultdict(tuple)
    with open(GFA_file, 'r') as file:
        for line in file:
            line = line.strip().split()
            if(line[0] == "L"):
                dir_node1 = 1 if line[2]=="+" else 0
                dir_node2 = 1 if line[4]=="+" else 0
                edge_dict[tuple([line[1], line[3]])] = (dir_node1, dir_node2)

    for u, v in global_digraph.edges():
        dir_tuple = edge_dict[tuple([u,v])]
        if dir_tuple:
            global_digraph.nodes[u]['direction'] = dir_tuple[0]
            global_digraph.nodes[v]['direction'] = dir_tuple[1]

    for idx, component in enumerate(nx.weakly_connected_components(G)):

        # 获取每个连通分量的起始节点
        start_node_list = [ node for node in list(component) if len(list(G.predecessors(node))) == 0]
        if len(start_node_list) > 1:
            print(f"连通分量{idx} 中起始节点个数大于1")
            return False
        elif start_node_list:
            start_node = start_node_list[0]
            print(f"连通分量{idx}, start_node : {start_node}")
        else:
            print(f"遍历连通分量{idx} ERROR...")

        # 得到线性路径;
        line_path_list = get_linear_path(G, start_node)
        line_path_All_list, line_path_All_filtered_list= list(), list()

        for i in range(len(line_path_list) - 1):
            shortest_path = nx.shortest_path(global_digraph, source=line_path_list[i], target=line_path_list[i+1])
            for node in shortest_path[:-1]:
                if node in global_digraph.nodes():
                    line_path_All_list.append((node, global_digraph.nodes[node]['direction']))
                else:
                    line_path_All_list.append((node, 1))

        if line_path_list[-1] in global_digraph.nodes():
            line_path_All_list.append((line_path_list[-1], global_digraph.nodes[line_path_list[-1]]['direction']))
        else:
            line_path_All_list.append((line_path_list[-1], 1))

        line_path_All_filtered_list = [ (ctg, ctg_dir) for ctg, ctg_dir in line_path_All_list if ctg in line_path_list]

        scaffold_list.append(line_path_All_filtered_list)

    return  scaffold_list


def get_AGP(scaffold_list, ctg_RE_dict):

    with open("subgraphGroup.agp", 'w') as file:
        for Group_idx, scaffold_utgs_list in enumerate(scaffold_list, 1):
            len_sum = 0
            for idx, (ctg, ctg_dir) in enumerate(scaffold_utgs_list):
                ctg_len = int(ctg_RE_dict[ctg][1])
                file.write(f"scaffold_{Group_idx}\t{len_sum+1}\t{ctg_len+len_sum}\t{(idx+1)*2-1}\t{'W'}\t{ctg}\t{1}\t{ctg_RE_dict[ctg][1]}\t{str('+' if str(ctg_dir)=='1' else '-')}\n")
                if idx != len(scaffold_utgs_list)-1:
                    file.write(f"scaffold_{Group_idx}\t{ctg_len+len_sum+1}\t{ctg_len+len_sum+100}\t{(idx+1)*2}\t{'U'}\t{100}\t{'scaffold'}\t{'yes'}\t{'proximity_ligation'}\n")
                    len_sum += ctg_len + 100

def Get_subgraph_scaffold_v2(gfa_file, RE_file, digraph_file, group_file, subgraph_file):

    # gfa_file = "C88.h.asm.bp.p_utg.gfa"
    # RE_file = "C88.RE_counts.txt"
    # digraph_file = "C88.digraph.csv"
    # group_file = "group2.txt"
    # subgraph_file = "group_ctgs_All.txt"

    # args = parser.parse_args()
    # gfa_file = args.graph_file
    # RE_file = args.RE_file
    # digraph_file = args.digraph_file
    # group_file = args.group_file
    # subgraph_file = args.subgraph_file

    ctgs_set = read_group(group_file)
    ctg_RE_dict = read_RE(RE_file)
    digraph, digraph_dict = read_digraph(digraph_file)
    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)

    subgraph_digraph, subgraph_connect_dict = get_subgraph_link(digraph_dict, subgraph_ctgs_dict, ctg_subgraph_dict)

    reduce_DAG = iter_node_nth_successors(digraph, subgraph_digraph,ctgs_set, subgraph_ctgs_dict, ctg_subgraph_dict, 15)

    reduce_DAG_clear = simple_scaffold_Graph(reduce_DAG)
    G1 = Nth_nei(reduce_DAG_clear, reverse = False)
    G2 = Nth_nei(G1, reverse = True)

    # 去环
    G_DAG = break_cycles_keep_low_weight(G2)
    for u in G_DAG.nodes():
        for v in G_DAG.nodes():
            if u != v and G_DAG.has_edge(u,v) and G_DAG.has_edge(v,u):
                G_DAG.remove_edge(u,v)

    nx.write_edgelist(G_DAG, path='edges_clear.csv', delimiter=',',data=True)
    scaffold_list = orient(G_DAG, digraph, gfa_file)
    get_AGP(scaffold_list, ctg_RE_dict)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Contigs are initially mounted based on subgraphs and information between subgraphs")
    parser.add_argument('-c', '--group_file', required=True,
                        help='<filepath> group file for ctgs or utgs')
    parser.add_argument('-subgraph', '--subgraph_file', required=True,
                        help='<filepath>all subgraph for GFA')
    parser.add_argument('-graph', '--graph_file', required=True,
                        help='<filepath>graph file for GFA')
    parser.add_argument('-digraph', '--digraph_file', required=True,
                        help='<filepath> digraph for GFA')
    parser.add_argument('-r', '--RE_file', required=True,
                        help='<filepath>Length and REs of contig or unitig')


    # gfa_file = "C88.h.asm.bp.p_utg.gfa"
    # RE_file = "C88.RE_counts.txt"
    # digraph_file = "C88.digraph.csv"
    # group_file = "group2.txt"
    # subgraph_file = "group_ctgs_All.txt"

    args = parser.parse_args()
    gfa_file = args.graph_file
    RE_file = args.RE_file
    digraph_file = args.digraph_file
    group_file = args.group_file
    subgraph_file = args.subgraph_file

    Get_subgraph_scaffold_v2(gfa_file, RE_file, digraph_file, group_file, subgraph_file)


