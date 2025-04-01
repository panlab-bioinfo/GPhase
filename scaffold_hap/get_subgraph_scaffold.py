from collections import defaultdict, deque
import csv
import networkx as nx
import argparse
import copy as cp
import matplotlib.pyplot as plt


def read_group(group_file):

    ctgs_list = list()
    with open(group_file, 'r') as file:
        for line in file:
            if line.startswith('u'):
                line = line.strip().split()
                ctgs_list.append(line[0])

    return ctgs_list

def read_digraph(digraph_file):

    digraph_dict = defaultdict()
    with open(digraph_file, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            if not line[0].startswith("s"):
                digraph_dict[tuple([line[0], line[1]])] = 1
    
    return digraph_dict


def read_l(l):
    hic_nei_dict = defaultdict(set)
    hic_links_dict = defaultdict()
    with open(l, 'r') as file:
        for line in file:
            if line.startswith("utg") or line.startswith("utig") :
                line = line.strip().split(',')
                hic_links_dict[tuple(sorted([line[0], line[1]]))] = float(line[2])
                hic_nei_dict[line[0]].add(line[1])
                hic_nei_dict[line[1]].add(line[0])
    return hic_links_dict, hic_nei_dict

def read_RE(REFile):
    ctg_RE_dict = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_dict[line[0]] = (line[1], line[2])
    return ctg_RE_dict

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

# 过滤子图： 如果子图中没有group中ctg则过滤
def filter_subgraph(subgraph_ctgs_dict, ctg_subgraph_dict, ctgs_list):

    filter_subgraph_ctgs_dict = defaultdict(list)
    filter_ctg_subgraph_dict = defaultdict()

    for subgraph_num, subgraph_ctgs in subgraph_ctgs_dict.items():
        intersection = set(subgraph_ctgs).intersection(set(ctgs_list))

        if len(intersection):
            filter_subgraph_ctgs_dict[subgraph_num] = subgraph_ctgs
        
        for ctg in subgraph_ctgs:
            filter_ctg_subgraph_dict[ctg] = subgraph_num
        
    
    return filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict

# 全部子图之间的连接
def get_subgraph_link(digraph_dict, subgraph_ctgs_dict, ctg_subgraph_dict):

    subgraph_connect_dict = defaultdict()
    for (ctg1, ctg2) in digraph_dict:

        subgraph1 = ctg_subgraph_dict.get(ctg1, None)
        subgraph2 = ctg_subgraph_dict.get(ctg2, None)

        if subgraph1 and subgraph2 and subgraph1 != subgraph2:
            subgraph_connect_dict[tuple([subgraph1, subgraph2])] = 1

    
    with open("subgraph_digraph.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['source', 'target'])

        for subgraph1, subgraph2 in subgraph_connect_dict:
            writer.writerow([subgraph1, subgraph2])
    
    return subgraph_connect_dict

# 找有向图N阶前驱
def n_hop_predecessors(G, node, N):
    """用 ancestors() + BFS 找 N 阶前驱"""
    ancestors = nx.ancestors(G, node)  # 所有前驱
    return {p for p in ancestors if nx.shortest_path_length(G, p, node) == N}

# 在包含group中ctg的子图上进行连接
def connect_subgraph(filter_subgraph_ctgs_dict, subgraph_connect_dict):

    subgraph_digraph = nx.DiGraph()

    # 添加边和权重
    for (node1, node2) in subgraph_connect_dict.keys():
        subgraph_digraph.add_edge(node1, node2)
    
    
    graph = nx.DiGraph()
    filter_subgraph_list = list(filter_subgraph_ctgs_dict.keys())

    for subgraph in filter_subgraph_list:
        graph.add_node(subgraph)


    # N_neighbors = 3
    # for node in graph.nodes():
    #     neighbors = list(n_hop_predecessors(subgraph_digraph, node, N_neighbors))
    #     neighbors_filter = [ node for node in neighbors if node in graph.nodes()]
    #     for nei_node in neighbors_filter:
    #         graph.add_edge(node, nei_node)

    # for node in graph.nodes():
    #     neighbors = list(n_hop_predecessors(subgraph_digraph, node, N_neighbors))
    #     neighbors_filter = [ node for node in neighbors if node in graph.nodes()]
    #     for nei_node in neighbors_filter:
    #         graph.add_edge(node, nei_node)


    Path_Length_T = 5

    for node1 in graph.nodes():
        for node2 in graph.nodes():
            if node1 != node2 :
                if node1 not in list(subgraph_digraph.nodes()) or node2 not in list(subgraph_digraph.nodes()):
                    continue
                # 检测当前graph是否存在路径，不存在则在subgraph_digraph中检测路径并得到路径中最后的子图与之链接
                if not nx.has_path(graph, node1, node2):
                    if nx.has_path(subgraph_digraph, node1, node2): 
                        
                        # paths = list(nx.all_simple_paths(subgraph_digraph, source=node1, target=node2))
                        # for path in paths:
                        path = nx.shortest_path(subgraph_digraph, source=node1, target=node2)
                        last_node_in_path = [ node for node in path if node in graph.nodes()][-2]
                        last_node_in_path_idx = list(path).index(last_node_in_path)
                        if len(path) - last_node_in_path_idx < Path_Length_T:
                            graph.add_edge(last_node_in_path, node2)

                else:
                    continue



    # 对每个连通分量进行check 并进行拓扑排序
    topo_order_list = list()
    for component in nx.weakly_connected_components(graph):
        break_points_set = set()

        if len(component) == 1:
            topo_order_list.append(component)
            continue

        subgraph_digraph = graph.subgraph(component).copy()
        # 执行拓扑排序
        while True:
            try:
                # 查找一个环
                cycle = nx.find_cycle(subgraph_digraph)
                # 删除环中的第一条边
                subgraph_digraph.remove_edge(*cycle[0])
            except nx.NetworkXNoCycle:
                # 没有环时退出循环
                break

        topo_order = list(nx.topological_sort(subgraph_digraph))
        print(f"topo_order:{topo_order}")

        # check : 相邻两点是否只存在一条路径, 如果没有路径或者存在多条路径则打断,并且无同源节点
        for i in range(len(topo_order)-1):

            # 两节点是否连通
            is_strongly_connected = nx.has_path(graph, topo_order[i], topo_order[i+1])
            if not is_strongly_connected:
                break_points_set.add(i+1)
                continue

            paths = list(nx.all_simple_paths(graph, topo_order[i], topo_order[i+1]))
            if len(paths) > 1:
                break_points_set.add(i+1)

            # 检测是否有相连的同源节点
            successors = list(graph.successors(topo_order[i]))
            predecessors = list(graph.predecessors(topo_order[i+1]))

            if len(successors) > 1 or len(predecessors) > 1:
                break_points_set.add(i+1)

        # 根据断点打断
        start = 0
        for break_points in break_points_set:
            topo_order_list.append(topo_order[start:break_points])
            start = break_points
        topo_order_list.append(topo_order[start:])


    # print(f"Final : {topo_order_list}")
    # check topo_order_list:
    topo_order_list_cp = cp.deepcopy(topo_order_list)
    for idx1, list1 in enumerate(topo_order_list):
        for idx2, list2 in enumerate(topo_order_list):
            if idx1 < idx2:
                overlap = set(list1) & set(list2)
                if len(overlap) > 0:
                    for subgraph in list(overlap):
                        topo_order_list_cp[idx1].remove(subgraph)
    
    print(f"Final : {topo_order_list_cp}")




    return topo_order_list_cp

def get_filter_subgraph_digraph(filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict, gfa_filePath, topo_order_list, digraph_dict):

    global_digraph = nx.DiGraph()

    # 添加边和权重
    for (node1, node2) in digraph_dict.keys():
        global_digraph.add_edge(node1, node2)
    
    for node in global_digraph.nodes():
        global_digraph.nodes[node]['direction'] = 1

    

    graphs_dict = dict()
    group_subgraph_dict, ctg_group_dict = defaultdict(list), defaultdict()
    subgraph_all_list = [ subgraph for subgraph_list in topo_order_list for subgraph in subgraph_list]
    for idx, subgraph_list in enumerate(topo_order_list):
        graph = nx.DiGraph()
        for subgraph in subgraph_list:
            group_subgraph_dict[f"group{idx}"].append(subgraph)
            for ctg in filter_subgraph_ctgs_dict[subgraph]:
                ctg_group_dict[ctg] = f"group{idx}"
                graph.add_node(ctg,direction=1)

        graphs_dict[f"group{idx}"] = graph

    # 构造有向图
    Path_Length_T = 5
    for graph in graphs_dict:
        nodes = list(graphs_dict[graph].nodes())
        for node1 in nodes:
            for node2 in nodes:
                if node1 != node2:
                    if node1 not in list(global_digraph.nodes()) or node2 not in list(global_digraph.nodes()):
                        continue
                    if not nx.has_path(graphs_dict[graph], node1, node2):
                        if nx.has_path(global_digraph, node1, node2):
                            path = nx.shortest_path(global_digraph, source=node1, target=node2)
                            last_node_in_path = [ node for node in path if node in nodes][-2]
                            last_node_in_path_idx = list(path).index(last_node_in_path)
                            if len(path) - last_node_in_path_idx < Path_Length_T:
                                graphs_dict[graph].add_edge(last_node_in_path, node2)

                    else:
                        continue


    edge_dict = defaultdict(tuple)
    with open(gfa_filePath, 'r') as file:
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
    
        
    return graphs_dict, global_digraph


def get_topological_sort(digraph, ctgs, hic_links_dict, hic_nei_dict, global_digraph):

    ctgs_dir_path_filter_dict = defaultdict(list)
    digraph_copy = cp.deepcopy(digraph)

    for u,v in digraph_copy.edges():
        if tuple(sorted([u,v])) in hic_links_dict and u in hic_nei_dict and v in hic_nei_dict:
            digraph_copy[u][v]['weight'] = hic_links_dict[tuple(sorted([u,v]))]
        else:
            digraph_copy[u][v]['weight'] = 0.0


    if len(ctgs) == 1:
        # 默认设置正向
        ctgs_dir_path_filter_dict[0] = [(ctgs[0],1)]
        return ctgs_dir_path_filter_dict


    # 执行拓扑排序
    while True:
        try:
            cycle = nx.find_cycle(digraph_copy)
            # print("找到环:", cycle)
            # topo_order = list(nx.topological_sort(digraph))
            edge_to_remove = min(cycle, key=lambda x: digraph_copy.edges[x]['weight'])
            digraph_copy.remove_edge(*edge_to_remove)

        except nx.NetworkXNoCycle:
            # print("所有环已处理完毕")
            break

    topo_order = list(nx.topological_sort(digraph_copy))

    ctgs_sort = [ ctg for ctg in topo_order if ctg in ctgs ]


    if len(ctgs_sort) != len(ctgs):
        ctgs_sort = list(set(graph.nodes()) & set(ctgs))


    
    
    # 检测相邻节点是否连通，若不连通则打断
    ctgs_sort_check = list()
    break_points_list = list()
    for i in range(len(ctgs_sort) - 1):
        # print(f"source:{ctgs_sort[i]}\ttarget:{ctgs_sort[i+1]}")
        try:
            if not nx.has_path(digraph, source=ctgs_sort[i], target=ctgs_sort[i+1]):
                break_points_list.append(i+1)
        except:
            break_points_list.append(i+1)


    # 根据断点打断
    start = 0
    for break_point in break_points_list:
        ctgs_sort_check.append(ctgs_sort[start:break_point])
        start = break_point
    ctgs_sort_check.append(ctgs_sort[start:])

    ctgs_sort_list = cp.deepcopy(ctgs_sort_check)


    for idx, ctgs_sort in enumerate(ctgs_sort_list):

        ctgs_dir_path = list()

        for i in range(len(ctgs_sort) - 1):
            shortest_path = nx.shortest_path(digraph, source=ctgs_sort[i], target=ctgs_sort[i+1])
            for node_j in shortest_path[:-1]:
                if node_j in global_digraph.nodes():
                    ctgs_dir_path.append((node_j, global_digraph.nodes[node_j]['direction']))
                else:
                    ctgs_dir_path.append((node_j, 1))
                
                
        if ctgs_sort[-1] in global_digraph.nodes():
            ctgs_dir_path.append((ctgs_sort[-1], global_digraph[ctgs_sort[-1]]))
        else:
            ctgs_dir_path.append((ctgs_sort[-1], 1))

    
        ctgs_dir_path_filter = [ (ctg, ctg_dir) for ctg, ctg_dir in ctgs_dir_path if ctg in ctgs_sort]
        

        ctgs_dir_path_filter_dict[str(idx)] = list(ctgs_dir_path_filter)
    # print(ctgs_dir_path_filter_dict)

    return ctgs_dir_path_filter_dict
    

def get_subgraph_group_inner_sort(graphs_dict, ctgs_list, ctg_RE_dict, global_digraph):

    filter_subgraph_sort_dir_dict = defaultdict(list)
    filter_subgraph_sort_dict = defaultdict(list)
    ctg_subgraphList_dict = defaultdict()
    subgraphGroup_RE_dict = defaultdict(tuple)
    
    subgraphGroup_idx = 1
    for graph in graphs_dict.values():
        # print(','.join(list(graph.nodes())))

        if len(graph.nodes()) == 0:
            continue

        ctgs = list(set(graph.nodes()) & set(ctgs_list))
        if len(ctgs) == 1:
            filter_subgraph_sort_dir_dict[subgraphGroup_idx] = [(ctgs[0], 1)]
            filter_subgraph_sort_dict[subgraphGroup_idx] = [ctgs[0]]
            subgraphGroup_idx += 1
            continue



        ctgs_dir_path_filter_dict = get_topological_sort(graph, ctgs, hic_links_dict, hic_nei_dict, global_digraph)
        
        
        for subgraphGroup, ctgs_dir_path_filter in ctgs_dir_path_filter_dict.items():
            ctgs_path_filter_list = [ pair[0] for pair in ctgs_dir_path_filter]
            filter_subgraph_sort_dir_dict[subgraphGroup_idx] = list(ctgs_dir_path_filter)
            filter_subgraph_sort_dict[subgraphGroup_idx] = list(ctgs_path_filter_list)
            subgraphGroup_idx += 1



    # with open(f"group.subgraph.sort.dir.cluster.txt", 'w') as file: 
    #     for subgraph, ctgs_dir_list in filter_subgraph_sort_dir_dict.items():
    #         ctgs_dir_trans_list = [ str(ctg)  + str('+' if str(ctg_dir)=='1' else '-') for ctg, ctg_dir in ctgs_dir_list]
    #         file.write(f"{subgraph}\t{len(ctgs_dir_trans_list)}\t{' '.join(ctgs_dir_trans_list)}\n")
    
    # with open(f"group.subgraph.sort.cluster.txt", 'w') as file: 
    #     for subgraph, ctgs_list in filter_subgraph_sort_dict.items():
    #         file.write(f"{subgraph}\t{len(ctgs_dir_list)}\t{' '.join(ctgs_list)}\n")

    # with open(f"subgraphGroup.counts_GATC.txt", 'w') as file: 
    #     for subgraph, pair in subgraphGroup_RE_dict.items():
    #         file.write(f"{subgraph}\t{pair[0]}\t{pair[1]}\n")

    return filter_subgraph_sort_dir_dict, filter_subgraph_sort_dict


def get_AGP(filter_subgraph_sort_dir_dict, ctg_RE_dict):

    with open("subgraphGroup.agp", 'w') as file:
        for Group_idx,(subgraphGroup, ctgs_dir) in enumerate(filter_subgraph_sort_dir_dict.items(), 1):
            len_sum = 0
            for idx, (ctg, ctg_dir) in enumerate(ctgs_dir):
                ctg_len = int(ctg_RE_dict[ctg][1])
                file.write(f"scaffold_{Group_idx}\t{len_sum+1}\t{ctg_len+len_sum}\t{(idx+1)*2-1}\t{'W'}\t{ctg}\t{1}\t{ctg_RE_dict[ctg][1]}\t{str('+' if str(ctg_dir)=='1' else '-')}\n")
                if idx != len(ctgs_dir)-1:
                    file.write(f"scaffold_{Group_idx}\t{ctg_len+len_sum+1}\t{ctg_len+len_sum+100}\t{(idx+1)*2}\t{'U'}\t{100}\t{'scaffold'}\t{'yes'}\t{'proximity_ligation'}\n")
                    len_sum += ctg_len + 100




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Contigs are initially mounted based on subgraphs and information between subgraphs")
    parser.add_argument('-c', '--group_file', required=True,
                        help='<filepath> group file for ctgs or utgs')
    parser.add_argument('-l', '--hic_link', required=True,
                        help='<filepath> hic links for ctgs or utgs')
    parser.add_argument('-subgraph', '--subgraph_file', required=True,
                        help='<filepath>all subgraph for GFA')
    parser.add_argument('-graph', '--graph', required=True,
                        help='<filepath>graph file for GFA')
    parser.add_argument('-digraph', '--digraph_file', required=True,
                        help='<filepath> digraph for GFA')
    parser.add_argument('-r', '--REFile', required=True,
                        help='<filepath>Length and REs of contig or unitig')

    # group_file = "group4.txt"
    # subgraph_file = "group_ctgs_All.txt"
    # gfa_filePath = "waxapple.asm.bp.p_utg.gfa"
    # subgraph_links_file = "subgraph_digraph.csv"
    # l = "group4.links.csv"
    # REFile = "alignment.counts_GATC.txt"

    args = parser.parse_args()
    group_file = args.group_file
    subgraph_file = args.subgraph_file
    gfa_filePath = args.graph
    REFile = args.REFile
    digraph_file = args.digraph_file
    hic_link_file = args.hic_link

    


    digraph_dict = read_digraph(digraph_file)
    ctg_RE_dict = read_RE(REFile)
    hic_links_dict, hic_nei_dict = read_l(hic_link_file)

    ctgs_list = read_group(group_file)

    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)
    filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict = filter_subgraph(subgraph_ctgs_dict, ctg_subgraph_dict, ctgs_list)

    subgraph_connect_dict = get_subgraph_link(digraph_dict, subgraph_ctgs_dict, ctg_subgraph_dict)

    topo_order_list = connect_subgraph(filter_subgraph_ctgs_dict, subgraph_connect_dict)
    graphs_dict, global_digraph = get_filter_subgraph_digraph(filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict, gfa_filePath, topo_order_list, digraph_dict)
    filter_subgraph_sort_dir_dict, filter_subgraph_sort_dict = get_subgraph_group_inner_sort(graphs_dict, ctgs_list, ctg_RE_dict, global_digraph)

    get_AGP(filter_subgraph_sort_dir_dict, ctg_RE_dict)


