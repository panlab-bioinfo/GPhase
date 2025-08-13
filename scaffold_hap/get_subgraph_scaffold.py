from collections import defaultdict, deque
import csv
import networkx as nx
import argparse
import copy as cp
import matplotlib.pyplot as plt
import networkit as nk


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
    
    return digraph_dict


def read_l(l):
    hic_nei_dict = defaultdict(set)
    hic_links_dict = {}

    with open(l, newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            if not row:
                continue
            if row[0].startswith(('utg', 'utig')):
                a, b, w = row[0], row[1], float(row[2])

                # 用 tuple 排序缓存减少运算
                key = (a, b) if a < b else (b, a)

                hic_links_dict[key] = w
                hic_nei_dict[a].add(b)
                hic_nei_dict[b].add(a)

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


def nx_to_nk_graph(nx_graph, weighted=False):

    # 1. 节点映射：nx节点名 ➜ index
    node_list = list(nx_graph.nodes())
    node_map = {node: idx for idx, node in enumerate(node_list)}
    rev_map = {idx: node for node, idx in node_map.items()}

    # 2. 初始化 NetworKit 图
    nk_graph = nk.graph.Graph(n=len(node_list), weighted=weighted, directed=nx_graph.is_directed())


    for u, v, data, in nx_graph.edges(data=True):
        if weighted:
            weight = data.get("weight", 1.0)
            nk_graph.addEdge(node_map[u], node_map[v], weight)
        else:
            nk_graph.addEdge(node_map[u], node_map[v])

    return nk_graph, node_map, rev_map


def has_multiple_paths(graph, source, target, limit=2):
    count = 0

    def dfs(node, visited):
        nonlocal count
        if count >= limit:
            return
        if node == target:
            count += 1
            return
        for neighbor in graph.successors(node):
            if neighbor not in visited:
                dfs(neighbor, visited | {neighbor})

    dfs(source, {source})
    return count >= limit


# 找有向图N阶前驱
def n_hop_predecessors(G, node, N):
    n_hop_predecessors_set = set()
    ancestors = nx.ancestors(G, node) 

    G_nk, node_map, rev_map = nx_to_nk_graph(G, weighted=True)

    for p in ancestors:
        source_node = p
        target_node = node
        source_idx = node_map[source_node]
        target_idx = node_map[target_node]

        dijkstra = nk.distance.Dijkstra(G_nk, source_idx, storePaths=True)
        dijkstra.run()

        distance = dijkstra.distance(target_idx)
        if int(distance) == N:
            n_hop_predecessors_set.add(p)
        
    return n_hop_predecessors_set

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

    graph_nodes_list = list(graph.nodes())
    subgraph_digraph_nk, node_map, rev_map = nx_to_nk_graph(subgraph_digraph, weighted=False)

    Path_Length_T = 5

    subgraph_nodes = set(subgraph_digraph.nodes())

    for node1 in graph_nodes_list:
        if node1 not in subgraph_nodes:
            continue

        source_idx = node_map[node1]
        dijkstra = nk.distance.Dijkstra(subgraph_digraph_nk, source_idx, storePaths=True)
        dijkstra.run()

        for node2 in graph_nodes_list:
            if node1 == node2 or node2 not in subgraph_nodes:
                continue

            # graph 无路径，但 subgraph_digraph 中存在路径
            if nx.has_path(graph, node1, node2):
                continue

            target_idx = node_map[node2]
            distance = dijkstra.distance(target_idx)

            if distance == float('inf'):
                continue

            path_indices = dijkstra.getPath(target_idx)
            path = [rev_map[i] for i in path_indices]

            # 找到 path 中倒数第二个出现在 graph 中的节点
            for i in range(len(path) - 2, -1, -1):
                if path[i] in graph_nodes_list:
                    last_node_in_path = path[i]
                    last_node_in_path_idx = i
                    break
            else:
                continue  # 没找到合适的连接点

            # 如果路径后段长度小于阈值，就添加边
            if len(path) - last_node_in_path_idx < Path_Length_T:
                graph.add_edge(last_node_in_path, node2)



    # 对每个连通分量进行check 并进行拓扑排序
    topo_order_list = list()
    for component in nx.weakly_connected_components(graph):
        break_points_set = set()

        if len(component) == 1:
            topo_order_list.append(component)
            continue

        subgraph_digraph = graph.subgraph(component).copy()


        while True:
            subgraph_digraph_nk, node_map, rev_map = nx_to_nk_graph(subgraph_digraph, weighted=False)
            scc = nk.components.StronglyConnectedComponents(subgraph_digraph_nk)
            scc.run()
            components = scc.getComponents()

            has_cycle = False

            for comp in components:
                if len(comp) > 1:
                    has_cycle = True
                    nodes = [rev_map[i] for i in comp]

                    for u in nodes:
                        for v in subgraph_digraph.successors(u):
                            if v in nodes:
                                subgraph_digraph.remove_edge(u, v)
                                break
                        else:
                            continue
                        break
                    break 
            if not has_cycle:
                break



        topo_order = list(nx.topological_sort(subgraph_digraph))
        # print(f"topo_order:{topo_order}")

        # check : 相邻两点是否只存在一条路径, 如果没有路径或者存在多条路径则打断,并且无同源节点
        for i in range(len(topo_order)-1):

            # 两节点是否连通
            is_strongly_connected = nx.has_path(graph, topo_order[i], topo_order[i+1])
            if not is_strongly_connected:
                break_points_set.add(i+1)
                continue

            # path number > 1
            if has_multiple_paths(graph, topo_order[i], topo_order[i+1]):
                break_points_set.add(i+1)

            #检测是否有相连的同源节点
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


    topo_order_list_cp = cp.deepcopy(topo_order_list)
    for idx1, list1 in enumerate(topo_order_list):
        for idx2, list2 in enumerate(topo_order_list):
            if idx1 < idx2:
                overlap = set(list1) & set(list2)
                if len(overlap) > 0:
                    for subgraph in list(overlap):
                        if subgraph in topo_order_list_cp[idx1]:
                            topo_order_list_cp[idx1].remove(subgraph)
    
    # print(f"Final : {topo_order_list_cp}")
    return topo_order_list_cp




def get_filter_subgraph_digraph(filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict, gfa_filePath, topo_order_list, digraph_dict):

    global_digraph = nx.DiGraph()
    

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

    G_nk, node_map, rev_map = nx_to_nk_graph(global_digraph, weighted=True)


    Path_Length_T = 5

    for graph_name, subgraph in graphs_dict.items():
        nodes = list(subgraph.nodes())
        node_set = set(nodes)
        reachable_dict = {n: nx.descendants(subgraph, n) for n in subgraph.nodes()}

        for i, node1 in enumerate(nodes):
            if node1 not in global_digraph:
                continue

            source_idx = node_map.get(node1)
            if source_idx is None:
                continue

            dijkstra = nk.distance.Dijkstra(G_nk, source_idx, storePaths=True)
            dijkstra.run()

            for j, node2 in enumerate(nodes):
                if i == j or node2 not in global_digraph:
                    continue

                # if nx.has_path(subgraph, node1, node2):
                #     continue 

                if node2 in reachable_dict[node1]:
                    continue

                target_idx = node_map.get(node2)
                if target_idx is None:
                    continue

                distance = dijkstra.distance(target_idx)
                if distance == float('inf'):
                    continue

                path_indices = dijkstra.getPath(target_idx)
                path = [rev_map[i] for i in path_indices]

                # 找 path 中倒数第二个在子图节点集的点
                internal_path = [n for n in path if n in node_set]
                if len(internal_path) < 2:
                    continue  # 无法形成合理的边

                last_node_in_path = internal_path[-2]
                last_node_in_path_idx = path.index(last_node_in_path)

                if len(path) - last_node_in_path_idx < Path_Length_T:
                    subgraph.add_edge(last_node_in_path, node2)


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

    digraph_nk, node_map, rev_map = nx_to_nk_graph(digraph, weighted=True)

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


    while True:
            digraph_copy_nk, node_map, rev_map = nx_to_nk_graph(digraph_copy, weighted=False)
            scc = nk.components.StronglyConnectedComponents(digraph_copy_nk)
            scc.run()
            components = scc.getComponents()

            has_cycle = False

            for comp in components:
                if len(comp) > 1:
                    has_cycle = True
                    nodes = [rev_map[i] for i in comp]

                    for u in nodes:
                        for v in digraph_copy.successors(u):
                            if v in nodes:
                                digraph_copy.remove_edge(u, v)
                                break
                        else:
                            continue
                        break
                    break 

            if not has_cycle:
                break

    topo_order = list(nx.topological_sort(digraph_copy))

    ctgs_sort = [ ctg for ctg in topo_order if ctg in ctgs ]


    if len(ctgs_sort) != len(ctgs):
        ctgs_sort = list(set(graph.nodes()) & set(ctgs))
    
    # 检测相邻节点是否连通，若不连通则打断;如存在同源节点则断开
    ctgs_sort_check = list()
    break_points_list = list()
    for i in range(len(ctgs_sort) - 1):
        try:
            if not nx.has_path(digraph, source=ctgs_sort[i], target=ctgs_sort[i+1]):
                break_points_list.append(i+1)
            
            successors = list(digraph.successors(ctgs_sort[i]))
            predecessors = list(digraph.predecessors(ctgs_sort[i+1]))

            if len(successors) > 1 or len(predecessors) > 1:
                break_points_set.add(i+1)


        except:
            break_points_list.append(i+1)


    # 根据断点打断
    start = 0
    for break_point in break_points_list:
        ctgs_sort_check.append(ctgs_sort[start:break_point])
        start = break_point
    ctgs_sort_check.append(ctgs_sort[start:])

    ctgs_sort_list = cp.deepcopy(ctgs_sort_check)

    path_cache = {}

    for idx, ctgs_sort in enumerate(ctgs_sort_list):

        ctgs_dir_path = []

        for i in range(len(ctgs_sort) - 1):
            source_node = ctgs_sort[i]
            target_node = ctgs_sort[i + 1]
            source_idx = node_map[source_node]
            target_idx = node_map[target_node]

            path_key = (source_idx, target_idx)

            # 检查缓存中是否已有路径
            if path_key in path_cache:
                path_indices = path_cache[path_key]
            else:
                dijkstra = nk.distance.Dijkstra(digraph_nk, source_idx, storePaths=True)
                dijkstra.run()
                 # 将所有可达路径存入 cache
                for possible_target_idx in range(digraph_nk.numberOfNodes()):
                    
                    if dijkstra.distance(possible_target_idx) < float('inf'):
                        path_cache[(source_idx, possible_target_idx)] = dijkstra.getPath(possible_target_idx)

                path_indices = path_cache.get(path_key, [])

            shortest_path = [rev_map[i] for i in path_indices]

            for node_j in shortest_path[:-1]:
                if node_j in global_digraph.nodes():
                    ctgs_dir_path.append((node_j, global_digraph.nodes[node_j].get('direction', 1)))
                else:
                    ctgs_dir_path.append((node_j, 1))

        # 处理末尾节点方向
        if ctgs_sort:
            last_ctg = ctgs_sort[-1]
            if last_ctg in global_digraph.nodes():
                ctgs_dir_path.append((last_ctg, global_digraph.nodes[last_ctg].get('direction', 1)))
            else:
                ctgs_dir_path.append((last_ctg, 1))

            ctgs_dir_path_filter = [(ctg, dir_) for ctg, dir_ in ctgs_dir_path if ctg in ctgs_sort]

            ctgs_dir_path_filter_dict[str(idx)] = ctgs_dir_path_filter


    return ctgs_dir_path_filter_dict
    

def get_subgraph_group_inner_sort(graphs_dict, ctgs_list, ctg_RE_dict, global_digraph, hic_links_dict, hic_nei_dict):

    filter_subgraph_sort_dir_dict = defaultdict(list)
    filter_subgraph_sort_dict = defaultdict(list)
    ctg_subgraphList_dict = defaultdict()
    subgraphGroup_RE_dict = defaultdict(tuple)
    
    subgraphGroup_idx = 1
    for graph in graphs_dict.values():

        graph_nodes_list = list(graph.nodes())

        if len(graph_nodes_list) == 0:
            continue

        ctgs = list(set(graph_nodes_list) & set(ctgs_list))
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


    return filter_subgraph_sort_dir_dict, filter_subgraph_sort_dict


def get_AGP(filter_subgraph_sort_dir_dict, ctg_RE_dict):

    with open("subgraph_sort.agp", 'w') as file:
        for Group_idx,(subgraphGroup, ctgs_dir) in enumerate(filter_subgraph_sort_dir_dict.items(), 1):
            len_sum = 0
            for idx, (ctg, ctg_dir) in enumerate(ctgs_dir):
                ctg_len = int(ctg_RE_dict[ctg][1])
                file.write(f"scaffold_{Group_idx}\t{len_sum+1}\t{ctg_len+len_sum}\t{(idx+1)*2-1}\t{'W'}\t{ctg}\t{1}\t{ctg_RE_dict[ctg][1]}\t{str('+' if str(ctg_dir)=='1' else '-')}\n")
                if idx != len(ctgs_dir)-1:
                    file.write(f"scaffold_{Group_idx}\t{ctg_len+len_sum+1}\t{ctg_len+len_sum+100}\t{(idx+1)*2}\t{'U'}\t{100}\t{'scaffold'}\t{'yes'}\t{'proximity_ligation'}\n")
                    len_sum += ctg_len + 100


def Get_subgraph_scaffold(digraph_file, REFile, hic_link_file, group_file, subgraph_file, gfa_filePath):
    digraph_dict = read_digraph(digraph_file)
    ctg_RE_dict = read_RE(REFile)
    hic_links_dict, hic_nei_dict = read_l(hic_link_file)

    ctgs_list = read_group(group_file)

    subgraph_ctgs_dict, ctg_subgraph_dict = read_subgraph(subgraph_file)
    filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict = filter_subgraph(subgraph_ctgs_dict, ctg_subgraph_dict, ctgs_list)

    subgraph_connect_dict = get_subgraph_link(digraph_dict, subgraph_ctgs_dict, ctg_subgraph_dict)

    topo_order_list = connect_subgraph(filter_subgraph_ctgs_dict, subgraph_connect_dict)
    graphs_dict, global_digraph = get_filter_subgraph_digraph(filter_subgraph_ctgs_dict, filter_ctg_subgraph_dict, gfa_filePath, topo_order_list, digraph_dict)
    filter_subgraph_sort_dir_dict, filter_subgraph_sort_dict = get_subgraph_group_inner_sort(graphs_dict, ctgs_list, ctg_RE_dict, global_digraph, hic_links_dict, hic_nei_dict)

    get_AGP(filter_subgraph_sort_dir_dict, ctg_RE_dict)




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


    args = parser.parse_args()
    group_file = args.group_file
    subgraph_file = args.subgraph_file
    gfa_filePath = args.graph
    REFile = args.REFile
    digraph_file = args.digraph_file
    hic_link_file = args.hic_link

    Get_subgraph_scaffold(digraph_file, REFile, hic_link_file, group_file, subgraph_file, gfa_filePath)