from collections import defaultdict, deque
import csv
import networkx as nx
import argparse



def read_gfa(gfa_filePath):

    utgs_list = list()
    # graph = nx.MultiDiGraph()
    graph = nx.Graph()
    with open(gfa_filePath, 'r') as file:
        for line in file:
            line = line.strip().split()
            if(line[0] == "S"):
                utgs_list.append(line[1])
                # graph.add_node(line[1], length = len(line[2]), seq = line[2], coverage = int(line[4][5:]), visit_count = 1, out=set(), enter=set(), type=None)
                graph.add_node(line[1], length = len(line[2]), seq = line[2], visit_count = 1, out=set(), enter=set(), type=None)
            elif(line[0] == "L"): # no self loop
                # gfa.add_edge(line[1],line[3])
                graph.add_edge(line[1],line[3],strand1=line[2], strand2=line[4], match=int(line[5][:-1]))
                if line[2] == '+' :
                    graph.nodes[line[1]]['out'].add(line[3]) 
                elif line[2] == '-':
                    graph.nodes[line[1]]['enter'].add(line[3])

    return graph, utgs_list



#将无向图转换为有向图
def trans_digraph(graph, digraph, start_node):
    # 双向传播
    pair = (start_node, 1)

    queue = deque([pair])  
 
    visited = {start_node}



    while queue:
        # + : 1 : out ; - : 0 : enter
        # 通过 比较邻居节点的交集确认连接节点的方向 

        (utg, utg_dir) = queue.popleft() 

        nei_list_1 = graph.nodes[utg]['enter'] if utg_dir == 0 else graph.nodes[utg]['out']
        nei_list_2 = graph.nodes[utg]['enter'] if utg_dir == 1 else graph.nodes[utg]['out']

        nei_dir_list_1 = set([ (nei_utg, 1) if utg in graph.nodes[nei_utg]['enter'] else (nei_utg, 0) for nei_utg in nei_list_1 ] )
        nei_dir_list_2 = set([ (nei_utg, 1) if utg in graph.nodes[nei_utg]['out'] else (nei_utg, 0) for nei_utg in nei_list_2 ] )

        queue.extend([ (nei_utg, nei_dir) for nei_utg, nei_dir in nei_dir_list_1 if nei_utg not in visited ])
        queue.extend([ (nei_utg, nei_dir) for nei_utg, nei_dir in nei_dir_list_2 if nei_utg not in visited ])
        visited.update([ nei_utg for nei_utg, nei_dir in nei_dir_list_1 ])
        visited.update([ nei_utg for nei_utg, nei_dir in nei_dir_list_2 ])


        for pair in nei_dir_list_1:
            digraph.add_edge(utg, pair[0])

        for pair in nei_dir_list_2:
            digraph.add_edge(pair[0], utg)

    return digraph


def run_trans_digraph(gfa_filePath, gfa_file, flag_output_csv=False):

    digraph = nx.DiGraph()
    graph, utgs_list = read_gfa(gfa_filePath)
    connected_components = nx.connected_components(graph)

    for i, component in enumerate(connected_components):
        if len(component) == 1:
            digraph.add_node(list(component)[0])
        else:
            digraph = trans_digraph(graph, digraph, list(component)[0])

    if flag_output_csv:
        with open(f"{gfa_file}.digraph.csv", 'w', newline='') as file:
            writer = csv.writer(file)
            # 写入CSV文件头部
            writer.writerow(['source', 'target'])

            # 遍历图中的所有边并写入CSV文件
            for u, v in digraph.edges():
                writer.writerow([u, v])

    return digraph
    

if __name__ == '__main__':
    

    parser = argparse.ArgumentParser(description="trans GFA", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gfa_file', required=True, help='gfa file')
    parser.add_argument('-o', '--output_prefix', required=True, help='output prefix')


    args = parser.parse_args() 
    gfa_file = args.gfa_file
    output_prefix = args.output_prefix

    run_trans_digraph(gfa_file, output_prefix)




