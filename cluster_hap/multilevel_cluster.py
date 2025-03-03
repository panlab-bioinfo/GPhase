import igraph as ig
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import networkx as nx
import community as louvain
from collections import defaultdict
# import argcomplete

def read_REs(REFile):
    ctg_RE_len = defaultdict(tuple)
    with open(REFile, 'r') as fp:
        for line in  fp:
            if line[0] == "#":
                continue
            line = line.strip().split()
            ctg_RE_len[line[0]] = (int(line[1]), int(line[2]))
    return ctg_RE_len

def run(csv_file, output_file, resolution,check, RE_file):

    df = pd.read_csv(csv_file)
    df['source'] = df['source'].astype(str)
    df['target'] = df['target'].astype(str)
    df['links'] = df['links'].astype(float)
    g = ig.Graph()  
    nodes = list(set(df['source']).union(set(df['target'])))
    g.add_vertices(nodes)

    # 添加边和节点，并设置权重
    for i, row in df.iterrows():
        source = row['source']
        target = row['target']
        weight = row['links']
        g.add_edge(source, target, weight=weight)



    # 图的联通分量个数
    # components = g.connected_components()
    # print(components)

    # g = components.subgraphs()[0]



    # 使用community_multilevel进行社区检测
    communities = g.community_multilevel(weights='weight',resolution=float(resolution)) 

    # palette = plt.get_cmap("tab20")
    # community_colors = [palette(i / len(communities)) for i in range(len(communities))]

    # # 创建颜色映射
    # node_colors = [community_colors[communities.membership[v]] for v in range(len(g.vs))]

    # 绘制图形并保存
    # g.vs['label'] = nodes
    # layout = g.layout_fruchterman_reingold(
    #     niter=1000,       # 迭代次数
    #     dim=2,
    #     # weights='weight'
    # )
    # visual_style = {
    #     "vertex_color": node_colors,
    #     "vertex_size": 20,
    #     "edge_width": 1,
    #     "layout": layout,
    #     "bbox": (1200, 800),
    #     "margin": 20
    # }
    # plot = ig.plot(g, **visual_style)
    # plot.save("clusters.png")



    if check:
        # 检查有效聚类簇数目
        # 阈值设置为平均聚类簇长度的 1/10
        ctg_RE_len = read_REs(RE_file)
        group_len_dict = defaultdict()

        for idx, community in enumerate(communities):
            utgs = [g.vs[i]['name'] for i in communities[idx]]
            utgs_len = sum([ int(ctg_RE_len[utg][1]) for utg in utgs])
            group_len_dict[idx] = utgs_len
        
        all_len = sum(group_len_dict.values())

        save_group = [ k for k,v in group_len_dict.items() if v > float(all_len)/10]

        with open(output_file, 'w') as file:
            flag = 0
            for idx, community in enumerate(communities):
                if idx in save_group:
                    flag += 1
                    utgs = [g.vs[i]['name'] for i in communities[idx]]
                    file.write(f"group{flag}\t{len(utgs)}\t{' '.join(utgs)}\n")



    else:
        with open(output_file, 'w') as file:
            for idx, community in enumerate(communities):
                utgs = [g.vs[i]['name'] for i in communities[idx]]
                file.write(f"group{idx+1}\t{len(communities[idx])}\t{' '.join(utgs)}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="multilevel clustering algorithm is used and cluster file is generated")
    parser.add_argument('-c', '--csv_file', required=True,
                        help='<filepath> csv file of the hic signal')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath> output file')
    parser.add_argument('-r', '--resolution', default=1,
                    help='<float> The resolution parameter of multilevel cluster algorithm, the larger the value, the more clusters with fewer nodes')
    parser.add_argument('--check', action='store_true', help='get the number of vaild clusters')
    parser.add_argument('--RE_file',help='File required when --check is enabled')

    # argcomplete.autocomplete(parser)
    args = parser.parse_args()

    if args.check and args.RE_file is None:
        parser.error("--RE_file is required when --check is enabled")

    run(args.csv_file, args.output_file, args.resolution, args.check, args.RE_file)
