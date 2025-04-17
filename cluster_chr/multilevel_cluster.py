import igraph as ig
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import networkx as nx
import community as louvain
# import argcomplete


def run(csv_file, output_file, resolution):

    df = pd.read_csv(csv_file)
    df['source'] = df['source'].astype(str)
    df['target'] = df['target'].astype(str)
    df['links'] = df['links'].astype(float)
    g = ig.Graph()  
    nodes = list(set(df['source']).union(set(df['target'])))
    g.add_vertices(nodes)

    for i, row in df.iterrows():
        source = row['source']
        target = row['target']
        weight = row['links']
        g.add_edge(source, target, weight=weight)

    communities = g.community_multilevel(weights='weight',resolution=float(resolution)) 

    # palette = plt.get_cmap("tab20")
    # community_colors = [palette(i / len(communities)) for i in range(len(communities))]
    # node_colors = [community_colors[communities.membership[v]] for v in range(len(g.vs))]
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

    # argcomplete.autocomplete(parser)
    args = parser.parse_args()
    run(args.csv_file, args.output_file, args.resolution)
