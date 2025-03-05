from collections import defaultdict
import numpy as np
import pandas as pd
from kneed import KneeLocator
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import statistics
import argparse


def find_best_knee(csv_file, output_prefix):

    try:
        # 先读取第一行，看看是否存在header
        with open(csv_file, 'r') as file:
            first_line = file.readline()
            has_header = 1 if first_line == "source,target,links\n" else 0
            
        # 如果有header
        if has_header:
            data = pd.read_csv(csv_file, sep=",")
        else:
            # 如果没有header，手动加上默认header
            data = pd.read_csv(csv_file, sep=",", header=None)
            data.columns = ["source", "target", "links"]

    except Exception as e:
        print(f"Error reading the CSV file: {e}")
        return None

    data = pd.read_csv(csv_file, sep=",")
    sort_data = data['links'].sort_values()
    sort_data.index = data.index


    x = sort_data.index
    y = list(sort_data)


    kl = KneeLocator(x, y, curve="convex", direction="increasing", online=True)

    knees_y_median = statistics.median(kl.all_knees_y)
    knees_y_mean = statistics.mean(kl.all_knees_y)


    output_file = open(f"{output_prefix}.result.txt", 'w')
    output_file.write(f"kl.knee_x:{kl.knee}\nkl.knee_y:{kl.knee_y}\nkl.knees_y_median:{knees_y_median}\nkl.knees_y_mean:{knees_y_mean}\n")
    output_file.write(f"all_knees_x:{kl.all_knees}\n")
    output_file.write(f"all_knees_y:{kl.all_knees_y}\n")
    output_file.close()


    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.plot(x, y, zorder=1)
    plt.axvline(x=kl.knee, color='orange', linestyle='--', label='best_knee_x', zorder=2)
    plt.axhline(y=kl.knee_y, color='orange', linestyle='--', label='best_knee_y', zorder=3)
    # plt.axhline(y=knees_y_median, color='red', linestyle='--', label='knees_y_median', zorder=4)
    # plt.axhline(y=knees_y_mean, color='blue', linestyle='--', label='knees_y_mean', zorder=5)

    yticks = plt.gca().get_yticks()
    new_yticks = list(yticks) + [kl.knee_y]
    plt.gca().set_yticks(new_yticks)
    plt.gca().set_yticklabels([f'{tick:.2f}' for tick in new_yticks])



    plt.legend()
    plt.savefig(f'{output_prefix}.png')

    return  float(kl.knee_y)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="find knees")
    parser.add_argument('-c', '--csv_file', required=True,
                        help='<filepath> csv file of the hic signal')
    parser.add_argument('-o', '--output_prefix', required=True,
                        help='<str> output file prefix')
    args = parser.parse_args()
    find_best_knee(args.csv_file, args.output_prefix)


