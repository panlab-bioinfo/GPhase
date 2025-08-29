import pandas as pd
import argparse

def process_chromap_pairs(input_file, output_prefix, cluster_q):

    df = pd.read_csv(
        input_file,
        sep='\t',
        comment='#',
        header=None,
        usecols=[1, 3, 5],  
        dtype={1: str, 3: str, 5: float}  
    )

    filtered = df[(df[1] != df[3]) & (df[5] >= cluster_q)]

    counts = filtered.groupby([1, 3], sort=False).size().reset_index(name='count')

    output_file = f"{output_prefix}.map.links.csv"
    counts.to_csv(output_file, index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter chromap pairs and generate links CSV.")
    parser.add_argument('-i', '--input_file', required=True,
                        help='Input .chromap.pairs file path')
    parser.add_argument('-o', '--output_prefix', required=True,
                        help='Prefix for output CSV (e.g., output.chromap.links.csv)')
    parser.add_argument('-q', '--cluster_q', type=float, default=1,
                        help='Threshold for filtering on column 6 (default: 1)')

    args = parser.parse_args()
    process_chromap_pairs(
        args.input_file,
        args.output_prefix,
        args.cluster_q,
    )
