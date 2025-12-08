import pickle
import pandas as pd
import argparse

# 设置命令行参数解析器
def parse_args():
    parser = argparse.ArgumentParser(description="Convert a pickle file to CSV.")
    parser.add_argument('input_pkl', type=str, help="Path to the input .pkl file")
    parser.add_argument('output_csv', type=str, help="Path to the output .csv file")
    return parser.parse_args()


def pkl_to_csv(input_pkl, output_csv):

    with open(input_pkl, 'rb') as f:
        data = pickle.load(f)


    if isinstance(data, dict) and all(isinstance(k, tuple) and len(k) == 2 and isinstance(v, int) for k, v in data.items()):

        df = pd.DataFrame(list(data.items()), columns=['utg_pair', 'value'])
        
        df[['utg1', 'utg2']] = pd.DataFrame(df['utg_pair'].tolist(), index=df.index)
        df = df.drop(columns=['utg_pair']) 
        
        df = df[['utg1', 'utg2', 'value']] 

    else:
        raise ValueError("Data is not in the expected format (dictionary with tuple keys and integer values).")


    df.to_csv(output_csv, index=False)
    print(f"File has been successfully converted to {output_csv}")


if __name__ == "__main__":
    args = parse_args()
    pkl_to_csv(args.input_pkl, args.output_csv)
