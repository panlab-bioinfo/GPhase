#!/usr/bin/env python
import argparse
from pop_cnv.pipeline.pipeline import Pipeline


class Main:
    def __init__(self):
        self.opts = self.get_opts()
        self.pipeline = Pipeline(self.opts)
        self.pipeline.run()

    @staticmethod
    def get_opts():
        group = argparse.ArgumentParser()
        group.add_argument('-g', '--genome', help='Input genome fasta file', required=True)
        group.add_argument('-b', '--bam', help='Directory of bam files', required=True)
        group.add_argument('-s', '--win_size', type=int, help='window size, default=1000', default=1000)
        group.add_argument('-l', '--list', help='List file contain 4 columns: chromosome_name start_position '
                                                'end_position gene_name',
                           required=True)
        group.add_argument('-w', '--workdir', help='Working directory', required=True)
        group.add_argument('--group', help='Group file contain 2 columns: sample_name group_name(without space)',
                           required=True)
        group.add_argument('--wild', help='Wild group id', required=True)
        group.add_argument('-t', '--threads', type=int, help='Threads', default=10)
        return group.parse_args()


if __name__ == "__main__":
    Main()
