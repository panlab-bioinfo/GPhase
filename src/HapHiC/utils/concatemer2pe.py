#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created at: 2026-01-29 09:51

import argparse
from collections import defaultdict
from copy import copy
from typing import List, DefaultDict, Optional

from pysam import AlignmentFile, AlignedSegment # type: ignore[import]


def convert_concatemer_to_paired_end_bam(
    bam: str,
    output: str,
    mapq: int,
    percent_identity: float,
    alignment_length: int,
    threads: int
) -> None:
    """
    Convert mapped concatemer long reads (Pore-C/CiFi) in BAM format to pairwise
    alignments (Hi-C-like) in BAM format.
    """

    def update_current_aln_list(
        current_aln_list: List[List[AlignedSegment]],
        aln: AlignedSegment
    ) -> None:
        """
        Update the current alignment list.
        """
        if aln.is_supplementary:
            aln.flag -= 2048
            current_aln_list.append([aln])
        elif aln.is_secondary:
            current_aln_list[-1].append(aln)
        else:
            current_aln_list.append([aln])

    def generate_pairwise_alignments(
        current_aln_list: List[List[AlignedSegment]],
        aln: AlignedSegment,
        fout: AlignmentFile
    ) -> None:
        """
        Generate Hi-C-like alignments from the current alignment list and
        write them to the output BAM file.
        """
        pri_and_sup_aln_num: int = len(current_aln_list)
        if pri_and_sup_aln_num >= 2:
            xa_tag_dict: DefaultDict[int, str] = defaultdict(str)
            for n, aln_list in enumerate(current_aln_list):
                if len(aln_list) > 1:
                    for sec_aln in aln_list[1:]:
                        strand: str = '+' if sec_aln.is_forward else '-'
                        xa_tag: str = (
                            f'{sec_aln.reference_name},{strand}{sec_aln.reference_start+1},'
                            f"{sec_aln.cigarstring},{sec_aln.get_tag('NM')};"
                        )
                        xa_tag_dict[n] += xa_tag
            k: int = 0
            for i in range(pri_and_sup_aln_num):
                for j in range(i+1, pri_and_sup_aln_num):
                    aln_i: AlignedSegment = current_aln_list[i][0]
                    aln_j: AlignedSegment = current_aln_list[j][0]
                    aln_i_copy: AlignedSegment = copy(aln_i)
                    aln_j_copy: AlignedSegment = copy(aln_j)
                    mock_read_name: str = f'{aln_i.query_name}_read{k}'
                    aln_i_copy.query_name = aln_j_copy.query_name = mock_read_name
                    # Add PAIRED, READ1, READ2, and MREVERSE flags.
                    aln_i_copy.flag += 65 + 2 * aln_j.flag
                    aln_j_copy.flag += 129 + 2 * aln_i.flag
                    # Set next_reference_name
                    aln_i_copy.next_reference_name = aln_j.reference_name
                    aln_j_copy.next_reference_name = aln_i.reference_name
                    # Set next_reference_start
                    aln_i_copy.next_reference_start = aln_j.reference_start
                    aln_j_copy.next_reference_start = aln_i.reference_start
                    # Set SA tag to None
                    aln_i_copy.set_tag('SA', None)
                    aln_j_copy.set_tag('SA', None)
                    if xa_tag_dict:
                        # Set XA tag for secondary alignments
                        if i in xa_tag_dict:
                            aln_i_copy.set_tag('XA', xa_tag_dict[i])
                        if j in xa_tag_dict:
                            aln_j_copy.set_tag('XA', xa_tag_dict[j])
                    fout.write(aln_i_copy)
                    fout.write(aln_j_copy)
                    k += 1
        current_aln_list.clear()
        update_current_aln_list(current_aln_list, aln)

    format_options: List[bytes] = [
        (
            f'filter=!flag.unmap && mapq >= {mapq} && '
            f'[de] <= {(100-percent_identity)/100} && '
            f'(endpos - pos + 1) >= {alignment_length}'
        ).encode()
    ]

    with (
        AlignmentFile(
            bam,
            mode='rb',
            threads=threads,
            format_options=format_options # type: ignore[arg-type]
        ) as fin,
        AlignmentFile(output, mode='wb', threads=threads, template=fin) as fout
    ):
        if fin.header.to_dict().get('HD', {}).get('SO', '') != 'unsorted':
            raise ValueError(
                'The input BAM file should be unsorted.'
            )
        current_read: Optional[str] = None
        current_aln_list: List[List[AlignedSegment]] = []
        for aln in fin:
            if aln.query_name == current_read:
                update_current_aln_list(current_aln_list, aln)
            else:
                generate_pairwise_alignments(current_aln_list, aln, fout)
                current_read = aln.query_name
        generate_pairwise_alignments(current_aln_list, AlignedSegment(), fout)


class CustomHelpFormatter(argparse.HelpFormatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, width=100, **kwargs)


def main() -> None:
    """
    The main function of concatemer2pe.
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description=(
            'Convert mapped concatemer long reads (Pore-C/CiFi) in BAM format to pairwise '
            'alignments (Hi-C-like) in BAM format.'
        ),
        prog='concatemer2pe',
        formatter_class=CustomHelpFormatter
    )
    parser.add_argument(
        'bam',
        help=(
            'Mapped concatemer long reads in BAM format generated by minimap2. Note that the '
            'BAM file should be unsorted.'
        )
    )
    parser.add_argument(
        '--output',
        default='paired.bam',
        help='Filename of the output BAM file, default: %(default)s.'
    )
    parser.add_argument(
        '--mapq',
        type=int,
        default=0,
        help=(
            'MAPQ cutoff, default: %(default)s (range: 0-60). Read pairs with both MAPQ >= this '
            'value will be kept'
        )
    )
    parser.add_argument(
        '--percent-identity',
        type=float,
        default=0,
        help=(
            'Percent identity cutoff, default: %(default)s (range: 0-100). Read pairs with both '
            'percent identity >= this value will be kept'
        )
    )
    parser.add_argument(
        '--alignment-length',
        type=int,
        default=0,
        help=(
            'Alignment length cutoff, default: %(default)s. Read pairs with both alignment length '
            '>= this value will be kept'
        )
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Threads for parsing BAM file, default: %(default)s'
    )
    args: argparse.Namespace = parser.parse_args()

    convert_concatemer_to_paired_end_bam(
        args.bam,
        args.output,
        args.mapq,
        args.percent_identity,
        args.alignment_length,
        args.threads
    )


if __name__ == '__main__':
    main()
