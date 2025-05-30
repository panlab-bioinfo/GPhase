import csv
import pysam
from collections import defaultdict
from multiprocessing import Pool, cpu_count

def process_read(read, bamfile):
    if read.mapping_quality <= 0:
        return None
    contig1 = bamfile.get_reference_name(read.reference_id)
    contig2 = bamfile.get_reference_name(read.next_reference_id)
    if contig1 is None or contig2 is None:
        return None
    return tuple(sorted([contig1, contig2]))


def filter_and_count_bam(input_bam, output_file):
    connections = defaultdict(int)
    seen_reads = set()

    with pysam.AlignmentFile(input_bam, 'rb') as bamfile:
        # Use multiprocessing for faster processing
        with Pool(cpu_count()) as pool:
            results = pool.starmap(
                process_read, [(read, bamfile) for read in bamfile if read.query_name not in seen_reads]
            )

        # Filter None results and update counts
        for result in filter(None, results):
            connections[result] += 1

    # Write results to CSV
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Contig1", "Contig2", "Connections"])
        for (contig1, contig2), count in connections.items():
            writer.writerow([contig1, contig2, count])

input_bam = "map.chromap.q0.bam"  # Change to your BAM file path
output_file = "output.csv"  # Change to your desired output file path
filter_and_count_bam(input_bam, output_file)
