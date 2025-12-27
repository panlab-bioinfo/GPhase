import sys
import pysam
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import tempfile
from multiprocessing import Pool
from contextlib import ExitStack


def rename_agp_duplicate_utg(agp_file, out_agp_file):
    """
    Process AGP file to detect truly duplicated contigs/utgs (used in different scaffolds or non-contiguously).
    Works regardless of contig name prefix.
    Renames duplicated instances and returns instance_map: {original_name: [new_name1, new_name2, ...]}
    """
    current_object = None
    current_name = None
    current_part = 0

    instance_counts = defaultdict(int)
    instance_map = defaultdict(list)

    with open(agp_file) as fin, open(out_agp_file, "w") as fout:
        for line in fin:
            if line.startswith("#") or line.strip() == "":
                fout.write(line)
                continue

            parts = line.strip().split("\t")
            if len(parts) < 6:
                fout.write(line)
                continue

            object_name = parts[0]
            comp_type = parts[4]
            comp_id = parts[5]
            part_num = int(parts[3])

            if comp_type.upper() in {"W", "D", "F", "A"}:
                is_new_instance = False
                if (object_name != current_object or
                    comp_id != current_name or
                    part_num != current_part + 1):
                    is_new_instance = True

                if is_new_instance:
                    instance_counts[comp_id] += 1

                count = instance_counts[comp_id]
                new_name = comp_id if count == 1 else f"{comp_id}_dup{count-1}"

                instance_map[comp_id].append(new_name)
                parts[5] = new_name

                current_object = object_name
                current_name = comp_id
                current_part = part_num

            fout.write("\t".join(parts) + "\n")

    return instance_map


# ====================== File Splitting Functions ======================

def split_pairs_file(pairs_file, num_chunks, temp_dir):
    """Split .pairs file into roughly equal chunks using hash-based distribution."""
    chunk_files = [os.path.join(temp_dir, f"chunk_{i:03d}.pairs") for i in range(num_chunks)]

    with open(pairs_file) as fin, ExitStack() as stack:
        writers = [stack.enter_context(open(cf, "w")) for cf in chunk_files]

        # Extract and preserve original header
        header = ""
        for line in fin:
            if line.startswith("#"):
                header += line
                continue
            break

        # Write standardized header to all chunks
        for w in writers:
            w.write(header + "# read_name contig1 pos1 contig2 pos2\n")

        # Distribute lines
        fin.seek(0)
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            idx = hash(line) % num_chunks
            writers[idx].write(line)

    return chunk_files


def split_bam_file(bam_file, num_chunks, temp_dir):
    """Split BAM into chunks with roughly equal valid proper pairs."""
    chunk_files = [os.path.join(temp_dir, f"chunk_{i:03d}.bam") for i in range(num_chunks)]

    # Get header
    bam_in = pysam.AlignmentFile(bam_file, "rb")
    header_dict = bam_in.header.to_dict()
    bam_in.close()

    writers = [pysam.AlignmentFile(cf, "wb", header=header_dict) for cf in chunk_files]

    chunk_idx = 0
    total_reads = 0
    valid_reads = 0

    for read in pysam.AlignmentFile(bam_file, "rb"):
        total_reads += 1
        if (read.is_secondary or read.is_supplementary or
            read.is_unmapped or read.mate_is_unmapped or
            not read.is_proper_pair or read.is_read2 or
            read.is_reverse == read.mate_is_reverse or
            read.next_reference_start < 0):
            continue

        writers[chunk_idx].write(read)
        valid_reads += 1
        chunk_idx = (chunk_idx + 1) % num_chunks

    for w in writers:
        w.close()

    print(f"BAM split: {total_reads:,} total alignments → {valid_reads:,} valid pairs → {num_chunks} chunks")
    return chunk_files


# ====================== Parallel Workers ======================

def process_pairs_chunk(task):
    chunk_file, instance_map, out_file = task
    written = 0
    with open(chunk_file) as fin, open(out_file, "w") as fout:
        header_written = False
        for line in fin:
            if line.startswith("#"):
                if not header_written:
                    fout.write("# read_name contig1 pos1 contig2 pos2\n")
                    header_written = True
                continue
            if not line.strip():
                continue

            parts = line.strip().split()
            if len(parts) < 5:
                continue

            read_name, contig1, pos1, contig2, pos2 = parts[:5]
            extra = parts[5:]

            copies1 = instance_map.get(contig1, [contig1])
            copies2 = instance_map.get(contig2, [contig2])

            for c1 in copies1:
                for c2 in copies2:
                    fout.write("\t".join([read_name, c1, pos1, c2, pos2] + extra) + "\n")
                    written += 1
    return written


def process_bam_chunk(task):
    chunk_file, instance_map, out_file = task
    written = 0
    with pysam.AlignmentFile(chunk_file, "rb") as bam, open(out_file, "w") as fout:
        fout.write("# read_name contig1 pos1 contig2 pos2\n")
        for read in bam:
            if (read.is_secondary or read.is_supplementary or
                read.is_unmapped or read.mate_is_unmapped or
                not read.is_proper_pair or read.is_read2 or
                read.is_reverse == read.mate_is_reverse or
                read.next_reference_start < 0):
                continue

            chr1 = read.reference_name
            chr2 = read.next_reference_name
            pos1 = read.reference_start + 1
            pos2 = read.next_reference_start + 1
            read_name = read.query_name

            copies1 = instance_map.get(chr1, [chr1])
            copies2 = instance_map.get(chr2, [chr2])

            for c1 in copies1:
                for c2 in copies2:
                    fout.write(f"{read_name}\t{c1}\t{pos1}\t{c2}\t{pos2}\n")
                    written += 1
    return written


# ====================== Main Parallel Processing ======================

def process_hic_parallel(hic_file, instance_map, out_pairs_file, num_processes=None):
    """Parallel Hi-C processing with automatic temporary file cleanup."""
    if num_processes is None:
        num_processes = max(1, os.cpu_count() - 1)

    # TemporaryDirectory automatically deletes everything on exit (even on exception)
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Using temporary directory: {temp_dir}")

        # 1. Split input
        if hic_file.endswith(".bam"):
            print("Splitting BAM file...")
            chunk_files = split_bam_file(hic_file, num_processes, temp_dir)
            worker = process_bam_chunk
        elif hic_file.endswith(".pairs"):
            print("Splitting .pairs file...")
            chunk_files = split_pairs_file(hic_file, num_processes, temp_dir)
            worker = process_pairs_chunk
        else:
            print("Error: Unsupported format. Use .bam or .pairs")
            sys.exit(1)

        # 2. Parallel duplication
        out_chunk_files = [os.path.join(temp_dir, f"dup_chunk_{i:03d}.pairs") for i in range(num_processes)]
        tasks = list(zip(chunk_files, [instance_map] * num_processes, out_chunk_files))

        total_written = 0
        print(f"Launching {num_processes} processes for duplication...")
        with Pool(num_processes) as pool:
            for written in pool.imap_unordered(worker, tasks):
                total_written += written
                print(f"Chunk completed → current total interactions: {total_written:,}")

        # 3. Merge into final file
        print("Merging results into final .pairs file...")
        with open(out_pairs_file, "w") as final_out:
            final_out.write("# read_name contig1 pos1 contig2 pos2\n")
            first = True
            for cf in out_chunk_files:
                with open(cf) as f:
                    for line in f:
                        if first and line.startswith("#"):
                            continue
                        final_out.write(line)
                first = False

        print(f"Hi-C processing completed!")
        print(f"Final file: {out_pairs_file}")
        print(f"Total duplicated interactions: {total_written:,}")

    # ← temp_dir and ALL its contents are automatically deleted here


def duplicate_fasta_sequences(fasta_file, instance_map, out_fasta_file):
    """Duplicate FASTA entries based on AGP renaming."""
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    out_records = []

    seen = set()
    for orig, names in instance_map.items():
        if orig not in seq_dict:
            print(f"Warning: {orig} not found in FASTA", file=sys.stderr)
            continue
        seq = seq_dict[orig].seq
        for name in names:
            if name in seen:
                print(f"Warning: duplicate name {name}", file=sys.stderr)
            seen.add(name)
            out_records.append(SeqRecord(seq, id=name, description=""))

    for uid, rec in seq_dict.items():
        if uid not in instance_map:
            out_records.append(rec)

    SeqIO.write(out_records, out_fasta_file, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Disambiguate duplicated contigs in assembly (AGP + FASTA + optional Hi-C). "
                    "Handles any contig prefix and cleans up all temporary files automatically."
    )
    parser.add_argument("agp_file", help="Input AGP file")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--no-hic", action="store_true", help="Skip Hi-C processing")
    parser.add_argument("--hic-file", help="Hi-C file (.bam or .pairs). Required unless --no-hic")
    parser.add_argument("--processes", type=int, default=16,
                        help="Number of processes (default: 16)")

    args = parser.parse_args()

    out_agp = args.output_prefix + ".agp"
    out_fa = args.output_prefix + ".fa"
    out_pairs = args.output_prefix + ".pairs"

    print("Step 1: Processing AGP and renaming duplicated contigs...")
    instance_map = rename_agp_duplicate_utg(args.agp_file, out_agp)

    if args.no_hic:
        print("Skipping Hi-C processing (--no-hic)")
    else:
        if not args.hic_file:
            parser.error("--hic-file required unless --no-hic is used")
        print("Step 2: Parallel Hi-C duplication (temporary files will be auto-deleted)...")
        process_hic_parallel(args.hic_file, instance_map, out_pairs, num_processes=args.processes)

    print("Step 3: Duplicating FASTA sequences...")
    duplicate_fasta_sequences(args.fasta_file, instance_map, out_fa)

    print("\n=== All Done! Temporary files have been automatically cleaned up ===")
    print("Output files:")
    print(f"  AGP:    {out_agp}")
    print(f"  FASTA:  {out_fa}")
    if not args.no_hic:
        print(f"  PAIRS:  {out_pairs}")


if __name__ == "__main__":
    main()