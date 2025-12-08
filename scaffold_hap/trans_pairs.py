import pandas as pd
import re
import argparse


def read_agp(agp_file):
    try:

        agp = pd.read_csv(agp_file, sep='\t', comment='#', header=None,names=['scaffold', 'start', 'end', 'part_num', 'type', 'component_id', 'component_start', 'component_end', 'orientation'])
        if agp.columns[0] != 'scaffold':
            agp.columns = ['scaffold', 'start', 'end', 'part_num', 'type', 'component_id', 'component_start', 'component_end', 'orientation']
        for col in ['start', 'end', 'component_start', 'component_end']:
            agp[col] = pd.to_numeric(agp[col], errors='coerce')
        agp = agp.dropna(subset=['start', 'end', 'component_start', 'component_end'])

        return agp
    except Exception as e:
        raise


def build_contig_to_scaffold_map(agp):
    contig_map = {}
    scaffold_lengths = {}
    
    for _, row in agp.iterrows():
        if row['type'] == 'W':
            contig_id = row['component_id']
            scaffold = row['scaffold']
            scaffold_start = row['start']
            component_start = row['component_start']
            component_end = row['component_end']
            orientation = row['orientation']
            contig_length = component_end - component_start + 1
            contig_map[contig_id] = {
                'scaffold': scaffold,
                'scaffold_start': scaffold_start,
                'orientation': orientation,
                'length': contig_length
            }
        

        scaffold = row['scaffold']
        end = row['end']
        scaffold_lengths[scaffold] = max(scaffold_lengths.get(scaffold, 0), end)
    
    return contig_map, scaffold_lengths


def convert_pairs(pairs_file, contig_map, scaffold_lengths, output_file):
    with open(pairs_file, 'r') as f_in, open(output_file, 'w') as f_out:
        scaffold_set = set()
        for line in f_in:

            if line.startswith('#'):
                match = re.match(r'#(?:#chromsize|chromosome|chr|contig|sequence):\s*(\S+)(?:\s*(\d+))?', line.strip(), re.IGNORECASE)
                if match:
                    contig_id = match.group(1)
                    if contig_id in contig_map:
                        scaffold = contig_map[contig_id]['scaffold']
                        scaffold_length = scaffold_lengths.get(scaffold, 0)
                        if scaffold not in scaffold_set:
                            f_out.write(f"##chromsize: {scaffold} {scaffold_length}\n")
                            scaffold_set.add(scaffold)
                else:
                    continue
                continue

            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue  
            read_id, contig1, pos1, contig2, pos2, mapq = fields[:6]
            try:
                pos1, pos2 = int(pos1), int(pos2)
            except ValueError:
                continue  


            if contig1 in contig_map:
                map1 = contig_map[contig1]
                scaffold1 = map1['scaffold']
                scaffold_pos1 = map1['scaffold_start'] + (pos1 - 1 if map1['orientation'] == '+' else map1['length'] - pos1)
            else:
                continue 

            if contig2 in contig_map:
                map2 = contig_map[contig2]
                scaffold2 = map2['scaffold']
                scaffold_pos2 = map2['scaffold_start'] + (pos2 - 1 if map2['orientation'] == '+' else map2['length'] - pos2)
            else:
                continue  

            new_fields = [read_id, scaffold1, str(int(scaffold_pos1)), scaffold2, str(int(scaffold_pos2)), mapq] 
            f_out.write('\t'.join(new_fields) + '\n')


def Trans_pairs(agp_file, pairs_file, output_file):

    agp = read_agp(agp_file)
    contig_map, scaffold_lengths = build_contig_to_scaffold_map(agp)
    
    convert_pairs(pairs_file, contig_map, scaffold_lengths, output_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert pairs file based on AGP file.")
    
    parser.add_argument("agp_file", help="Path to the input AGP file (e.g., yahs_iter2_scaffolds_final.agp)")
    parser.add_argument("pairs_file", help="Path to the input pairs file (e.g., converted.pairs)")
    parser.add_argument("output_file", help="Path to the output pairs file (e.g., converted2.pairs)")

    args = parser.parse_args()
    Trans_pairs(args.agp_file, args.pairs_file, args.output_file)