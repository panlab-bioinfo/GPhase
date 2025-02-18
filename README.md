# Installation

```
git clone https://gitee.com/qztanging/anhic.git
conda env create -f anhic/environment.yaml
conda activate anhic
/path/to/anhic/pipeline/anhic_pipeline.sh -h
```

# Align Hi-C data to the assembly

```
chromap -i -r asm.fa -o index
chromap --preset hic -x index -r asm.fa \
    -1 hic_R1.fq.gz -2 hic_R2.fq.gz \
    --remove-pcr-duplicates -t 80 \
    -o map.chromap.pairs

```

# Run HapHiC scaffolding pipeline

1. `asm.fa` :  your genome assembly file in FASTA format.
2. `p_utg.gfa` : the GFA file representing the assembly graph.
3. `collapse_num.txt` : the file that number information for collapse unitigs.
4. `map.chromap.pairs` : the mapping file used to map the Hi-C reads.
5. `n_chr` : the number of chromosomes.
6. `n_hap` : the number of haplotypes.
7. `p` : The prefix for the output files.

```
/path/to/anhic/pipeline/anhic_pipeline.sh \
 -f asm.fa \
 -g genome.bp.p_utg.gfa \
 -c collapse_num.txt \
 -m map_file.pairs \
 --n_hap 12 \
 --n_hap 4 \
 -p output_prefix
```
