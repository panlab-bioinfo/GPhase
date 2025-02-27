# Installation

```
git clone https://gitee.com/qztanging/anhic.git
conda env create -f anhic/anhic_environment.yaml
conda activate anhic
/path/to/anhic/pipeline/anhic_pipeline.sh -h
```

# Aligning Hi-C data to assembly

```
chromap -i -r asm.fa -o index
chromap --preset hic -x index -r asm.fa \
    -1 hic_R1.fq.gz -2 hic_R2.fq.gz \
    --remove-pcr-duplicates -t 80 \
    -o map.chromap.pairs

```

# Estimating of the number of contig collapses based on HiC data and popCNV
By default, popCNV_pipeline.sh will use the file ending with fastaq.gz in the directory as HiFi data. [popCNV](https://github.com/sc-zhang/popCNV)

1. `asm.fa` : your genome assembly file in FASTA format (Unitigs).
2. `p` : The prefix for the output files.
3. `t` : The number of threads.
```
/path/to/anhic/pipeline/popCNV_pipeline.sh \
-f asm.fa \
-p output_prefix \
-t 32
```
`collapse_num.txt` : popcnv/06.genes.round.cn



# Running the AnHiC scaffolding pipeline
| Scaffold Pipeline : subGraph_scaffold + YaHs [YaHS](https://github.com/c-zhou/yahs) + HapHiC [HapHiC](https://github.com/zengxiaofei/HapHiC)

1. `asm.fa` :  your genome assembly file in FASTA format (Unitigs).
2. `p_utg.gfa` : the GFA file representing the assembly graph.
3. `collapse_num.txt` : the file that number information for collapse unitigs (popCNV output file: popcnv/06.genes.round.cn).
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
 --n_chr 12 \
 --n_hap 4 \
 -p output_prefix
```
