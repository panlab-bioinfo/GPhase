Gphase: A phasing assembly tool using assembly graphs and Hi-C data
---
# Installation

```
git clone https://gitee.com/qztanging/Gphase.git
conda env create -f Gphase/gphase_environment.yaml
conda activate gphase
/path/to/Gphase/pipeline/gphase_pipeline.sh -h
```

# Aligning Hi-C data to assembly
#### tips: When the default MAPQ parameters of chromap do not work well, you can appropriately lower the MAPQ-threshold parameter value to obtain more Hi-C alignment data.

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
/path/to/Gphase/pipeline/popCNV_pipeline.sh \
-f asm.fa \
-p output_prefix \
-t 32
```
`collapse_num.txt` : popcnv/06.genes.round.cn

# Installing agptools
[agptools](https://github.com/WarrenLab/agptools)
```
cd /path/to/Gphase/src/agptools
pip install .
```


# Running the Gphase scaffolding pipeline
| Scaffolding Pipeline : subGraph_scaffold + YaHs [YaHS](https://github.com/c-zhou/yahs) + HapHiC [HapHiC](https://github.com/zengxiaofei/HapHiC)

1. `asm.fa` :  your genome assembly file in FASTA format (Unitigs).
2. `p_utg.gfa` : the GFA file representing the assembly graph.
3. `collapse_num.txt` : the file that number information for collapse unitigs (popCNV output file: popcnv/06.genes.round.cn).
4. `map.chromap.pairs` : the mapping file used to map the Hi-C reads.
5. `n_chr` : the number of chromosomes.
6. `n_hap` : the number of haplotypes.
7. `p` : The prefix for the output files.

```
/path/to/Gphase/gphase pipeline\
 -f asm.fa \
 -g genome.bp.p_utg.gfa \
 -c collapse_num.txt \
 -m map_file.pairs \
 --n_chr 12 \
 --n_hap 4 \
 -p output_prefix
```

# Output file
path: gphase_output/gphase_final
file:   
1. `*.scaffolds.fasta` : final genome assembly fasta file
2.  `*.final.agp` : final genome assembly agp file



