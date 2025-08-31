# GPhase: A Phasing Assembly Tool Leveraging an Assembly Graph and Hi-C Data

GPhase leverages an assembly graph and Hi-C data to facilitate genome assembly phasing, automatically resolves and assigns collapsed sequences, and fills assembly gaps based on the graph structure.
---
# Installation
To install GPhase, follow these steps:
```
git clone https://github.com/panlab-bioinfo/GPhase.git
cd GPhase
git submodule update --init  --recursive
conda env create -f gphase_environment.yml
conda activate gphase
./gphase -h
```

# Aligning Hi-C data to assembly
To align Hi-C reads to the genome assembly, use chromap or other alignment tools, such as bwa. When using chromap, if the default MAPQ parameters do not produce satisfactory results, you can lower the --MAPQ-threshold value to include more Hi-C alignment information. When using other alignment software, you need to sort the bam file.
```
chromap -i -r asm.fa -o index
chromap --preset hic -x index -r asm.fa -q 0 \
    -1 HiC_1.fq.gz -2 HiC_2.fq.gz \
    --remove-pcr-duplicates -t 64 --SAM -o map.chromap.sam
samtools view -@ 64 -bh -o map.chromap.bam map.chromap.sam
```

# Estimating of the number of contig collapses based on HiC data and popCNV
The popCNV_pipeline.sh script estimates the copy number of collapsed contigs collapse based on HiFi data or NGS data using the popCNV software. The file used by popCNV for GPhase input is `collapse_num.txt` : popcnv/06.genes.round.cn. For details, see [popCNV](https://github.com/sc-zhang/popCNV)
```
/path/to/GPhase/pipeline/popCNV_pipeline.sh \
-f asm.fa \
-r reads.fq.gz \
-p output_prefix \
-t 32
```


# Running the GPhase scaffolding pipeline
1. `asm.fa` :  Genome assembly file in FASTA format (unitigs).
2. `p_utg.gfa` : Assembly graph file in GFA format.
3. `collapse_num.txt` : File with contig collapse information (from popCNV: popcnv/06.genes.round.cn).
4. `map.chromap.bam` : BAM file with mapped Hi-C reads.
5. `n_chr` : Number of chromosomes.
6. `n_hap` : Number of haplotypes.
7. `p` : Prefix for output files.
```
/path/to/GPhase/gphase pipeline\
 -f asm.fa \
 -g genome.bp.p_utg.gfa \
 -c collapse_num.txt \
 -m map.chromap.bam \
 --n_chr 12 \
 --n_hap 4 \
 -p output_prefix
```

# Output file
- `gphase_final.agp` : unitig level assembly result agp file
- `gphase_final.fasta` : unitig level assembly result fasta file
- `gphase_final_rescue.agp` :  unitig level assembly result agp file after rescue
- `gphase_final_ctg2utg.txt` : correspondence between unitig and contig
- `gphase_final_contig.fasta` : Contig-level fasta sequence
- `gphase_final_contig.agp` : contig level assembly result agp file
- `gphase_final_contig_scaffold.fasta` : contig level assembly result fasta file