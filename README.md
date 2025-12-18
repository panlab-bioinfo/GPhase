# GPhase: A Phasing Assembly Tool Leveraging an Assembly Graph and Hi-C/Pore-C Data

GPhase leverages an assembly graph and Hi-C/Pore-C data to facilitate genome assembly phasing, automatically resolves and assigns collapsed sequences, and fills assembly gaps based on the graph structure.
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

# Mapping Hi-C data to assembly
GPhase supports multiple data types, including Hi-C and Pore-C. It also supports their pairs(pa5) and bam(BAM) format mapping files.

Hi-C reads you can using Chromap or other mapping tools, such as BWA, can be used. When using Chromap, if the default MAPQ parameters do not produce satisfactory results, the `--MAPQ-threshold` value can be lowered to include more Hi-C mapping information. When using other mapping software, the BAM files need to be sorted.
```
chromap -i -r asm.fa -o index
chromap --preset hic -x index -r asm.fa -q 0 \
    -1 HiC_1.fq.gz -2 HiC_2.fq.gz \
    --remove-pcr-duplicates -t 64 -o map.chromap.pairs
```
To process Pore-C data, you can use [PPL Toolbox](https://github.com/versarchey/PPL-Toolbox). GPhase provides a script to run PPL Toolbox. You can quickly run it to get the final pairs file `map.PPL.pairs` and input it into GPhase.
```
/path/to/GPhase/pipeline/PPL_pipeline.sh -j /path/to/PPL-Toolbox.jar \
-f asm.fa \
-r reads.fq.gz \
-p PPL
```


# Estimating of the number of contig collapses based on HiFi data and popCNV
The popCNV_pipeline.sh script estimates the copy number of collapsed contigs collapse based on HiFi data using the popCNV software. The file used by popCNV for GPhase input is `collapse_num.txt` : popcnv/06.genes.round.cn. For details, see [popCNV](https://github.com/sc-zhang/popCNV)
```
/path/to/GPhase/pipeline/popCNV_pipeline.sh \
-f asm.fa \
-p output_prefix \
-t 32 -r reads.fq.gz
```


# Running the GPhase scaffolding pipeline
1. `asm.fa` :  Genome assembly file in FASTA format (unitigs).
2. `p_utg.gfa` : Assembly graph file in GFA format.
3. `collapse_num.txt` : File with contig collapse information (from popCNV: popcnv/06.genes.round.cn).
4. `map.chromap.bam` : pairs/bam file with mapped Hi-C reads.
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
GPhase will output a folder named gphase_output, which will generate the following four folders in sequence.
- `preprocessing` : Data preprocessing
- `cluster_chr` : Results of chromosome clustering
- `cluster_hap` : Haplotype clustering results within each chromosome
- `scaffold_hap` : Scaffolding results for each haplotype within each chromosome

# Final assembly result
- `gphase_final.agp` : unitig level assembly result agp file
- `gphase_final.fasta` : unitig level assembly result fasta file
- `gphase_final_rescue.agp` :  unitig level assembly result agp file after rescue
- `gphase_final_ctg2utg.txt` : correspondence between unitig and contig
- `gphase_final_contig.fasta` : Contig-level fasta sequence
- `gphase_final_contig.agp` : contig level assembly result agp file
- `gphase_final_contig_scaffold.fasta` : contig level assembly result fasta file

# Test dataset
To help you quickly verify the functionality of the software, we provide a small test dataset. This dataset contains input data that demonstrates the core functionality of the software. You can download it from this link https://drive.google.com/drive/folders/1M_ZlSHBTDwtCHGrUI6uMCVutfIweECaY?usp=sharing

Use the following command to run the test dataset
```
tar -zxvf test_dataset.tar.gz
export PATH=$PATH:/path/to/GPhase
bash run_gphase.sh
```
# Contact
This software is developed by Professor Wei-Hua Pan's team at the Shenzhen Institute of Genome Research, Chinese Academy of Agricultural Sciences. 

If you have any questions or concerns while using the software, please submit an issue in the repository or contact us through the following methods:
### Email:
#### Prof. Pan: panweihua@caas.cn
#### Du Wenjie: duwenjie1024@163.com