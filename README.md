# GPhase: A Phasing Assembly Tool Leveraging an Assembly Graph and Hi-C/Pore-C Data

GPhase leverages an assembly graph and Hi-C/Pore-C data to facilitate genome assembly phasing, automatically resolves and assigns collapsed sequences, and fills assembly gaps based on the graph structure.
---
# Table of contents
* [Installation](#installation)
* [Mapping Hi-C data to assembly](#mapping-hi-c-data-to-assembly)
* [Estimating of the number of contig collapses based on HiFi data and popCNV](#estimating-of-the-number-of-contig-collapses-based-on-hifi-data-and-popcnv)
* [Running the GPhase scaffolding pipeline](#running-the-gphase-scaffolding-pipeline)
* [Output file](#output-file)
* [Final assembly result](#final-assembly-result)
* [Generate a Hi-C heatmap](#generate-a-hi-c-heatmap)
* [Tips](#tips)
* [Test dataset](#test-dataset)
* [Contact](#contact)
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
    --remove-pcr-duplicates -t 64 --SAM -o map.chromap.sam
samtools view -@ 64 -bh map.chromap.sam -o map.chromap.bam
```
To process Pore-C data, you can use [PPL Toolbox](https://github.com/versarchey/PPL-Toolbox). GPhase provides a script to run PPL Toolbox. You can quickly run it to get the final pairs file `map.PPL.pairs` and input it into GPhase.
```
/path/to/GPhase/pipeline/PPL_pipeline.sh -j /path/to/PPL-Toolbox.jar \
-g asm.fa \
-f reads.fq.gz \
-o PPL
```


# Estimating of the number of contig collapses based on HiFi data and popCNV
The popCNV_pipeline.sh script estimates the copy number of collapsed contigs collapse based on HiFi data using the popCNV software. The file used by popCNV for GPhase input is `collapse_num.txt` : popcnv/06.genes.round.cn. For details, see [popCNV](https://github.com/sc-zhang/popCNV)
```
/path/to/GPhase/pipeline/popCNV_pipeline.sh \
-f asm.fa \
-p output_prefix \
-t 32 \
-r reads.fq.gz
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
For more parameters, please refer to the `gphase pipeline -h`.
Below are some of the more important optional parameters:
- `--cluster_q` : the HiC/Pore-C data mapping quality score threshold (MAPQ) used during phasing is 1 by default. This applies when the input is a BAM file.
- `--scaffold_q` : the HiC/Pore-C data mapping quality score threshold (MAPQ) used during scaffolding is 0 by default. This applies when the input is a BAM file.
- `--hap_pm ` : the threshold for the intensity parameter of homologous sequence identification is 0.7 by default. If the heterozygosity of the assembled species is high, 0.6 can be used; if the heterozygosity of the species is low, 0.8 can be used.
- `-p` : the prefix for the output file should be specified as "only [a-zA-Z0-9.] allowed", meaning only the character ".", numbers, and uppercase and lowercase letters are allowed.

# Output file
GPhase will output a folder named gphase_output, which will generate the following four folders in sequence.
- `preprocessing` : Data preprocessing
- `cluster_chr` : Results of chromosome clustering
- `cluster_hap` : Haplotype clustering results within each chromosome
- `scaffold_hap` : Scaffolding results for each haplotype within each chromosome

# Final assembly result
The final assembly result file is located in the scaffold_hap folder and mainly contains the following:
- `gphase_final.agp` : unitig level assembly result agp file
- `gphase_final.fasta` : unitig level assembly result fasta file
- `gphase_final_rescue.agp` :  unitig level assembly result agp file after rescue
- `gphase_final_ctg2utg.txt` : correspondence between unitig and contig
- `gphase_final_contig.fasta` : Contig-level fasta sequence
- `gphase_final_contig.agp` : contig level assembly result agp file
- `gphase_final_contig_scaffold.fasta` : contig level assembly result fasta file

# Generate a Hi-C heatmap
Since the collapsed contig sequences have been duplicated, it is necessary to add markers to the collapsed unitigs in the AGP, fasta, and Hi-C mapping file to distinguish them when generating the Hi-C heatmap. There are two possible solutions:
### Rename fasta and agp + Remapping Hi-C + Generate Hi-C  heatmap

In short, the process first adds a fixed suffix to the duplicated collapsed unitigs (in both AGP and FASTA files), then remaps the Hi-C data using the mapQ:0 parameter (retaining multiple mappings), and finally generates a Hi-C heatmap using Juicer.

```
# rename agp and unitigs
python /Path/to/GPhase/scaffold_hap/rename_collapse_agp_pairs_fasta.py \
gphase_final.agp asm.fa rename --no-hic 

# chromap remapping
chromap -i -r rename.fa -o reindex
chromap --preset hic -x reindex -r rename.fa -q 0 \
    -1 HiC_1.fq.gz -2 HiC_2.fq.gz \
    --remove-pcr-duplicates -t 64 -o remap.chromap.pairs

# Generate Hi-C heatmap
bash /Path/to/GPhase/scaffold_hap/juicebox.sh \
-f rename.fa \
-a rename.agp \
-p remap.chromap.pairs \
-o final_hic -g /Path/to/GPhase

```
### Rename fasta,agp and pairs/bam + Generate Hi-C heatmap
In short, the process first adds a fixed suffix to the duplicated collapse unitigs (in AGP, FASTA, and pairs/bam files), and finally uses Juicer to generate the Hi-C heatmap.
```
# rename agp, unitigs and bam
python /Path/to/GPhase/scaffold_hap/rename_collapse_agp_pairs_fasta.py \
gphase_final.agp asm.fa rename --hic-file map.chromap.bam

# Generate Hi-C heatmap
bash /Path/to/GPhase/scaffold_hap/juicebox.sh \
-f rename.fa \
-a rename.agp \
-p rename.pairs \
-o final_hic -g /Path/to/GPhase
```



# Tips
1. The `cluster_q` and `scaffold_q` parameters are only enabled when the input mapping file format is BAM. If using pairs, the `mapQ` parameter of the mapping software (e.g., Chromap) can be adjusted, but it is not recommended to set `mapQ` to 0, as this will affect the accuracy of the phasing due to multiple-mapping.
2. When assembling `polyploids`, it is recommended to use `unitig-level` assembly `sequences` and `graph` for phasing assembly. Generally, unitig results in fewer errors compared to contig. Furthermore, using unitig allows for the utilization of more assembly graph information, leading to better assembly results.
3. GPhase can largely solve the problem of sequence collapse during assembly, but it cannot solve the problem of `large fragments collapsing` in haplotypes.


# Test dataset
To help you quickly verify the installation and use of the software, we provide a small test dataset. This dataset contains input data that demonstrates the core functionality of the software. You can download it from this link https://drive.google.com/drive/folders/1M_ZlSHBTDwtCHGrUI6uMCVutfIweECaY?usp=sharing

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