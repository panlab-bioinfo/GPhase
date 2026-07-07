# GPhase: A Phasing Assembly Tool Leveraging an Assembly Graph and Hi-C/Pore-C/Omni-C Data

GPhase leverages an assembly graph and Hi-C/Pore-C data to facilitate genome assembly phasing, automatically resolves and assigns collapsed sequences, and fills assembly gaps based on the graph structure.
---

## Table of Contents

- [Installation](#installation)
- [Step1: Mapping Hi-C/Pore-C data to assembly](#step1-mapping-hi-cpore-c-data-to-assembly)
- [Step2: Estimating the number of contig collapses based on HiFi data and popCNV](#step2-estimating-the-number-of-contig-collapses-based-on-hifi-data-and-popcnv)
- [Step3: Running the GPhase scaffolding pipeline](#step3-running-the-gphase-scaffolding-pipeline)
- [Output directory structure](#output-directory-structure)
- [Final assembly files](#final-assembly-files)
- [Visualization and curation in Juicebox](#visualization-and-curation-in-juicebox)
- [Tips](#tips)
- [Test dataset](#test-dataset)
- [Contact](#contact)

## Installation

To install GPhase, follow these steps:
```
# conda
git clone --depth 1 https://github.com/panlab-bioinfo/GPhase.git
cd GPhase
conda env create -f gphase_environment.yml
conda activate gphase
./gphase -h

# docker
docker pull --platform linux/amd64 tanging1024/gphase:latest
docker run --platform linux/amd64 -v /your/data/path:/your/data/path -w /your/data/path tanging1024/gphase:latest gphase -h

# singularity
singularity pull gphase.sif docker://tanging1024/gphase:latest
singularity exec --bind /your/data/path:/your/data/path gphase.sif gphase

```

> [!WARNING]
>
> GPhase requires the raw, unprocessed hifiasm unitig FASTA (`*.p_utg.fa`) as the assembly input. This FASTA must correspond exactly to the hifiasm primary unitig GFA (`*.p_utg.gfa`), with matching sequence IDs and graph records. Do not rename, reorder, filter, polish, purge, scaffold, or otherwise modify the primary unitigs before running GPhase. The target genome should also have sufficient heterozygosity; if the genome is nearly homozygous, haplotype phasing will have little biological meaning.

## Step1: Mapping Hi-C/Pore-C data to assembly

GPhase supports multiple long-range data types, including Hi-C, Pore-C and Omni-C. It also supports mapping files in pairs (PA5) and BAM formats.

1. For Hi-C reads, you can use Chromap or other mapping tools, such as BWA. When using Chromap, if the default MAPQ parameters do not produce satisfactory results, the `--MAPQ-threshold` value can be lowered to include more Hi-C mapping information. When using other mapping software, the BAM files need to be sorted.
```
chromap -i -r p_utg.fa -o index
chromap --preset hic -x index -r p_utg.fa -q 0 \
    -1 HiC_1.fq.gz -2 HiC_2.fq.gz \
    --remove-pcr-duplicates -t 64 --SAM -o map.chromap.sam
samtools view -@ 64 -bh map.chromap.sam -o map.chromap.bam
```
2. For contact-pair long reads, including Pore-C and CiFi, we recommend using `gphase contact-pair` as the default workflow. This command performs read mapping internally with minimap2, then converts the alignments to a Hi-C-like BAM file `map.concatemer2pe.bam` with `concatemer2pe.py`. That BAM can be passed directly to `gphase pipeline -m`. Use `-x map-ont` for Pore-C reads and `-x map-hifi` for CiFi reads. The `-o` option specifies the output directory prefix, and the final BAM path is `<prefix>/map.concatemer2pe.bam`. For all `contact-pair` parameters, see [doc/README.md#gphase-contact-pair](doc/README.md#gphase-contact-pair).
```
/path/to/GPhase/gphase contact-pair \
    p_utg.fa \
    reads1.fq.gz reads2.fq.gz \
    -x map-ont \
    -o contact_pair \
    -t 32
```

The previous Pore-C workflow based on [PPL Toolbox](https://github.com/versarchey/PPL-Toolbox) is still available as a backup workflow. If you need to reproduce the results reported in the paper, you can use this PPL-based process to generate the final pairs file `map.PPL.pairs` and then input it into GPhase. For all `ppl` parameters, see [doc/README.md#gphase-ppl](doc/README.md#gphase-ppl).
```
/path/to/GPhase/gphase ppl \
    -g p_utg.fa \
    -f reads.fq.gz \
    -o PPL
```

## Step2: Estimating the number of contig collapses based on HiFi data and popCNV

The `gphase popcnv` command estimates the copy number of collapsed contigs based on HiFi data using the popCNV software. The file used by popCNV for GPhase input is `collapse_num.txt` : `popcnv/06.genes.round.cn`. For details, see [popCNV](https://github.com/sc-zhang/popCNV). For all `popcnv` parameters, see [doc/README.md#gphase-popcnv](doc/README.md#gphase-popcnv).
```
/path/to/GPhase/gphase popcnv \
    -f p_utg.fa \
    -p output_prefix \
    -t 32 \
    -r reads.fq.gz
```

## Step3: Running the GPhase scaffolding pipeline

```
/path/to/GPhase/gphase pipeline \
    -f p_utg.fa \
    -g p_utg.gfa \
    -c collapse_num.txt \
    -m map.chromap.bam \
    --n_chr 12 \
    --n_hap 4 \
    -p output_prefix \
    --rescue \
    --min_len 50
```

Required parameters:
- `-f` : Genome assembly file in FASTA format (unitigs).
- `-g` : Assembly graph file in GFA format.
- `-c` : File with contig collapse information (from popCNV: `popcnv/06.genes.round.cn`).
- `-m` : Hi-C/Pore-C/Omni-C/CiFi mapping file in `.bam` or `.pairs` format. For contact-pair data, the recommended input is `contact_pair/map.concatemer2pe.bam`, where `contact_pair` is the `-o` output directory/prefix used by `gphase contact-pair`; for paper reproduction, `PPL/map.PPL.pairs` can also be used.
- `--n_chr` : Number of chromosomes.
- `--n_hap` : Number of haplotypes.
- `-p` : Prefix for output files. Only the character `.`, numbers, and uppercase/lowercase letters are allowed (`[a-zA-Z0-9.]`).

Below are some of the more important optional parameters:
- `--cluster_q` : Hi-C/Pore-C/Omni-C/CiFi mapping quality score threshold (MAPQ) used during clustering. The default is `1`. This applies when the input is a BAM file.
- `--scaffold_q` : Hi-C/Pore-C/Omni-C/CiFi mapping quality score threshold (MAPQ) used during scaffolding. The default is `0`. This applies when the input is a BAM file.
- `--hap_pm` : The threshold for the intensity parameter of homologous sequence identification. The default is `0.7`. If the heterozygosity of the assembled species is high, `0.6` can be used; if the heterozygosity of the species is low, `0.8` can be used.
- `--chr_pm` : Similarity threshold for chromosome-level partig clustering. The default is `0.95`.
- `--nor_hic` : Normalization mode for 3C link connections. Choices are `no`, `ratio`, and `length`; default is `ratio`.
- `--min_len` : Minimum scaffold length in kb in HapHiC sorting. The default is `50`.
- `--thread` : Number of parallel processes used in scaffolding. The default is `12`.

For more parameters, please refer to `gphase pipeline -h` or [doc/README.md#gphase-pipeline](doc/README.md#gphase-pipeline).

## Output directory structure

GPhase will create a folder named `gphase_output` containing the following four subfolders:
- `preprocessing` : Data preprocessing.
- `cluster_chr` : Results of chromosome clustering.
- `cluster_hap` : Haplotype clustering results within each chromosome.
- `scaffold_hap` : Scaffolding results for each haplotype within each chromosome.

## Final assembly files

The final assembly result files are located in the `scaffold_hap` folder. GPhase produces chromosome-level assemblies where each scaffold is composed of contiguous sequences. There are two types of results based on the granularity of the contiguous sequences:

### Results based on unitigs

Unitigs are unbranched paths in the assembly graph. They are typically short but correctly phased. The following files contain results where adjacent unitigs remain separate:
- `gphase_final.unitig.fasta` : Unitig sequences (derived from `p_utg.fa`; collapsed unitigs are duplicated and renamed) in FASTA format.
- `gphase_final.unitig.scaffold.agp` : Chromosome-level scaffold structure of unitigs in AGP format.
- `gphase_final.unitig.scaffold.fasta` : Chromosome-level scaffold sequences assembled from unitigs in FASTA format.

### Results based on contigs

Contigs are formed by stitching adjacent unitigs after phasing and scaffolding — equivalent to gap closing along the correct path identified from branched paths in the assembly graph. The following files contain results where adjacent unitigs have been merged into contigs:
- `gphase_final.contig.fasta` : Contig sequences (obtained by merging adjacent unitigs) in FASTA format.
- `gphase_final.contig.scaffold.agp` : Chromosome-level scaffold structure of contigs in AGP format.
- `gphase_final.contig.scaffold.fasta` : Chromosome-level scaffold sequences assembled from contigs in FASTA format.

The correspondence between unitigs and contigs can be found in `gphase_final_ctg2utg.txt`.

## Visualization and curation in Juicebox

The input for this workflow can be either the unitig-based results (`gphase_final.unitig.fasta` and `gphase_final.unitig.scaffold.agp`) or the contig-based results (`gphase_final.contig.fasta` and `gphase_final.contig.scaffold.agp`), depending on your needs. Unitigs offer lower contiguity but rarely contain phasing errors, making them simpler to curate in Juicebox. Contigs provide higher contiguity, but if phasing errors are present, correcting them first requires breaking contigs at misjoined sites, which is difficult to perform accurately.

### Remapping long-range data

In the following instructions, `<prefix>.fasta` and `<prefix>.scaffold.agp` refer to the input files chosen above (e.g., `gphase_final.unitig` or `gphase_final.contig`).

Before visualizing in Juicebox, the long-range data must be remapped to the selected sequences. The method depends on the data type:

**For Hi-C and Omni-C reads**, use Chromap to produce a pairs file:

```
# Build index
chromap -i -r <prefix>.fasta -o reindex

# Remap Hi-C/Omni-C reads
chromap --preset hic -x reindex -r <prefix>.fasta -q 1 \
    -1 HiC_1.fq.gz -2 HiC_2.fq.gz \
    --remove-pcr-duplicates -t 64 -o remap.chromap.pairs
```

**For Pore-C and CiFi reads**, use `gphase contact-pair`, which internally maps reads with minimap2 and converts the alignments to a Hi-C-like BAM file. Use `-x map-ont` for Pore-C reads and `-x map-hifi` for CiFi reads. The output BAM is placed at `remapping/map.concatemer2pe.bam`:

```
# Remap Pore-C/CiFi reads
/path/to/GPhase/gphase contact-pair \
    <prefix>.fasta \
    reads1.fq.gz reads2.fq.gz \
    -x map-ont \
    -q 2 \
    -o remapping \
    -t 32
```

> [!NOTE]
>
> During visualization, Hi-C/Omni-C alignments are filtered at MAPQ ≥ 1 by default (`-q 1`), and Pore-C/CiFi alignments at MAPQ ≥ 2 (`-q 2`), to exclude unreliable multi-mapped reads. If the resulting contact heatmap appears too sparse, you may lower these thresholds. Note that lower thresholds can introduce a large number of spurious contact signals.

### Generating .assembly and .hic files

Use `juicebox.sh` to generate `.assembly` and `.hic` files from the assembly and the remapped long-range data. Replace `<mapfile>` with `remap.chromap.pairs` for Hi-C/Omni-C, or `remapping/map.concatemer2pe.bam` for Pore-C/CiFi:

```
bash /Path/to/GPhase/scaffold_hap/juicebox.sh \
    -f <prefix>.fasta \
    -a <prefix>.scaffold.agp \
    -p <mapfile> \
    -o juicebox \
    -g /Path/to/GPhase
```

### Curating in Juicebox

1. Open Juicebox Assembly Tools (JBAT).
2. Load `juicebox.hic` and `juicebox.assembly`.
3. Manually correct misassemblies as needed.
4. Export the modified results as `juicebox.review.assembly`.

### Regenerating the curated assembly

After manual curation, use `juicer post` to regenerate the scaffold AGP and FASTA files by rearranging `<prefix>.fasta` according to the chromosome assignment, order and orientation specified in `juicebox.review.assembly`:

```
/path/to/GPhase/src/HapHiC/utils/juicer post \
    -o reviewed \
    juicebox.review.assembly \
    juicebox.liftover.agp \
    <prefix>.fasta
```

- `reviewed.FINAL.agp` : Scaffold structure after manual curation in AGP format.
- `reviewed.FINAL.fa` : Scaffold sequences after manual curation in FASTA format.

## Tips

1. The `cluster_q` and `scaffold_q` parameters are only enabled when the input mapping file format is BAM. If using pairs, the `mapQ` parameter of the mapping software (e.g., Chromap) can be adjusted, but it is not recommended to set `mapQ` to 0, as this will affect the accuracy of the phasing due to multiple mapping.
2. When assembling `polyploids`, it is recommended to use `unitig-level` assembly `sequences` and `graph` for phasing assembly. Generally, unitigs results in fewer errors compared to contig. Furthermore, using unitigs incorporates more assembly graph information, leading to better assembly results.
3. GPhase can largely resolve sequence collapse during assembly, but cannot resolve the large fragment collapse in haplotypes.

## Test dataset

To help you quickly verify the installation and use of the software, we provide a small test dataset. This dataset contains input data that demonstrates the core functionality of the software. You can download it from this link https://drive.google.com/drive/folders/1M_ZlSHBTDwtCHGrUI6uMCVutfIweECaY?usp=sharing

Use the following command to run the test dataset:
```
tar -zxvf test_dataset.tar.gz
export PATH=$PATH:/path/to/GPhase
bash run_gphase.sh
```

## Contact

This software is developed by Professor Wei-Hua Pan's team at the Shenzhen Institute of Genome Research, Chinese Academy of Agricultural Sciences.

If you have any questions or concerns while using the software, please submit an issue in the repository or contact us through the following methods:
### Email:
#### Prof. Pan: panweihua@caas.cn
#### Du Wenjie: duwenjie1024@163.com
