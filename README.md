# A-multi-omicsRuan
Scripts used for data analysis in "A multi-omics integrative analysis based on CRISPR screens re-defines the pluripotency regulatory network in ESCs"

We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.

Most scripts rely on the use of basic read aligment and processing tools:

- FastQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- MultiQC https://multiqc.info/
- Trim_galore https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
- Samtools http://samtools.sourceforge.net/
- BEDtools https://bedtools.readthedocs.io/en/latest/
- Hisat2 https://daehwankimlab.github.io/hisat2/
- Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/
- DeepTools https://deeptools.readthedocs.io/en/develop/index.html
- R https://www.r-project.org/
- MAGeCK https://sourceforge.net/projects/mageck/



## WGCNA

WGCNA.R

## MAGeCK

Brie library+non-target(mageck).csv

MAGeCK.sh

## Pipelines for ChIP-seq and RNA-seq

RNAseq_Alignment_Pipeline.txt

ChIPseq_Alignment_Pipeline.txt

## **Calculate Z-Score**

Calculate Z-Score.md

## **Identification of the histone modification status**

bamCoverage&Deeptool.sh

## **Adasamping predicts target genes**

Adasamping predicts target genes.md

## **Analysis of gene set enrichment pattern**

Analysis in GBMLGG.R

## Module expression activity analysis & Survival analysis

Analysis in GBMLGG.sh




