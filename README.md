# SVHack_simulatemosaic

## Contributers

Meet SpikeVar and TykeVar:  

[<img src="images/Spike and Tyke image4.jpg" width="200"/>](workflow1.png)

## Background

In the context of individual genome comparison, mutations that appear within a small fraction of the population (<5%) are considered uncommon variants[[1]](#1). When assessing a population of cells from a tissue of the same individual in turn, uncommon variants only present in a small fraction of the cells are defined as a mosaic variants (MVs)[[2]](#2)[PMID: 32371413]. Recent studies have shown that there are potential disease association of for certain MVs[[3]](#3)[PMID: 32371413]. However, MVs are a challenging to detect because they are mixed in with data from the non-mutated cells and present in the same sequencing file. Therefore, several pipelines have been developed or adjusted to extract mosaic single nucleotide, structural or indel variants from whole genome sequencing data such as Sniffles[[4]](#4) [PMID: 29713083], DeepMosaic[[5](#5) [PMID: 36593400], Mutect2[[6]](#6) [PMID: 20644199], DeepVariant[[7]](#7) [PMID: 30247488]. To benchmark and validate the efficiency and accuracy of these methods, sequencing files with known MVs are necessary. We developed two simulation workflows called SpikeVar (Spike in Known Exogenous Variants) and TykeVar (Track in Your Key Endogenous Variants), which output sequencing read files with artificial MVs and a ground truth annotation file for the MVs. SpikeVar accomplishes this by spiking in real reads from a sample at user defined ratio into the sequencing file from a second sample. In contrast, TykeVar creates a list of random mutations and modifies a fraction of existing reads to match the user defined MV frequency.

## Method Description 

### 1. SpikeVar - Generation of sequencing data with a low frequencing of reads from another sample
[<img src="images/SpikeReads_flowchart.png" width="500"/>](workflow1.png)

Here we generate a workflow that can automatically spike in one sample at a given concentration (e.g. 5% at variant allele frequency, VAF) into another sample. We will demonstrate this over spiking in HG002 in to HG0733 for the purpose of demonstration. The downside is that the so generated mixed bam file will include 4 haplotype structures that we cannot correct. The challenging part is further that certain variants (e.g. HG002) will not be presented at the targeted VAF. For example, heterozygous variants wont be represented by 5% VAF but rather at ~2.5% VAF. To account for this we re-genotype variants and report only variants that should be identifiable at the user defined threshold or higher. 

## 2. TykeVar- Creation of sequencing data with a subset of modified reads
[<img src="images/TykeReads_flowchart.png" width="500"/>](Simulate_Mosaic_Simulation_on_reads_flowchart.png)

Here we modifiy reads directly at their reference position to include artifical mutations. In contrast to the above approach we dont introduce new haplotypes with this. The disadvantage however, is that more complex mutations (e.g. rearrangements, duplication or very long structural variants) will not be able to be introduced to the data itself, since the size of the reads is limited.

To get started, please refer to the [Tyke README](scripts/Tyke/README.md).

## Installation

## Example Implementation

## References

<a id="1">[1]</a>   
Sariya S, Lee JH, Mayeux R, Vardarajan BN, Reyes-Dumeyer D, Manly JJ, Brickman AM, Lantigua R, Medrano M, Jimenez-Velazquez IZ, Tosto G. Rare Variants Imputation in Admixed Populations: Comparison Across Reference Panels and Bioinformatics Tools. Front Genet. 2019;10:239. Epub 20190403. doi: 10.3389/fgene.2019.00239. PubMed PMID: 31001313; PMCID: PMC6456789.  
<a id="2">[2]</a>   
Miller CR, Lee K, Pfau RB, Reshmi SC, Corsmeier DJ, Hashimoto S, Dave-Wala A, Jayaraman V, Koboldt D, Matthews T, Mouhlas D, Stein M, McKinney A, Grossman T, Kelly BJ, White P, Magrini V, Wilson RK, Mardis ER, Cottrell CE. Disease-associated mosaic variation in clinical exome sequencing: a two-year pediatric tertiary care experience. Cold Spring Harb Mol Case Stud. 2020;6(3). Epub 20200612. doi: 10.1101/mcs.a005231. PubMed PMID: 32371413; PMCID: PMC7304353.


