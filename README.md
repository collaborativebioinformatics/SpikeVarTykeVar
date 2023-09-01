# SVHack_simulatemosaic

## Contributers

Meet SpikeVar and TykeVar:  

[<img src="images/Spike and Tyke image4.jpg" width="200"/>](workflow1.png)

## Background

In the context of individual genome comparison, mutations that appear within a small fraction of the population (<5%) are considered uncommon variants[<sup>1</sup>](#1). When assessing a population of cells from a tissue of the same individual in turn, uncommon variants only present in a small fraction of the cells are defined as a mosaic variants (MVs)[<sup>2</sup>](#2). Recent studies have shown that there are potential disease association of for certain MVs[<sup>2</sup>](#2). However, MVs are a challenging to detect because they are mixed in with data from the non-mutated cells and present in the same sequencing file. Therefore, several pipelines have been developed or adjusted to extract mosaic single nucleotide, structural or indel variants from whole genome sequencing data such as Sniffles[<sup>3</sup>](#3), DeepMosaic[<sup>4</sup>](#4), Mutect2[<sup>5</sup>](#5), DeepVariant[<sup>6</sup>](#6). To benchmark and validate the efficiency and accuracy of these methods, sequencing files with known MVs are necessary. We developed two simulation workflows called SpikeVar (Spike in Known Exogenous Variants) and TykeVar (Track in Your Key Endogenous Variants), which output sequencing read files with artificial MVs and a ground truth annotation file for the MVs. SpikeVar accomplishes this by spiking in real reads from a sample at user defined ratio into the sequencing file from a second sample. In contrast, TykeVar creates a list of random mutations and modifies a fraction of existing reads to match the user defined MV frequency.

## Method Description 

### 1. SpikeVar - Generation of sequencing data with a low frequencing of reads from another sample
[<img src="images/SpikeVar_flowchart_updated.png" width="500"/>](workflow1.png)

The SpikeVar workflow outputs a mixed sequencing read dataset in .bam format containing reads from one dominant sample and reads from another sample spiked in at a user defined ratio corresponding to the simulated mosaic variant allele frequency (VAF) together with a .vcf file annotating the confirmed mosaic variant locations within the mixed dataset. The SpikeVarDatasetCreator takes aligned sequencing reads from sample A and sample B as the initial input. In this step, a spike-in methodology is applied to strategically introduce x% of mutations from one sample to another using <insert tool>. Accordingly, sample A is first down-sampled to retain 100-x% of its original reads, then sample B is down-sampled to x% considering the coverage differences between the samples. Using <insert tool>, both down-sampled datasets are then merged to create a mixed dataset that represents a sequence read dataset with mosaic variants, including structural variations (SVs), single nucleotide variations (SNVs), and insertions/deletions (indels).  

The SpikeVarReporter then determines VAFs for each variant in the mixed dataset using <insert tool> based on the mixed variant locations derived by merging the .vcf files from sample A and sample B using <insert tool>. Variants with VAFs exceeding or equal to the introduced mutations (i.e., x%) are then selected to create a truth set for benchmarking using <insert tool>.  
 
To assess a mosaic variant caller’s sensitivity and accuracy, the same mixed dataset is used to call mosaic variants. The output mosaic variant locations and VAFs are then compared to the truth set for validation.  

## 2. TykeVar- Creation of sequencing data with a subset of modified reads
[<img src="images/TykeVar_flowchart_updated.png" width="500"/>](Simulate_Mosaic_Simulation_on_reads_flowchart.png)


To get started, please refer to the [Tyke README](scripts/Tyke/README.md).

## Installation

## Example Implementation

### SpikeVar

Here, we use the SpikeVar workflow to automatically spike in sample HG002 at a 5% concentration into sample HG0733, to result in a 5% mosaic variant allele frequency (VAF). A downside is that the generated mixed .bam file will include 4 haplotype structures which cannot be corrected for. Furthermore, certain variants (e.g. HG002 variants) will not be presented at the targeted VAF. Forexample, heterozygous variants will not be represented by 5% VAF but rather at ~2.5% VAF. To account for this we re-genotype variants and report only variants that should be identifiable at the user-defined threshold or higher VAF.   

### TykeVar

Here, we use the TykeVar workflow to modifiy reads of HG002 directly at their reference position by including artifical mutations to represent at variant allele frequency of 5%. In contrast to the above approach we do not introduce new haplotypes with this. However, more complex mutations (e.g. rearrangements, duplication or very long structural variants) will not be able to be introduced to the data itself, since the size of the reads is limited.

## References

<a id="1">[1]</a> Sariya S, Lee JH, Mayeux R, Vardarajan BN, Reyes-Dumeyer D, Manly JJ, Brickman AM, Lantigua R, Medrano M, Jimenez-Velazquez IZ, Tosto G. Rare Variants Imputation in Admixed Populations: Comparison Across Reference Panels and Bioinformatics Tools. Front Genet. 2019;10:239. Epub 20190403. doi: 10.3389/fgene.2019.00239. PubMed PMID: 31001313; PMCID: PMC6456789.  

<a id="2">[2]</a> Miller CR, Lee K, Pfau RB, Reshmi SC, Corsmeier DJ, Hashimoto S, Dave-Wala A, Jayaraman V, Koboldt D, Matthews T, Mouhlas D, Stein M, McKinney A, Grossman T, Kelly BJ, White P, Magrini V, Wilson RK, Mardis ER, Cottrell CE. Disease-associated mosaic variation in clinical exome sequencing: a two-year pediatric tertiary care experience. Cold Spring Harb Mol Case Stud. 2020;6(3). Epub 20200612. doi: 10.1101/mcs.a005231. PubMed PMID: 32371413; PMCID: PMC7304353.  

<a id="3">[3]</a> Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A, Schatz MC. Accurate detection of complex structural variations using single-molecule sequencing. Nat Methods. 2018;15(6):461-8. Epub 20180430. doi: 10.1038/s41592-018-0001-7. PubMed PMID: 29713083; PMCID: PMC5990442.  

<a id="4">[4]</a> Yang X, Xu X, Breuss MW, Antaki D, Ball LL, Chung C, Shen J, Li C, George RD, Wang Y, Bae T, Cheng Y, Abyzov A, Wei L, Alexandrov LB, Sebat JL, Network NBSM, Gleeson JG. Control-independent mosaic single nucleotide variant detection with DeepMosaic. Nat Biotechnol. 2023;41(6):870-7. Epub 20230102. doi: 10.1038/s41587-022-01559-w. PubMed PMID: 36593400; PMCID: PMC10314968.  

<a id="5">[5]</a> McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297-303. Epub 20100719. doi: 10.1101/gr.107524.110. PubMed PMID: 20644199; PMCID: PMC2928508.  

<a id="6">[6]</a> Poplin R, Chang PC, Alexander D, Schwartz S, Colthurst T, Ku A, Newburger D, Dijamco J, Nguyen N, Afshar PT, Gross SS, Dorfman L, McLean CY, DePristo MA. A universal SNP and small-indel variant caller using deep neural networks. Nat Biotechnol. 2018;36(10):983-7. Epub 20180924. doi: 10.1038/nbt.4235. PubMed PMID: 30247488.




