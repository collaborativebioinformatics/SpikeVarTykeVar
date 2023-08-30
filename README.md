# SVHack_simulatemosaic

The goal of this project is to simulate mosaic SNV and SV to evaluate existing tools. For this purpose we are implementing two simulation workflows. Mosaic variant are in generally defined as mutaitons that occure only in a tiny fraction of your reads. This can be for example in 5% of your reads which is also represented as 5% variant allele frequency (VAF).

![Workflow](workflow1.png)

## 1. Spike in experiments
Here we generate a workflow that can automatically spike in one sample at a given concentration (e.g. 5% VAF) into another sample. We will demonstrate this over spiking in HG002 in to HG0733 for the purpose of demonstration. The downside is that the so generated mixed bam file will include 4 haplotype structures that we cannot correct. The challaning part is further that certain variants (e.g. HG002) will not be presented at the targeted VAF. For example, heterozygous variants wont be represented by 5% bu rather at ~2.5%. To account for this we regenotype variants and report only variants that should be identifiable at the user defined threshold or higher. 

![Workflow2](Simulate_Mosaic_Simulation_on_reads_flowchart.png)

<img src="Simulate_Mosaic_Spike_in_framework_flowchart.png" width="50%" style="display: block; margin: auto;" />
<center>
`Simulate Mosaic - Spike in framework flowchart
</center>

## 2. Modification experiments
Here we are modifiying reads directly at their reference position to include artifical mutatioins. In contrast to the above approach we dont introduce new haplotypes with this. The disadvantage however, is that more complex mutations (e.g. rearrangements, duplication or very long SV) wont be able to introduce to the data itself. 
