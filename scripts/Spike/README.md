Installing two conda environments with python 3.10

```bash
conda create -n mosaicSim python=3.10
conda activate mosaicSima
pip install -r requirements.txt
```

To spike-in sample B into sample A:

Dependencies:

* mosdepth 0.3.2
* samtools 1.15.1
* Python 3.6.8
* bcftools

Example:
```
sh spike-in.sh <path to sampleA.bam> <path to sampleB.bam> <spike-in ratio> <path to samtools binary> <path to mosdepth binary> <output dirpath> <path to script calculate_ratio.py>
```


### Datasets

##### Long Reads

**HG002**

BAM

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam

VCF

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz


**HG00733**

BAM

In house. File name: NHGRI_UCSC_panel-HG00733-grch37_sorted.bam

VCF

In house. File name: HG00733_truvari_li_grch37_SV_with_IDs.vcf

##### Short Reads

**HG002**

BAM

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam

**HG007**

BAM 

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG007_NA24695-hu38168_mother/CompleteGenomics_normal_cellsDNA/mom_NA24695_GS000037477-ASM/BAM/chr20_mapping_sorted_header.bam

VCF

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG007_NA24695_mother/latest/GRCh37/HG007_GRCh37_1_22_v4.2.1_benchmark.vcf.gz