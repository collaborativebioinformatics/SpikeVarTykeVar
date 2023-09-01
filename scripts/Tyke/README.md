# TykeVAR

## Installation

This tool is written fully in python. To install the relevant dependencies, run
```
pip install -r $REPO_ROOT/requirements.txt
```
where `$REPO_ROOT` is the root folder of the repository.

## How to use it

### TykeVar

#### Generate simulated VCF

The VCF simulator generates a random set of mosaic variants (SNVs and SVs). The variants
can be parameterized with VAF, number of variants to simulate and the size of the variations.
The generated file is in the Sniffles VCF format.

The variants generated here act as an input into the read editor step (described below) which
generates modified reads with the variants inserted into them. The same VCF file is also the
ground truth for evaluating mosaic variant callers.

```
python vcfgen.py <path_to_bam> <path_to_ref> <output_path_prefix>

e.g. python vcfgen.py chr22.bam hs37d5.fa chr22
The above generates a chr22SV.vcf and chr22SNV.vcf file
```

#### Generate edited reads based on simulated VCF

```
python main.py -v <SIMULATED_VCF> -b <BAM> -r <REF> -o <OUTPUT_FASTQ>
```

This command above takes in the VCF which determines which variants to introduce into the reads.
The BAM file is used to find the reads which overlap with variant locations. Only a subset of the reads
corresponding to a particular variant location are edited. This is determined by the allele frequency.
The output FASTQ file has the edited reads. The query name of each read is kept the same.

## Example Implementation

### TykeVar

Here, we use the TykeVar workflow to modifiy reads of HG002 directly at their reference position by including artifical mutations to represent at variant allele frequency of 5%. In contrast to the above approach we do not introduce new haplotypes with this. However, more complex mutations (e.g. rearrangements, duplication or very long structural variants) will not be able to be introduced to the data itself, since the size of the reads is limited.


#### Fetch data
In order to simulate and edit reads, the pipeline first needs an initial set of aligned reads and a reference. For our demonstration, we will use the GIAB datasets.

Reads - `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.2.4_2020-01-22/HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam` and `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.2.4_2020-01-22/HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam.bai`

Reference - `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/hs37d5.fa.gz`

#### Generate variants and modified reads

First we decompress the FASTA file.
```
gunzip hs37d5.fa.gz -c hs37d5.fa
```

Then we simulate variants
```
python vcfgen.py HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam hs37d5.fa hg002
```

Generate a set of modified reads with inserted variants.
```
python main.py -v hg002SV.vcf -b HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam -r hs37d5.fa -o hg002_modified_reads.fastq
```

#### Re-align modified reads and merge them
Once the new reads are generated, they need to be re-aligned and re-inserted back into
the dataset by replacing the original reads.

TODO:

#### Run your favorite mosaic variant caller

TODO:
