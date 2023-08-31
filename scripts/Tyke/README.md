# Tyke

## Installation

This tool is written fully in python. To install the relevant dependencies, run
```
pip install -r $REPO_ROOT/requirements.txt
```
where `$REPO_ROOT` is the root folder of the repository.

## Generate simulated VCF

...

## Generate edited reads based on simulated VCF

```
python main.py -v <SIMULATED_VCF> -b <BAM> -r <REF> -o <OUTPUT_FASTQ>
```

This command above takes in the VCF which determines which variants to introduce into the reads.
The BAM file is used to find the reads which overlap with variant locations. Only a subset of the reads
corresponding to a particular variant location are edited. This is determined by the allele frequency.
The output FASTQ file has the edited reads. The query name of each read is kept the same.
