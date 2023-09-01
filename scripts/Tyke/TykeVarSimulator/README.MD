The VCF simulator generates a random set of mosaic variants (SNVs and SVs). The variants can be parameterized with VAF, 
number of variants to simulate and the size of the variations. The generated file is in the Sniffles VCF format.

The variants generated here act as an input into the read editor step (described below) which generates modified reads 
with the variants inserted into them. The same VCF file is also the ground truth for evaluating mosaic variant callers.

python vcfgen.py <path_to_bam> <path_to_ref> <output_path_prefix>
