The VCF simulator generates a random set of mosaic variants (SNVs and SVs). The variants can be parameterized with VAF, 
number of variants to simulate and the size of the variations. The generated file is in the Sniffles VCF format.

The variants generated here (samples included for SV and SNV) act as an input into the read editor step which generates modified reads 
with the variants inserted into them. The same VCF file is also the ground truth for evaluating mosaic variant callers.

python TykeVarSimulator.py <path_to_bam> <path_to_ref> <output_path_prefix>

To modify runtime parameters, go to header of the script and change variables according to included comments.
