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
