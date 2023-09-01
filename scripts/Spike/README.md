Installing two conda environments with python 2 and 3

Python 2.7 is needed for svtyper
```bash
conda create -n mosaicSim27 python=2.7
conda activate mosaicSim
pip install git+https://github.com/hall-lab/svtyper.git
```
Python 3.10 is needed for the rest of the pipeline
```bash
conda create -n mosaicSim python=3.10
conda activate mosaicSim
pip install -r requirements.txt
```

To spike-in sample B into sample A:
```bash
sh spike-in.sh <path to sampleA.bam> <path to sampleB.bam> <spike-in ratio> <path to samtools binary> <path to mosdepth binary> <output dirpath> <path to script calculate_ratio.py>
```
