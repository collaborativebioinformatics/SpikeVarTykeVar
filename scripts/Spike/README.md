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