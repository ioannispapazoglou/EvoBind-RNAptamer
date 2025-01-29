# EvoBind-RNAptamer

An evolution-algorithm-based platform that simulates natural selection in silico to explore the vast RNA sequence space and identify optimal aptamers that bind to given protein targets.

## INSTALLATION

Third-party dependencies installation:

1. RhoFold+  :  Download the official git repo inside ./strpred/ directory. Follow the instructions in ./strpred/setup4nsga2.txt
  
2. HDock     :  Download the pre-compiled executables and place in ./3pbins/ directory. Follow the instructions in ./3pbins/download3rdparties.txt
  
3. RASP      :  Download the pre-compiled executables and place in ./3pbins/ directory. Follow the instructions in ./3pbins/download3rdparties.txt

All scripts should run under the provided conda environment:

conda env create -f evobind.yml

Another option is to run the scripts under the RhoFold+ conda environmnet (after installing the extra pymoo dependency), given in: https://github.com/ml4bio/RhoFold

Usage instructions and installation steps for Miniconda / Anaconda refer to the official web page: https://docs.anaconda.com/free/miniconda/index.html

## Inference

python3 inference.py --proteintarget ./test/ang2.pdb --outputfolder ./output/ 

To experiment abroad the default parameters: python3 inference.py --help

Easy jupyter notebook inference is also available.

#### Notes:

1.	Is adviced for protein target structures to be prepared (protonated / fully modelled). Proteins must contain OXT atom, otherwise will raise fatal error in complex minimization.
2.	Currenty the algorithm supports only solvable proteins, as membranes and other system components tend to raise errors or inaccuracies to the results.

## DurdagiLab message

Thank you for your interest in EvoBind-RNAptamer. We are constantly working on improving our code. 
