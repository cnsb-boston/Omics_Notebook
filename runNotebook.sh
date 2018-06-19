#!/bin/bash -l

# 
# load require modules
module load java/1.8.0_92
module load java/1.8.0_151
module load gcc/5.3.0
module load R/3.4.3
module load python/3.6.2

# open cytoscape
/bin/sh "/project/cnsbomic/Tools/Cytoscape_v3.5.1/Cytoscape" -R 1234 &>/dev/null &

python "/project/cnsbomic/Omics_Notebook/Notebook.py"


