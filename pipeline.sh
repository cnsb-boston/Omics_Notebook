#!/bin/bash -l

# run pipeline
module load java/1.8.0_92
module load java/1.8.0_151
module load gcc/5.3.0
module load R/3.4.3
/bin/sh "/project/cnsbomic/Tools/Cytoscape_v3.5.1/Cytoscape" -R 1234 &
Rscript /project/cnsbomic/Omics_Notebook/Pipeline.R
kill %1
exit
