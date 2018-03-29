#!/bin/bash -l

# run pipeline

# load require modules
module load java/1.8.0_92
module load java/1.8.0_151
module load gcc/5.3.0
module load R/3.4.3

# open cytoscape
/bin/sh "/project/cnsbomic/Tools/Cytoscape_v3.5.1/Cytoscape" -R 1234 &

# get the pipeline directory
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

# run R pipeline and pass directory
Rscript $SCRIPTPATH/Pipeline.R $SCRIPTPATH

#close cytoscape and end interactive session
kill %1
exit
