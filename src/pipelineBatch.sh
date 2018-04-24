#!/bin/bash -l

# run pipeline through batch (non-interactive)

# load require modules
module load java/1.8.0_92
module load java/1.8.0_151
module load gcc/5.3.0
module load R/3.4.3

# get the pipeline directory
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

# run R pipeline and pass directory
Rscript $SCRIPTPATH/src/Pipeline.R $SCRIPTPATH
