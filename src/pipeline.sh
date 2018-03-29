#!/bin/bash
# run pipeline as batch job
qsub -P cnsbomic -pe omp 16 -l gpus=0.125 -b y "/project/cnsbomic/Omics_Notebook/src/pipelineBatch.sh"