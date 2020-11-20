import subprocess
import os
import sys

# Get options on how to run

# Get the location of the Notebook script
notebook_path = os.path.dirname(os.path.abspath(__file__))

# Make files paths to Pipeline scripts
pipeline_path = os.path.normpath(notebook_path + "/src" + "/Pipeline.py")
r_path = os.path.normpath(notebook_path + "/src" + "/Pipeline.R")

# Run pipeline scripts
p = subprocess.Popen(args=["python3", pipeline_path], stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
analysis_dir = str(p.communicate()[0], 'utf-8').rstrip()


# 1. Run natively
subprocess.Popen(args=["Rscript", r_path, notebook_path, analysis_dir])

# 2. Run with Docker
# subprocess.Popen(args=["docker run -it --rm -u docker -v ", notebook_path,":/home:rw -v ",analysis_dir,":/data:rw bblum/omics_notebook Rscript '/home/src/Pipeline.R' '/home' '/data'"])

# 3. Run with Singularity
# singularity_image = "/project/cnsbomic/Omics_Notebook/ON2.simg" # PATH/TO/SINGULARITY IMAGE
# subprocess.Popen(args=["singularity run --bind ",notebook_path,":'/home':rw --bind ",analysis_dir,":'/data':rw  ",singularity_image," Rscript '/home/src/Pipeline.R' '/home' '/data'"])
