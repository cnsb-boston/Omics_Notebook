import subprocess
import os
import sys

# Get the location of the Notebook script
notebook_path = os.path.dirname(os.path.abspath(__file__))

# Make files paths to Pipeline scripts
pipeline_path = os.path.normpath(notebook_path + "/src" + "/Pipeline.py")
r_path = os.path.normpath(notebook_path + "/src" + "/Pipeline.R")

# Run pipeline scripts
p = subprocess.Popen(args=["python3", pipeline_path], stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
analysis_dir = str(p.communicate()[0], 'utf-8').rstrip()

subprocess.Popen(args=["Rscript", r_path, notebook_path, analysis_dir])
