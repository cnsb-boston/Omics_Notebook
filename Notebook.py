import subprocess
import os

# Get the location of the Notebook script
notebook_path = os.path.dirname(os.path.abspath(__file__))

# Make files paths to Pipeline scripts
pipeline_path = os.path.normpath(notebook_path + "/src" + "/Pipeline.py")
r_path = os.path.normpath(notebook_path + "/src" + "/Pipeline.R")

# Run pipeline scripts
subprocess.run(["python", pipeline_path])

subprocess.run(["Rscript", r_path, notebook_path])

