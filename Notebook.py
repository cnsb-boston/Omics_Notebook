import subprocess
import time

print("!!!!!!WARNING!!!!!!\nYou may be using an outdated Notebook.sh\n This analysis should still run, but please take a fresh copy of the script from /project/cnsbomic/Omics_Notebook/Notebook.sh\n")
time.sleep(5)
subprocess.call("/project/cnsbomic/Omics_Notebook/Notebook.sh", shell=True)
