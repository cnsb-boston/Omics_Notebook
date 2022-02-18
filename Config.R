# notebook_path is the path to the Omics_Notebook installation when this file is sourced by Notebook.R
libdir="/project/cnsbomic/RTools/ON4.1" # extra libraries to map into the container
singularity_img=paste0(notebook_path,"/ON4.simg") # ignored if sigularity is not used or image is provided on the command line
startdir="/projectnb/cnsbomic" # default starting directory for the GUI file selection
