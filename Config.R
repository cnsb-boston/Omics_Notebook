# notebook_path is the path to the Omics_Notebook installation when this file is sourced by Notebook.R
libdir="/project/cnsb2022/RTools/ON4.1" # extra libraries to map into the container
singularity_img=paste0(notebook_path,"/ON4.simg") # ignored if singularity is not used or image is provided on the command line
docker_img="cnsbboston/omicsnotebook" # ignored if docker is not used or image is provided on the command line
startdir="/projectnb/cnsb2022" # default starting directory for the GUI file selection
env_vars=c("R_LIBS=/usr/local/lib/R/local-library") # environment variables to use during the analysis
