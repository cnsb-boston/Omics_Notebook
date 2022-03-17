############################################
##
##  SET FILE PATHS FOR LOCAL INSTALL
##
############################################

## R libraries path
#libraries_path = '/project/cnsbomic/Tools/R'

## Pandoc path
#pandoc_path = 

## Inherit paths
#inherit_paths = "TRUE"

# Loaded from Parameters

############################################

# Get Notebook directory - if running pipeline
if(length(commandArgs(trailingOnly=TRUE))>0) args <- commandArgs(trailingOnly=TRUE)

notebook_dir <- file.path(args[1], "src");
analysis_dir <- file.path(args[2]);
param_file=file.path(analysis_dir,if(length(args)>2)args[3] else "Parameters.R")

# source run variables
setwd(analysis_dir);
source(param_file);

if(analysis_dir != working_dir){
  write(paste("annotation_filename <- '",file.path(analysis_dir,basename(annotation_filename)), "';", sep=""),
        file = param_file, append = T)
  write(paste("working_dir <- '", analysis_dir, "';", sep=""),
        file = param_file, append = T)
}


# Set local file paths, if not needed, set BUSCC to false
if(inherit_paths==TRUE) {
  #.libPaths(libraries_path);
   #Sys.setenv(RSTUDIO_PANDOC=pandoc_path);
}

# load rmarkdown
library(rmarkdown);

# run normalization Shiny
if(shinyNorm==TRUE){
  library(shiny);
  source(file.path(notebook_dir, "Pipeline_Norm.R"));
}

# run notebook
rmarkdown::render(file.path(notebook_dir,'Notebook.Rmd'),
                  envir=new.env(),
                  knit_root_dir=analysis_dir,
                  intermediates_dir=analysis_dir,
                  params=list(param_file=param_file),
                  output_file=paste(gsub("\\.","",make.names(project_name)),"_",gsub("-","",Sys.Date()),".html", sep=""),
                  output_dir=analysis_dir);

