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
run_directories <- commandArgs(trailingOnly=FALSE)

index <- length(run_directories)

notebook_dir <- file.path(run_directories[(index-1)], "src");
analysis_dir <- file.path(run_directories[index]);


# source run variables
setwd(analysis_dir);
source(file.path(analysis_dir, "Parameters.R"));

# Set local file paths, if not needed, set BUSCC to false
if(inherit_paths==TRUE) {
  .libPaths(libraries_path);
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
                  output_file=paste(gsub("\\.","",make.names(project_name)),"_",gsub("-","",Sys.Date()),".html", sep=""),
                  output_dir=analysis_dir);

