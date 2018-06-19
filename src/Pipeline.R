# Get Notebook directory
run_directories <- commandArgs(trailingOnly=FALSE)

notebook_dir <- file.path(run_directories[6], "src");
analysis_dir <- file.path(run_directories[7]);

# source run variables
source(file.path(analysis_dir,"Variables.R"));

# set for BU SCC
if(BUSCC==TRUE) {
  .libPaths('/project/cnsbomic/Tools/R');
   Sys.setenv(RSTUDIO_PANDOC='/usr/local/apps/rstudio-0.98.1103/rstudio-0.98.1103/bin/pandoc'); 
} 

# load rmarkdown
library(rmarkdown);

# run notebook
rmarkdown::render(file.path(notebook_dir,'Notebook.Rmd'),
                  knit_root_dir=analysis_dir,
                  output_file=paste("Analysis_Notebook_",Sys.Date(),".html", sep=""),
                  output_dir=working_dir);

