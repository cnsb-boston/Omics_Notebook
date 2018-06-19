# Get Notebook directory
notebook_dir_sys <- tail(commandArgs(trailingOnly=TRUE), n=1);
notebook_dir_sys <- file.path(notebook_dir_sys, "src");

# source run variables
source(file.path(notebook_dir_sys,"Variables.R"));

# set for BU SCC
if(BUSCC==TRUE) {
  .libPaths('/project/cnsbomic/Tools/R');
  Sys.setenv(RSTUDIO_PANDOC='/usr/local/apps/rstudio-0.98.1103/rstudio-0.98.1103/bin/pandoc'); 
} 

# load rmarkdown and 
library(rmarkdown);

# run notebook
rmarkdown::render(file.path(notebook_dir,'Notebook.Rmd'),
                  knit_root_dir=notebook_dir,
                  output_file=paste("Analysis_Notebook_",Sys.Date(),".html", sep=""),
                  output_dir=working_dir);

