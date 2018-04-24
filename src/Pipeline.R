# Get Notebook directory
notebook_dir_sys <- tail(commandArgs(trailingOnly=TRUE), n=1);
notebook_dir_sys <- paste(notebook_dir_sys, "/src", sep="");

# source run variables
source(paste(notebook_dir_sys,"Variables.R", sep="/"));

# set CNSB library path
.libPaths('/project/cnsbomic/Tools/R');

# load rmarkdown and set rstudio version
library(rmarkdown);
Sys.setenv(RSTUDIO_PANDOC='/usr/local/apps/rstudio-0.98.1103/rstudio-0.98.1103/bin/pandoc');

# run notebook
rmarkdown::render(paste(notebook_dir,'Notebook.Rmd', sep="/"),
                  knit_root_dir=notebook_dir,
                  output_file=paste("Analysis_Notebook_",Sys.Date(),".html", sep=""),
                  output_dir=working_dir);

