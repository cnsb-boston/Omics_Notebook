# Get Notebook directory
notebook_dir_sys <- tail(commandArgs(trailingOnly=TRUE), n=1);

# source run variables
source(paste(notebook_dir_sys,"Variables.R", sep="/"));

# set CNSB library path
.libPaths('/project/cnsbomic/Tools/R');

# load rmarkdown and set rstudio version
library(rmarkdown);
Sys.setenv(RSTUDIO_PANDOC='/usr/local/apps/rstudio-0.98.1103/rstudio-0.98.1103/bin/pandoc');

# run notebook
rmarkdown::render(paste(notebook_dir,'Proteomics_Notebook.Rmd', sep="/"),
                  knit_root_dir=notebook_dir,
                  output_file=paste("Proteomics_Notebook_",Sys.Date(),".html", sep=""),
                  output_dir=working.dir);

