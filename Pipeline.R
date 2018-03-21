# Where are the data files
source("/project/cnsbomic/Omics_Notebook/Variables.R");
.libPaths('/project/cnsbomic/Tools/R');
library(rmarkdown);
Sys.setenv(RSTUDIO_PANDOC='/usr/local/apps/rstudio-0.98.1103/rstudio-0.98.1103/bin/pandoc');
rmarkdown::render('/project/cnsbomic/Omics_Notebook/MaxQuant_Notebook.Rmd',
                  output_file=paste("MaxQuant_Notebook_",Sys.Date(),".html", sep=""),
                  output_dir=working.dir);
