#Source this file from an R session to debug the most recent analysis
library(OmicsNotebook)
source("Parameters.R")
dirs=list.files(".",paste0(project_name,"_[0-9]+$"))
g=readRDS(file.path(dirs[length(dirs)],"1_Files","Data_g_state.RDS"))
g=OmicsNotebook:::setup_analysis_dir(g, Sys.getenv("TMPDIR"))

# basic calls that don't rely on files
calls=c("g.notebook.setup","g.make.omicsList","g.make.contrasts","g.norm","g.norm.to.first","g.use.site.norm","g.param.norm","g.combine.met","g.limma","g.fc")
g$calls=calls[calls %in% g$calls]
