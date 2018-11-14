

#####################################
##
##  Libraries installation
##
####################################

install.packages("rlang")
install.packages("Rcpp")

install.packages("knitr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggrepel")
install.packages("calibrate")
install.packages("RColorBrewer")
install.packages("viridis")
install.packages("reshape2")
install.packages("gridExtra")
install.packages("lattice")
install.packages("enrichR")
install.packages("VennDiagram")
install.packages("heatmaply")
install.packages("openxlsx")
install.packages("UpSetR")
install.packages("corrplot")

install.packages("rmarkdown")

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biobase")
biocLite("limma")
biocLite("ComplexHeatmap")
biocLite("preprocessCore")
biocLite("edgeR")
biocLite("Glimma")

# FOR COMBAT BATCH CORRECTION
biocLite("sva")
#######

# FOR TXT RAW DATA REPORT
install.packages("PTXQC")

# FOR METABOLITES ONLY
metanr_packages <- function(){
  metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", "randomForest", "caTools", "e1071", "som", "impute", "pcaMethods", "RJSONIO", "ROCR", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "pheatmap", "SSPA", "sva", "Rcpp", "pROC", "data.table", "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", "caret", "lattice", "igraph", "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", "siggenes", "plotly")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){
    source("https://bioconductor.org/biocLite.R")
    biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
    print(c(new_pkgs, " packages added..."))
  }
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()
install.packages("devtools")
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR")#, build_vignettes=TRUE)


biocLite("hmdbQuery") 
#######

# FOR GSEA SECTION
install.packages("RCurl")
install.packages("httr")
install.packages("RJSONIO")
#######

# FOR ENRICHMENT MAP
biocLite("RCy3")
#######




