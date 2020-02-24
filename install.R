
## Unbundle packrat libraries
install.packages("packrat")

packrat::unbundle("Omics_Notebook-2019-11-04.tar.gz", where="./")

# 
# 
# 
# #####################################
# ##
# ##  Libraries installation
# ##
# ####################################
# 
# install.packages("rlang")
# install.packages("Rcpp", dependencies=TRUE)
# install.packages("R.utils")
# install.packages("knitr")
# install.packages("ggplot2", dependencies=TRUE)
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("ggrepel")
# install.packages("calibrate")
# install.packages("RColorBrewer")
# install.packages("viridis")
# install.packages("reshape2")
# install.packages("gridExtra")
# install.packages("lattice")
# install.packages("VennDiagram")
# install.packages("heatmaply")
# install.packages("openxlsx")
# install.packages("UpSetR")
# install.packages("corrplot")
# install.packages("uwot")
# install.packages("cowplot")
# 
# install.packages("rmarkdown")
# install.packages("knitrBootstrap")
# 
# install.packages("BiocManager")
# 
# BiocManager::install("Biobase")
# BiocManager::install("limma")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("preprocessCore")
# BiocManager::install("edgeR")
# BiocManager::install("Glimma")
# 
# 
# BiocManager::install("fgsea")
# 
# 
# # FOR COMBAT BATCH CORRECTION
# BiocManager::install("sva")
# #######
# 
# # FOR TXT (Maxquant) RAW DATA REPORT
# install.packages("PTXQC")
# 
# # # FOR METABOLITES ONLY
# metanr_packages <- function(){
#   metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", "randomForest", "caTools", "e1071", "som", "impute", "pcaMethods", "RJSONIO", "ROCR",
#                  "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "pheatmap", "SSPA", "sva", "Rcpp", "pROC", "data.table",
#                  "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", "caret", "lattice", "igraph", "gplots",
#                  "KEGGgraph", "reshape", "RColorBrewer", "tibble", "siggenes", "plotly", "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel", "metap", "reshape2", "scales")
#   list_installed <- installed.packages()
#   new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
#   if(length(new_pkgs)!=0){
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install(new_pkgs)#, version = "3.9")
#     print(c(new_pkgs, " packages added..."))
#   }
#   if((length(new_pkgs)<1)){ print("No new packages added...") }
# }
# metanr_packages()
# BiocManager::install("ncdf4")
# 
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cnsb-boston/MetaboAnalystR", build_vignettes=FALSE)
# 
# 
# #######
# 
# ## ENRICHMENT
# 
# # For Enrichr
# install.packages("enrichR")
# 
# # Kinase Enrichment KSEA
# #install.packages("KSEAapp") created fork to add output directory option
# library(devtools)
# devtools::install_github("cnsb-boston/KSEAapp")
# 
# # FOR GSEA SECTION
# install.packages("RCurl")
# install.packages("httr")
# install.packages("RJSONIO")
# 
# # FOR ENRICHMENT MAP
# BiocManager::install("RCy3")
# #######
# 
# 


