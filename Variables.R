# Where are the data files
working.dir <- "/projectnb/cnsbomic/ben/MQ_pipeline"

# have you already run this analysis and just want a new contrast?
newcontrastonly <- FALSE

# did you provide the txt folder from Maxquant for raw data QC?
txtFolder <- FALSE 

## Data types
isthereProteo <- TRUE
istherePhospho <- TRUE

## The following are time consuming and should be run only if needed.
enrichr_section <- TRUE #include enrichR
# GSEA
gsea_section <- TRUE # Include GSEA?
run_gsea <- TRUE # if GSEA has already been run, set to false and include results in directory
enrichment_map <- TRUE

## Files
maxq.global.filename <- "proteinGroups.txt"
maxq.phospho.filename <- "Phospho (STY)Sites.txt"
annotation.filename <- "Annotation.xlsx"

# Options
map.color <- "RdYlBu" # "viridis", "RdBu", "RdYlBu" or other ColorBrewer palette
zero.percent <- 0.25 # features must be detected in atleast this fraction of samples
adjpcutoff <- 0.05 # FDR cutoff for differential analysis
norm.method <- "loess" # 'quantile' or 'loess'

# Databases for enrichR search
search_databases <- c("KEGG_2016", "GO_Biological_Process_2017b", "GO_Cellular_Component_2017b",
                      "GO_Molecular_Function_2017b", "HMDB_Metabolites", 
                      "MSigDB_Oncogenic_Signatures");
