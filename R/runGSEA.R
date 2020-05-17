#-------------------------------------------------
#' Wrapper fnction for GSEA
#'
#' This runs enrichment analysis with fgsea
#' 
#' @param genes, list of genes
#' @param type Type of data to be processed, or name for the Omics set
#' @param search_data list of enrichr databases to search
#' @param outputpath output file path for plots
#' @param run_seperate Whether or not to run up and down genes seperately
#'
#' @return 
#' 
#' @examples
#' 
#' @export
runGSEA <- function (rnk, gmt,
                     analysisName="GSEA", pathway_plots= F,
                     out_path=getwd()   ) {
  
  try({
    
    gsea_images_path <- file.path(out_path, "Images")
    gsea_absval_path <- file.path(out_path, "AbsoluteValueRank")
    if( dir.exists(gsea_images_path) == FALSE ) { dir.create(gsea_images_path) }
    if( dir.exists(gsea_absval_path) == FALSE ) { dir.create(gsea_absval_path) }
    
  ranked_features <- read.delim(rnk)
  file.copy(from=rnk, to=file.path(out_path,basename(rnk)) )
  ranked_vector <- ranked_features[,"rank"]
  names(ranked_vector) <- toupper(ranked_features[,"GeneName"])
  pathways_gmt <- gmtPathways(gmt)
  pathways_gmt <- sapply(pathways_gmt, toupper)
  
  suppressWarnings({
    fgsea_results <- fgsea(pathways=pathways_gmt, 
                           stats=ranked_vector,
                           minSize=15, 
                           maxSize=500, 
                           nperm=1000)
  })
  
  fgsea_out <- fgsea_results[,c(1,1,2,3)]
  colnames(fgsea_out) <- c("Term", "Description", "p.Val", "FDR")
  fgsea_out[,"Phenotype"] <- sign(fgsea_results[,"NES"])
  fgsea_out[,"Genes"] <- apply(fgsea_results[,"leadingEdge"], 1, function(x) paste(as.character(unlist(x["leadingEdge"])), collapse=","))
  fgsea_out[,"Genes"] <- apply(fgsea_out, 1, function(x) gsub(" ", "", x["Genes"]))
  fgsea_out[,"NES"] <- (fgsea_results[,"NES"])
  fgsea_out[,"ES"] <- (fgsea_results[,"ES"])
  fgsea_out[,"Gene_Hits"] <- apply(fgsea_results, 1, function(x) length(unlist(x["leadingEdge"])) )
  fgsea_out[,"Gene_Total"] <- fgsea_results[,"size"]
  fgsea_out <- fgsea_out[order(fgsea_out[,"p.Val"], decreasing=F)]
  output_filename <- file.path(out_path, paste("fgsea_", analysisName, ".txt", sep=""))
  write.table(fgsea_out, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
  
  output_filename <- file.path(gsea_images_path, paste("fgsea_", analysisName, ".pdf", sep=""))
  topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(output_filename, width=20, height=8)
  plotGseaTable(pathways_gmt[topPathways], ranked_vector, fgsea_results, gseaParam = 0.5)
  dev.off()
  
  if(pathway_plots){
  topPathways <- fgsea_results[order(pval, decreasing=F), pathway]
  if( length(topPathways) > 50 ) { topPathways <- topPathways[1:50]}

  output_filename <- file.path(gsea_images_path, paste("fgsea_", analysisName, "_allPathways.pdf", sep=""))
  pdf(output_filename, width=4, height=3)
  for(pathway_index in 1:length(topPathways)){ try({
    print(plotEnrichment(pathways_gmt[[topPathways[pathway_index] ]], ranked_vector) +
            labs(title=topPathways[pathway_index]) + theme(plot.title=element_text(size=6)) )
  }, silent=TRUE) }
  dev.off()
  }
  
  ranked_vector <- abs(ranked_features[,"rank"])
  names(ranked_vector) <- toupper(ranked_features[,"GeneName"])
  
  suppressWarnings({
    fgsea_results <- fgsea(pathways=pathways_gmt, 
                           stats=ranked_vector,
                           minSize=15, 
                           maxSize=500, 
                           nperm=1000)
  })
  
  fgsea_out <- fgsea_results[,c(1,1,2,3)]
  colnames(fgsea_out) <- c("Term", "Description", "p.Val", "FDR")
  fgsea_out[,"Phenotype"] <- sign(fgsea_results[,"NES"])
  fgsea_out[,"Genes"] <- apply(fgsea_results[,"leadingEdge"], 1, function(x) paste(as.character(unlist(x["leadingEdge"])), collapse=","))
  fgsea_out[,"Genes"] <- apply(fgsea_out, 1, function(x) gsub(" ", "", x["Genes"]))
  fgsea_out[,"NES"] <- (fgsea_results[,"NES"])
  fgsea_out[,"ES"] <- (fgsea_results[,"ES"])
  fgsea_out[,"Gene_Hits"] <- apply(fgsea_results, 1, function(x) length(unlist(x["leadingEdge"])) )
  fgsea_out[,"Gene_Total"] <- fgsea_results[,"size"]
  fgsea_out <- fgsea_out[order(fgsea_out[,"p.Val"], decreasing=F)]
  output_filename <- file.path(gsea_absval_path, paste("fgsea_", analysisName, "_AbsVal.txt", sep=""))
  write.table(fgsea_out, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
  
}) }
