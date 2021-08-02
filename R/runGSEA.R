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
                           nperm=10000)
  })
  
  fgsea_out <- fgsea_results[,c(1,1,2,3)]
  colnames(fgsea_out) <- c("Term", "Description", "p.Val", "FDR")
  fgsea_out[,"Phenotype"] <- sign(fgsea_results[,"ES"])
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
                           nperm=10000)
  })
  
  fgsea_out <- fgsea_results[,c(1,1,2,3)]
  colnames(fgsea_out) <- c("Term", "Description", "p.Val", "FDR")
  fgsea_out[,"Phenotype"] <- sign(fgsea_results[,"ES"])
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

#-------------------------------------------------
#' Fetch a new GMT file
#'
#' Fetch a new GMT file for a given species if it does not yet exist at the target destination
#' 
#' @param dest_path path to store the gmt file in
#' @param .species species to use for the download
#'
#' @return the gmt file path
#' 
#' @examples
#' 
#' @export
fetchGMT <- function(dest_path,.species="Other"){
  if(grepl("Mouse|Human", .species)){  suppressWarnings({ suppressMessages({
    # Only if you need a new GMT file
    gmt_url = sub("(Mouse|Human).*","http://download.baderlab.org/EM_Genesets/current_release/\\1/symbol/",.species)
    filenames = getURL(gmt_url)   #list all the files on the server
    tc = textConnection(filenames)
    contents = readLines(tc)
    close(tc)
    #get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA), start with gmt file that has pathways only
    rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",contents, perl = TRUE)
    gmt_file = unlist(regmatches(contents, rx))
    dest_gmt_file <- file.path(dest_path,gmt_file)
    if(!file.exists(dest_gmt_file)) { download.file(paste0(gmt_url,gmt_file),destfile=dest_gmt_file) }
  }) }) } else { try({
    geneset_lookup <-  read.delim(file.path(gsub("src", "data", notebook_dir),"geneset_table.txt"))
    gmt_file <- as.character(geneset_lookup[geneset_lookup[,"Species"]==.species,"GeneSet"])
    dest_gmt_file <-file.path(dest_path, gmt_file)
    if(!file.exists(dest_gmt_file)) { file.copy(from=file.path( gsub("src", "data", notebook_dir), "species_genesets",gmt_file),
                                                to= dest_gmt_file ) }
  }) }
  dest_gmt_file
}
