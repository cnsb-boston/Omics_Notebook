#-------------------------------------------------
#' Wrapper fnction for EnrichR
#'
#' This runs enrichment analysis with enrichr
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
runEnrichR <- function (genes, type, search_dat=search_databases, 
                        outputpath=enrichr_working_path, run_seperate=FALSE) {

  mergerows = function(enrich_list, id=NULL){
    dplyr::bind_rows(enrich_list[sapply(enrich_list,nrow)>0], .id=id)
  }

  mergedb = function(enrich_list, pheno){
    ret = mergerows(enrich_list, id="databases")
    ret$Phenotype = pheno
    ret
  }

  suppressMessages({
  if(run_seperate){
    elist <- list()
    genes_up <- genes[genes[,"logFC"]>0,"Gene"]
    genes_down <- genes[genes[,"logFC"]<0,"Gene"]
    
    if(length(genes_up)>0){
      invisible(capture.output({ enriched_up <- enrichR::enrichr(unique(as.character(genes_up)), databases = search_dat); }))
      elist <- c(elist, list(mergedb(enriched_up, "+1")))
    }
    if(length(genes_down)>0){
      invisible(capture.output({ enriched_down <- enrichR::enrichr(unique(as.character(genes_down)), databases = search_dat); }))
      elist <- c(elist, list(mergedb(enriched_down, "-1")))
    }
    enriched<- mergerows(elist)
    
  } else {
    invisible(capture.output({ enriched <- enrichR::enrichr(unique(as.character(genes_down)), databases = search_dat); }))
    enriched <- mergedb(enriched, "+1")
  }
  })
  enriched <- enriched[order(enriched[,"Adjusted.P.value"]),]
  
  formatted <- enriched[,c("Term", "Term")];
  colnames(formatted) <- c("Term","Description")
  formatted[,c("p.Val","FDR")] <- format(enriched[,c("P.value", "Adjusted.P.value")], scientific=FALSE);
  formatted[,"Phenotype"] <- enriched[,"Phenotype"]
  enriched <- dplyr::mutate(enriched, Genes=stringr::str_replace_all(Genes, ";", ","))
  formatted[,"Genes"] <- enriched[,"Genes"]
  formatted[,c("Gene_Hits", "Gene_Total")] <- stringr::str_split_fixed(enriched[,"Overlap"], "/",2)
  formatted[,"Gene_Ratio"] <- as.numeric(formatted[,"Gene_Hits"])/as.numeric(formatted[,"Gene_Total"])
  formatted[,c("Combined.Score")] <- format(enriched[,c("Combined.Score")], scientific=FALSE);
  formatted[,"Database"] <- enriched[,"databases"]
  #colnames(formatted) <- c("Term","Description","p.Val","FDR","Phenotype","Genes","Z.score","Combined.Score","Gene_Hits","Gene_Total")
  
  output_filename <- file.path(outputpath,paste("enrichr_", type,".txt", sep='')); 
  write.table(formatted, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
  
  drawEnrichment(enrichment_results=formatted[1:20,], type=type, label="enrichr", outputpath=outputpath);
}
