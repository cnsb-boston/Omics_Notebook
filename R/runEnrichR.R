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
  suppressMessages({
  if(run_seperate){
    enriched<- data.frame();
    genes_up <- genes[genes[,"logFC"]>0,"Gene"]
    genes_down <- genes[genes[,"logFC"]<0,"Gene"]
    
    if(length(genes_up)>0){
      invisible(capture.output({ enriched_up <- enrichr(unique(as.character(genes_up)), databases = search_dat); }))
      enriched_up <- bind_rows(enriched_up, .id="databases")
      enriched_up[,"Phenotype"] <- "+1"
      enriched <- rbind(enriched, enriched_up)
    }
    if(length(genes_down)>0){
      invisible(capture.output({ enriched_down <- enrichr(unique(as.character(genes_down)), databases = search_dat); }))
      enriched_down <- bind_rows(enriched_down, .id="databases")
      enriched_down[,"Phenotype"] <- "-1"
      enriched <- rbind(enriched, enriched_down)
    }
    enriched<- enriched[-1,]
    
  } else {
    invisible(capture.output({ enriched <- enrichr(unique(as.character(genes_down)), databases = search_dat); }))
    enriched <- bind_rows(enriched_down, .id="databases")
    enriched[,"Phenotype"] <- "+1"
  }
  })
  enriched <- enriched[order(enriched[,"Adjusted.P.value"]),]
  
  formatted <- enriched[,c("Term", "Term")];
  colnames(formatted) <- c("Term","Description")
  formatted[,c("p.Val","FDR")] <- format(enriched[,c("P.value", "Adjusted.P.value")], scientific=FALSE);
  formatted[,"Phenotype"] <- enriched[,"Phenotype"]
  enriched <- enriched %>% mutate(Genes=str_replace_all(Genes, ";", ","))
  formatted[,"Genes"] <- enriched[,"Genes"]
  formatted[,c("Gene_Hits", "Gene_Total")] <- str_split_fixed(enriched[,"Overlap"], "/",2)
  formatted[,"Gene_Ratio"] <- as.numeric(formatted[,"Gene_Hits"])/as.numeric(formatted[,"Gene_Total"])
  formatted[,c("Combined.Score")] <- format(enriched[,c("Combined.Score")], scientific=FALSE);
  formatted[,"Database"] <- enriched[,"databases"]
  #colnames(formatted) <- c("Term","Description","p.Val","FDR","Phenotype","Genes","Z.score","Combined.Score","Gene_Hits","Gene_Total")
  
  output_filename <- file.path(outputpath,paste("enrichr_", type,".txt", sep='')); 
  write.table(formatted, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
  
  drawEnrichment(enrichment_results=formatted[1:20,], type=type, label="enrichr", outputpath=outputpath);
}
