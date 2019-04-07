#-------------------------------------------------
#' Function to draw enrichment results
#'
#' This function draws plots for enrichment results
#' 
#' @param enrichment_results
#' @param type Type of data to be processed, or name for the Omics set
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#' 
#' @export
drawEnrichment <- function(enrichment_results, type, label, outputpath=output_contrast_path ) {
  
  enrichment_results <- data.frame(enrichment_results)
  enrichment_results[,"Gene_Hits"] <- as.numeric(enrichment_results[,"Gene_Hits"])
  enrichment_results[,"Gene_Ratio"] <- as.numeric(enrichment_results[,"Gene_Ratio"])
  enrichment_results[,"FDR"] <- as.numeric(enrichment_results[,"FDR"])
  enrichment_results[,"Phenotype"] <- as.factor(enrichment_results[,"Phenotype"])
  if ("NES" %in% colnames(enrichment_results)){
    enrichment_results[,"Enrichment_Score"] <- as.numeric(enrichment_results[,"NES"])
  } else {
    enrichment_results[,"Enrichment_Score"] <- -log10(enrichment_results[,"FDR"])*sign(as.numeric(enrichment_results[,"Phenotype"]))
  }
  enrichment_results[,"Term"] <- as.character( enrichment_results[,"Term"])
  enrichment_results[,"Term"] <- factor(enrichment_results[,"Term"], levels=enrichment_results[order(enrichment_results[,"Gene_Ratio"]),"Term"] )
  
  plot <- ggplot(enrichment_results, aes(x=Gene_Ratio, y=Term) ) + 
    geom_point(aes(size=Gene_Hits, color=Enrichment_Score)) + 
    scale_color_continuous(low="blue", high="red") +
    scale_x_continuous(limits=c(0, (max(enrichment_results[,"Gene_Ratio"])*1.2))) +
    theme_linedraw() + 
    labs(title=paste("Enrichment Results: ",type, sep=''),
         x="Gene Ratio", y="Enrichment Terms");
  
  output_filename <- file.path(outputpath,paste(label,"_results_", type,".pdf", sep='')); 
  pdf(output_filename, width=10, height=6)
  print(plot)
  dev.off()
}