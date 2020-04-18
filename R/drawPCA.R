#-------------------------------------------------
#' Draw PCA
#'
#' This function generates PCA plots based on data
#' 
#' @param eset an ExpressionSet object with omcis data
#' @param type Type of data to be processed, or name for the Omics set.
#' @param outputpath output file path for plots
#'
#' @return a pca plot
#' 
#' @examples
#'     drawPCA(eset)
#' 
#' @export
drawPCA <- function(eset, x_axis="PC1", y_axis="PC2", type, outputpath=output_plots_path, 
                    outputfile=output_files_path, show_sample_names=TRUE, .species=species) {
  
  # perform PC analysis
  PC_data <- prcomp(t(exprs(eset)))
  percent_variance <- summary(PC_data)$importance["Proportion of Variance",]*100
  
  # Make colors for groups
  make_colors<-TRUE;
  if ("ColorsHex" %in% colnames(pData(eset)) ){ if(checkColor(pData(eset)[,"ColorsHex"])){
    colors_dots <- c(unique(as.character(pData(eset)[,"ColorsHex"])))
    names(colors_dots) <- unique(pData(eset)[,"Group"])
    make_colors<-FALSE;
  } }
  
  # Make PCA plot
  
  
  pca_graph <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis]))+ 
    { if( "Group2" %in% colnames(pData(eset)) ) { geom_point(aes(shape=pData(eset)$Group2,colour=pData(eset)$Group)); }} +
    { if( "Group2" %in% colnames(pData(eset)) ) { scale_shape_manual(values=1:length(unique(pData(eset)$Group2))) } else { geom_point(aes(colour=pData(eset)$Group)); }} + 
    theme_bw() + theme(legend.title=element_blank()) + 
    theme_bw() + theme(legend.title=element_blank()) + 
    { if (!make_colors) scale_color_manual(values=colors_dots);  } +
    labs(title=paste("PCA: \n ",type,sep=""),
       x=paste(x_axis, sprintf(" (%2.0f%%)", percent_variance[x_axis]), sep=""),
       y=paste(y_axis, sprintf(" (%2.0f%%)", percent_variance[y_axis]), sep=""));
  
  gg_dist_1 <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,x_axis], fill=pData(eset)$Group)) + 
    geom_density(alpha=0.4, size=0.2) + ylab(paste(x_axis, "Density", sep="\n") ) + theme_bw() +
    theme(legend.position="none", axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank(), 
          axis.title=element_text(size=8) )
  gg_dist_2 <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,y_axis], fill=pData(eset)$Group)) + 
    geom_density(alpha=0.4, size=0.2) + ylab(paste(y_axis, "Density", sep="\n") ) + theme_bw() +
    theme(legend.position="none", axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank(), 
          axis.title=element_text(size=8) )
    
  try({
  pca_graph2 <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,y_axis], y=PC_data$x[,"PC3"]))+ 
    { if( "Group2" %in% colnames(pData(eset)) ) { geom_point(aes(shape=pData(eset)$Group2,colour=pData(eset)$Group)); }} +
    { if( "Group2" %in% colnames(pData(eset)) ) { scale_shape_manual(values=1:length(unique(pData(eset)$Group2))) }  else { geom_point(aes(colour=pData(eset)$Group)); }} + 
    theme_bw() + theme(legend.title=element_blank()) + 
    { if (!make_colors) scale_color_manual(values=colors_dots);  } +
    labs(title=paste("PCA: \n ",type,sep=""),
       x=paste(y_axis, sprintf(" (%2.0f%%)", percent_variance[y_axis]), sep=""),
       y=paste("PC3", sprintf(" (%2.0f%%)", percent_variance["PC3"]), sep=""));
  
  gg_dist_3 <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,"PC3"], fill=pData(eset)$Group)) + 
    geom_density(alpha=0.4, size=0.2) + ylab(paste("PC3", "Density", sep="\n") ) + theme_bw() +
    theme(legend.position="none", axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank(), 
          axis.title=element_text(size=8) )
  })
  
  # Identify most variable hits for loadings graph
  sig_hits <- c( rownames(PC_data$rotation[order(PC_data$rotation[,x_axis], decreasing=T)[1:12],]),
                 rownames(PC_data$rotation[order(PC_data$rotation[,x_axis], decreasing=F)[1:12],]),
                 rownames(PC_data$rotation[order(PC_data$rotation[,y_axis], decreasing=T)[1:12],]),
                 rownames(PC_data$rotation[order(PC_data$rotation[,y_axis], decreasing=F)[1:12],]) )
  sig_data <- data.frame(PC_data$rotation[sig_hits,])
  
  pca_loading_graph <- ggplot(data=as.data.frame(PC_data$rotation), aes(x=(PC_data$rotation[,x_axis]), y=PC_data$rotation[,y_axis])) +
    geom_point(colour="grey", size=0.7) +
    geom_text_repel(data=sig_data, aes(x=sig_data[,x_axis], y=sig_data[,y_axis]),label=rownames(sig_data), colour="black", size=2) +
    theme_bw() + theme(legend.title=element_blank()) +
    labs(x=paste(x_axis, sep=""), y=paste(y_axis, sep=""),
         title=paste("PC Factor Loadings \n",type, sep=""))

  # Graph output
  output_filename <- file.path(outputpath, paste("PCAplots_",type,".pdf",sep=""))
  pdf(output_filename, width=4, height=4)
  print(pca_graph+theme(legend.position="none") )
  print(pca_graph+theme(legend.position="none")+geom_text_repel(aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis], label=rownames(PC_data$x)), colour="black", size=3) )
  grid.arrange(g_legend(pca_graph) )
  try({
    print(pca_graph2+theme(legend.position="none") )
    print(pca_graph2+theme(legend.position="none")+geom_text_repel(aes(x=PC_data$x[,"PC2"], y=PC_data$x[,"PC3"], label=rownames(PC_data$x)), colour="black", size=3) )
  })
  print(pca_loading_graph)
  plot(PC_data, type = 'l', main=paste("PCs vs. Variance",sep=""))
  
  piece1 <- gg_dist_1 + theme(axis.title.x=element_blank(), plot.margin=unit(c(0.5, -0.3, 0, 0.6), "cm") )
  piece2 <- gg_dist_2 + theme(axis.title.y=element_blank(), plot.margin=unit(c(-0.3, 0.2, 0.35, 0), "cm") ) + coord_flip() 
  piece3 <- pca_graph+theme(legend.position="none",plot.title=element_blank(), plot.margin=unit(c(0,0,0.5,0.5), "cm") )
  
  print(plot_grid( plot_grid( piece1, piece3, ncol=1, rel_heights=c(1,4)), 
             plot_grid(NULL, piece2, ncol=1, rel_heights=c(1,4)), 
             ncol=2, rel_widths=c(4,1)) )
  
  piece1 <- gg_dist_2 + theme(axis.title.x=element_blank(), plot.margin=unit(c(0.5, -0.3, 0, 0.6), "cm") )
  piece2 <- gg_dist_3 + theme(axis.title.y=element_blank(), plot.margin=unit(c(-0.3, 0.2, 0.35, 0), "cm") ) + coord_flip()
  piece3 <- pca_graph2+theme(legend.position="none",plot.title=element_blank(), plot.margin=unit(c(0,0,0.5,0.5), "cm") )
  
  print(plot_grid( plot_grid( piece1, piece3, ncol=1, rel_heights=c(1,4)), 
             plot_grid(NULL, piece2, ncol=1, rel_heights=c(1,4)), 
             ncol=2, rel_widths=c(4,1)) )
  
  
  dev.off()
  
  # Save analysis
  write.table(PC_data$rotation, file=file.path(outputfile, paste("PCA_loadings_", type,".txt", sep="")), quote=F, sep="\t" )
  
  # Run enrichment based on top 3 PCs
  if( .species!="Other" & ("Gene" %in% colnames(fData(eset))) ) { try({
  
    gsea_working_dir <- "PCA_Loading_GSEA"
    gsea_working_path <- file.path(outputfile, gsea_working_dir)
    if( dir.exists(gsea_working_path) == FALSE ) { dir.create(gsea_working_path) }
  
    analysis_names <- c(paste(type, "PC1", sep="_"), paste(type, "PC2", sep="_"), paste(type, "PC3", sep="_") )

    if(grepl("Mouse", .species) | grepl("Human", .species)){  suppressWarnings({ suppressMessages({
      # Only if you need a new GMT file
      if(grepl("Mouse", .species)){gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/"}
      if(grepl("Human", .species)){gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"}
      filenames = getURL(gmt_url)   #list all the files on the server
      tc = textConnection(filenames)
      contents = readLines(tc)
      close(tc)
      #get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA), start with gmt file that has pathways only
      rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",contents, perl = TRUE)
      gmt_file = unlist(regmatches(contents, rx))
      dest_gmt_file <- file.path(gsea_working_path,gmt_file)
      if(!file.exists(dest_gmt_file) ){ download.file(paste(gmt_url,gmt_file,sep=""),destfile=dest_gmt_file) }
    }) }) } else { try({
      geneset_lookup <-  read.delim(file.path(gsub("src", "data", notebook_dir),"geneset_table.txt"))
      gmt_file <- as.character(geneset_lookup[geneset_lookup[,"Species"]==species,"GeneSet"])
      dest_gmt_file <-file.path(gsea_working_path, gmt_file)
      if(!file.exists(dest_gmt_file) ){ file.copy(from=file.path( gsub("src", "data", notebook_dir), "species_genesets",gmt_file),
                                                  to= dest_gmt_file )  }
    }) }
  
    for (i in 1:length(analysis_names)){
      analysis_name <- analysis_names[[i]]
      ranked_vector <- PC_data$rotation[,i]
      names(ranked_vector) <- fData(eset)[,"Gene"]
      pathways_gmt <- gmtPathways(dest_gmt_file)
      suppressWarnings({
        fgsea_results <- fgsea(pathways=pathways_gmt, 
                             stats=ranked_vector,
                             minSize=15, 
                             maxSize=500, 
                             nperm=10000)
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
      output_filename <- file.path(gsea_working_path, paste("fgsea_", analysis_name, "_loadings.txt", sep=""))
      write.table(fgsea_out, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
    }
  
  }, silent=T) } 
  
  return(pca_graph)
}

#-------------------------------------------------
#' Draw UMAP
#'
#' This function generates UMAP plots
#' 
#' @param eset an ExpressionSet object with omcis data
#' @param type Type of data to be processed, or name for the Omics set.
#' @param outputpath output file path for plots
#'
#' @return a plot
#' 
#' @examples
#'     drawUMAP(eset)
#' 
#' @export
drawUMAP <- function(eset, type, outputpath=output_plots_path) {
  
  data <- t(exprs(eset))
  data.labels <- pData(eset)$Group
  
  use_size <- 15
  if(length(data.labels)<15){
    use_size <- ( length(data.labels)/2 )
  }
  data.umap <- uwot::umap(data, n_neighbors= use_size ) 

  
  # Make colors for groups
  make_colors<-TRUE;
  if ("ColorsHex" %in% colnames(pData(eset)) ){ if(checkColor(pData(eset)[,"ColorsHex"])){
    colors_dots <- c(unique(as.character(pData(eset)[,"ColorsHex"])))
    names(colors_dots) <- unique(pData(eset)[,"Group"])
    make_colors<-FALSE;
  } }
  
  # Make plot
  umap_graph <- ggplot(data=as.data.frame(data.umap), aes(x=data.umap[,1], y=data.umap[,2]))+ 
    { if( "Group2" %in% colnames(pData(eset)) ) { geom_point(aes(shape=pData(eset)$Group2,colour=pData(eset)$Group)); } } +
    { if( "Group2" %in% colnames(pData(eset)) ) { scale_shape_manual(values=1:length(unique(pData(eset)$Group2))) } else { geom_point(aes(colour=pData(eset)$Group)); }} + 
    theme_bw() + theme(legend.title=element_blank()) + 
    { if (!make_colors) scale_color_manual(values=colors_dots);  } +
    labs(title=paste("UMAP: \n ",type,sep="") );
  
  # Graph output
  output_filename <- file.path(outputpath, paste("UMAPplot_",type,".pdf",sep=""))
  pdf(output_filename, width=4, height=4)
  
  print(umap_graph+theme(legend.position = "none"))
  grid.arrange(g_legend(umap_graph) )
  
  dev.off()
  
}


