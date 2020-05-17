#-------------------------------------------------
#' Draw XY Correlation Plot
#'
#' This function generates correlation plots
#' 
#' @param item_list a list of tibbles
#' @param item_name a name for the items
#' @param outputpath output file path for plots
#' @param file_name suffix for plot filename
#' 
#' @examples
#' 
#' @export
drawXYCorr <- function(item_list, item_name, outputpath=output_plots_path, file_name="", subset_genes=FALSE){
  
  output_filename<-file.path(outputpath,paste("Correlation_Plots_",item_name,"_",file_name,".pdf", sep=""))
  pdf(output_filename, width=3, height=3 )
  for (j in 1:(length(item_list)-1)){
    for (k in (j+1):length(item_list)){
      # make data to plot
      plot_data <- merge(item_list[[j]], item_list[[k]], by.x=item_name, by.y=item_name)
      colnames(plot_data) <- c(item_name, names(item_list)[j], names(item_list)[k])
      # density colors
      x<-densCols(plot_data[,2], plot_data[,3],colramp=colorRampPalette(c("black","white")))
      plot_data$Density <- col2rgb(x)[1,] + 1L
      # make plot
      corr_coef <- cor(plot_data[,2], plot_data[,3], method="pearson");
      label<-paste("italic(r) == ",round(corr_coef, digits=2), sep="");
      plot <- ggplot(data=plot_data, aes(x=plot_data[,2], y=plot_data[,3], color=Density)) + 
        geom_point() + scale_color_viridis(direction=-1) + labs(x=names(item_list)[j], y=names(item_list)[k]) +
        labs(title=paste("Avg. ",item_name, " ", file_name,sep=""))+#,": \n",names(item_list)[j]," vs. ",names(item_list)[k],sep="" )) +
        theme_bw() + annotate("text",x=-Inf, y=Inf, label=label, color="black", vjust=1.5, hjust=-0.4 , parse=TRUE) +
        geom_smooth(method=lm, se=FALSE, color="red")
      print(plot+theme(legend.position="none"))
      grid.arrange(g_legend(plot))
    }
  }
  
  # Repeat for subset genes
  if( class(subset_genes)!="logical"){
    for( i in 1:length(subset_genes) ){ try({
    for (j in 1:(length(item_list)-1)){ # FILTER GENES BASED ON SUBSET
      for (k in (j+1):length(item_list)){
        # make data to plot
        plot_data <- merge(item_list[[j]], item_list[[k]], by.x=item_name, by.y=item_name)
        colnames(plot_data) <- c(item_name, names(item_list)[j], names(item_list)[k])
        plot_data <- plot_data[which(plot_data[,item_name] %in% subset_genes[[i]]),]
        # density colors
        x<-densCols(plot_data[,2], plot_data[,3],colramp=colorRampPalette(c("black","white")))
        plot_data$Density <- col2rgb(x)[1,] + 1L
        # make plot
        corr_coef <- cor(plot_data[,2], plot_data[,3], method="pearson");
        label<-paste("italic(r) == ",round(corr_coef, digits=2), sep="");
        plot <- ggplot(data=plot_data, aes(x=plot_data[,2], y=plot_data[,3], color=Density)) + 
          geom_point() + scale_color_viridis(direction=-1) + labs(x=names(item_list)[j], y=names(item_list)[k]) +
          labs(title=paste("Avg. ",item_name, " ", file_name," in ",names(subset_genes)[i],": \n",
                           names(item_list)[j]," vs. ",names(item_list)[k],sep="" )) +
          theme_bw() + annotate("text",x=-Inf, y=Inf, label=label, color="black", vjust=1.5, hjust=-0.4 , parse=TRUE) +
          geom_smooth(method=lm, se=FALSE, color="red")
        print(plot+theme(legend.position="none"))
        grid.arrange(g_legend(plot))
      }
    }
    }) }
  }
  dev.off()
}

