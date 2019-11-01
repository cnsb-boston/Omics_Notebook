#-------------------------------------------------
#' Draw Range vs Rank
#'
#' This function generates venn diagrams
#' 
#' @param matrix
#' @param outputpath output file path for plots
#' 
#' @examples
#' 
#' @export
drawRange <- function(matrix, outputpath=output_plots_path, file_name=""){
  
  require(scales)
  
  matrix <- matrix[rowMeans(matrix)!=0,]
  
  output_filename<-file.path(outputpath,paste("RangePlot","_",file_name,".pdf", sep=""))
  pdf(output_filename, width=3, height=3 )
  
  plot_data <- data.frame(Abundance = rowMeans(matrix),
                          Rank = rank(-rowMeans(matrix)) )
                                
  plot <- ggplot(data=plot_data, aes(x=Rank, y=Abundance)) + 
        geom_point(size=0.8) + 
        labs(title=paste("Range ",file_name, sep="")) +
        scale_y_continuous(trans=log10_trans(),
                           breaks=trans_breaks("log10", function(x) 10^x),
                           labels=trans_format("log10", math_format(10^.x)) ) +
        theme_bw() + theme(plot.margin=margin(10,20,10,10), axis.text.x=element_text(angle=45, hjust=1))
  print(plot+theme(legend.position="none"))
  
  dev.off()
}

