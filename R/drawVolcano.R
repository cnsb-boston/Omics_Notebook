#-------------------------------------------------
#' Draw Volcano Plots
#'
#' This function makes a volcano plot
#' 
#' @param dat is the result of topTable()
#' @param type Type of data to be processed, or a name for the Omics analysis
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#'  
#' 
#' @export
drawVolcano <- function(dat, type, outputpath=output_contrast_path){ 
  # Gene labels
  label_names <- c(rownames(dat[dat$logFC>0,][order(dat$P.Value, decreasing=FALSE),])[1:6],
                   rownames(dat[dat$logFC<0,][order(dat$P.Value, decreasing=FALSE),])[1:6],
                   rownames(dat[order(dat$logFC, decreasing=FALSE),])[1:4],
                   rownames(dat[order(dat$logFC, decreasing=TRUE),])[1:4] )
  
  label_names <- unique(label_names)
  
  # Dot coloring
  dat$Significance <- "NS"
  dat$Significance[(dat$P.Value< 0.05 & sign(dat$logFC)>0) ] <- "Up"
  dat$Significance[(dat$P.Value< 0.05 & sign(dat$logFC)<0) ] <- "Down"
  dat$Signicance <- factor(dat$Significance, levels="NS", "Up", "Down")
  
  # Plot
  plot1 <- ggplot(data=data.frame(dat)) + geom_point(aes(x=logFC, y=-log10(P.Value),colour=(dat$Significance))) +
    scale_colour_manual(values=c(NS="grey", Up="Red", Down="blue"))   +
    labs(title=paste(type,"Volcano Plot", sep=" ")) + theme_bw() + 
    theme(legend.title=element_blank()) +
    { if( "Gene" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]), aes(x=logFC, y=-log10(P.Value), label=Gene))
      else if( "mz" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]), aes(x=logFC, y=-log10(P.Value), label=round(mz, digits=2)))
      else geom_text_repel(data=data.frame(dat[label_names,]), aes(x=logFC, y=-log10(P.Value), label=feature_identifier))
    }
  
  output_filename <- file.path(outputpath, paste("volcano_",type,".pdf", sep=''));
  pdf(output_filename);
  plot(plot1)
  tmp<-dev.off();
}
drawMDPlot <- function(dat, type, outputpath=output_contrast_path){ 
  # Gene labels
  label_names <- c(rownames(dat[order(dat$logFC, decreasing=FALSE),])[1:10],
                   rownames(dat[order(dat$logFC, decreasing=TRUE),])[1:10] )
  
  # Dot coloring
  dat$Significance <- "NS"
  dat$Significance[dat$logFC>0.5 ] <- "Up"
  dat$Significance[dat$logFC<-0.5 ] <- "Down"
  dat$Signicance <- factor(dat$Significance, levels="NS", "Up", "Down")
  
  # Plot
  plot1 <- ggplot(data=data.frame(dat)) + geom_point(aes(x=Mean, y=logFC,colour=(dat$Significance))) +
    scale_colour_manual(values=c(NS="grey", Up="Red", Down="blue"))   +
    labs(title=paste(type,"MD Plot", sep=" ")) + theme_bw() + 
    theme(legend.title=element_blank()) +
    { if( "Gene" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]), aes(x=Mean, y=logFC, label=Gene))
      else if( "mz" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]), aes(x=Mean, y=logFC, label=round(mz, digits=2)))
      else geom_text_repel(data=data.frame(dat[label_names,]), aes(x=Mean, y=logFC, label=feature_identifier))
    }

output_filename <- file.path(outputpath, paste("MDplot_",type,".pdf", sep=''));
pdf(output_filename);
plot(plot1)
tmp<-dev.off();
}

