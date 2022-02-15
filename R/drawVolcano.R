getLabelVars <- function(dat){
  if( "Gene" %in% colnames(dat)){
    textsize <- c(2,2)
    labname <- rep("Gene",2)
  } else if( "mz" %in% colnames(dat)) {
    textsize <- c(2,1.5)
    labname <- c("round(mz, digits=2)","feature_identifier")
  } else {
    textsize <- c(1,1.5)
    labname <- rep("feature_identifier",2)
  }
  list(textsize=textsize,labname=labname)
}

drawVPlots <- function(dat, xvar, yvar, title, label_names, colorby, subset_rows){
  lv <- getLabelVars(dat)
  suppressWarnings({
    plot1 <- ggplot(data=data.frame(dat)) + 
      geom_point(aes_string(x=xvar, y=yvar, colour=colorby),size=0.7, pch=19) +
      scale_colour_manual(values=c(NS="Grey", Up="Red", Down="Blue"))   +
      labs(title=title) + theme_bw() + 
      theme(legend.title=element_blank())
    plot2 <- plot1 + ggrepel::geom_text_repel(data=data.frame(dat[label_names,]),size=lv$textsize[1],
                                     aes_string(x=xvar, y=yvar, label=lv$labname[1])) +
                     geom_point(data=data.frame(dat[label_names,]) ,size=0.7, pch=21,
                                aes_string(x=xvar, y=yvar))
   
    plot(plot1+theme(legend.position="none"))
    plot(plot2+theme(legend.position="none"))
    if(class(subset_rows)!="logical" ){
      for( j in 1:length(subset_rows)){ try({
        if(length(subset_rows[[j]])>100){subset_rows[[j]]<- subset_rows[[j]][1:100]}
        if(length(subset_rows[[j]])>0){
          label_names2 <- subset_rows[[j]]
          plot(plot1 +
                 ggrepel::geom_label_repel(data=data.frame(dat[label_names,]),size=lv$textsize[2], 
                                  aes_string(x=xvar, y=yvar, label=lv$labname[2]),
                                  label.size=NA, label.padding = .1, na.rm=T, fill=alpha(c("white"),0.5) )+
                 theme(legend.position="none") +
                 geom_point(data=data.frame(dat[label_names2,]), size=0.7, pch=21,
                            aes_string(x=xvar, y=yvar)) )
        }
      }) }
    }
  })
  plot1
}

#-------------------------------------------------
#' Draw Volcano Plots
#'
#' This function makes a volcano plot
#' 
#' @param dat is the result of topTable()
#' @param type Type of data to be processed, or a name for the Omics analysis
#' @param subset_rows
#' @param outputpath output file path for plots
#' 
#'
#'
#' 
#' @examples
#'  
#' 
#' @import ggplot2
#' @export
drawVolcano <- function(dat, type, subset_rows=F,
                        top_values=0.05, top_fc=0,
                        outputpath=output_contrast_path){ 
  # Gene labels
  label_names <- c(rownames(dat[dat$logFC>0,][order(dat$P.Value, decreasing=FALSE),])[1:6],
                   rownames(dat[dat$logFC<0,][order(dat$P.Value, decreasing=FALSE),])[1:6],
                   rownames(dat[order(dat$logFC, decreasing=FALSE),])[1:4],
                   rownames(dat[order(dat$logFC, decreasing=TRUE),])[1:4] )
  
  label_names <- unique(label_names)
  
  # Dot coloring
  for(v in c("P.Value","adj.P.Val")){
    s <- paste0(v,".sig")
    dat[,s] <- "NS"
    dat[,s] [( dat[,v]<=top_values & dat$logFC>=top_fc ) ] <- "Up"
    dat[,s] [( dat[,v]<=top_values & dat$logFC<=(-1*top_fc) ) ] <- "Down"
    dat[,s] <- factor(dat[,s], levels=c("NS", "Up", "Down"))
  }

  # Plot
  output_filename <- file.path(outputpath, paste0(type,"_volcano",".pdf"));
  pdf(output_filename, width=3.5, height=3.5);
  for(v in c("P.Value","adj.P.Val")){
    plot1 <- drawVPlots(dat, xvar="logFC", yvar=paste0("-log10(",v,")"), title=paste(type,"Volcano Plot",sep="\n"), label_names=label_names, colorby=paste0(v,".sig"), subset_rows=subset_rows)
  }

  gridExtra::grid.arrange(g_legend(plot1))
  tmp<-dev.off();
}

drawMDPlot <- function(dat, type, subset_rows=FALSE, outputpath=output_contrast_path, cutoff=1){ 
  # Gene labels
  label_names <- c(rownames(dat[order(dat$logFC, decreasing=FALSE),])[1:10],
                   rownames(dat[order(dat$logFC, decreasing=TRUE),])[1:10] )
  
  # Dot coloring
  dat$Significance <- "NS"
  dat$Significance[dat$logFC>cutoff ] <- "Up"
  dat$Significance[(dat$logFC<(-1*cutoff)) ] <- "Down"
  dat$Significance <- factor(dat$Significance, levels=c("NS", "Up", "Down"))
  
  # Plot
  output_filename <- file.path(outputpath, paste0(type,"_MDplot",".pdf"));
  pdf(output_filename, width=3.5, height=3.5);
  plot1 <- drawVPlots(dat, xvar="Mean", yvar="logFC", title=paste(type,"MD Plot", sep="\n"), label_names=label_names, colorby="Significance", subset_rows=subset_rows)
  gridExtra::grid.arrange(g_legend(plot1))
  tmp<-dev.off();
}

