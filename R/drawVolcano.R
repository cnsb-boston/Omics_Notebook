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
#' @export
drawVolcano <- function(dat, type, subset_rows=F,
                        title_add="", top_rows="",
                        outputpath=output_contrast_path){ 
  # Gene labels
  label_names <- c(rownames(dat[dat$logFC>0,][order(dat$P.Value, decreasing=FALSE),])[1:6],
                   rownames(dat[dat$logFC<0,][order(dat$P.Value, decreasing=FALSE),])[1:6],
                   rownames(dat[order(dat$logFC, decreasing=FALSE),])[1:4],
                   rownames(dat[order(dat$logFC, decreasing=TRUE),])[1:4] )
  
  label_names <- unique(label_names)
  
  # Dot coloring
  dat$Significance <- "NS"
  dat$Significance[( (rownames(dat)%in%top_rows) & sign(dat$logFC)>0) ] <- "Up"
  dat$Significance[( (rownames(dat)%in%top_rows) & sign(dat$logFC)<0) ] <- "Down"
  dat$Signicance <- factor(dat$Significance, levels="NS", "Up", "Down")
  
  # Plot
  plot1 <- ggplot(data=data.frame(dat)) + 
    geom_point(aes(x=logFC, y=-log10(P.Value),colour=(dat$Significance)),size=0.7, pch=19) +
    scale_colour_manual(values=c(NS="grey", Up="Red", Down="blue"))   +
    labs(title=paste(type,"Volcano Plot\n ", title_add, sep=" ")) + theme_bw() + 
    theme(legend.title=element_blank())
  plot2 <- plot1 +
    { if( "Gene" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]),size=2, aes(x=logFC, y=-log10(P.Value), label=Gene))
      else if( "mz" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]),size=2, aes(x=logFC, y=-log10(P.Value), label=round(mz, digits=2)))
      else geom_text_repel(data=data.frame(dat[label_names,]), aes(x=logFC, y=-log10(P.Value),size=1, label=feature_identifier))
    } + geom_point(data=data.frame(dat[label_names,]), aes(x=logFC, y=-log10(P.Value) ),size=0.7, pch=21)
  
  output_filename <- file.path(outputpath, paste(type,"_volcano",".pdf", sep=''));
  pdf(output_filename, width=3, height=3);
  plot(plot1+theme(legend.position="none"))
  plot(plot2+theme(legend.position="none"))
  if(class(subset_rows)!="logical" ){
    for( j in 1:length(subset_rows)){ try({
      if(length(subset_rows[[j]])>100){subset_rows[[j]]<- subset_rows[[j]][1:100]}
      if(length(subset_rows[[j]])>0){
      if( "Gene" %in% colnames(dat)){
        plot(plot2+
               geom_text_repel(data=data.frame(dat[subset_rows[[j]],]),size=2,direction="y", aes(x=logFC, y=-log10(P.Value), label=Gene))+theme(legend.position="none") +
               geom_point(data=data.frame(dat[subset_rows[[j]],]), aes(x=logFC, y=-log10(P.Value) ),size=0.7, pch=21) ) 
      } else{
        plot(plot2 + 
               geom_text_repel(data=data.frame(dat[subset_rows[[j]],]),direction="y", aes(x=logFC, y=-log10(P.Value),size=1, label=feature_identifier))+theme(legend.position="none") +
               geom_point(data=data.frame(dat[subset_rows[[j]],]), aes(x=logFC, y=-log10(P.Value) ),size=0.7, pch=21) ) 
      }
      }
    }) }
  }
  grid.arrange(g_legend(plot1))
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
  dat$Signicance <- factor(dat$Significance, levels="NS", "Up", "Down")
  
  # Plot
  plot1 <- ggplot(data=data.frame(dat)) + geom_point(aes(x=Mean, y=logFC,colour=(dat$Significance)), size=0.8) +
    scale_colour_manual(values=c(NS="grey", Up="Red", Down="blue"))   +
    labs(title=paste(type,"\nMD Plot", sep="")) + theme_bw() + 
    theme(legend.title=element_blank())
  plot2 <- plot1 +
    { if( "Gene" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]),size=2, aes(x=Mean, y=logFC, label=Gene))
      else if( "mz" %in% colnames(dat)) geom_text_repel(data=data.frame(dat[label_names,]),size=2, aes(x=Mean, y=logFC, label=round(mz, digits=2)))
      else geom_text_repel(data=data.frame(dat[label_names,]),size=2, aes(x=Mean, y=logFC, label=feature_identifier))
    } + geom_point(data=data.frame(dat[label_names,]), aes(x=Mean, y=logFC ),size=0.8, pch=21)

output_filename <- file.path(outputpath, paste(type,"_MDplot",".pdf", sep=''));
pdf(output_filename, width=3, height=3);
plot(plot1+theme(legend.position="none"))
plot(plot2+theme(legend.position="none"))
if(class(subset_rows)!="logical"){
  for( j in 1:length(subset_rows)){ try({
    if(length(subset_rows[[j]])>100){subset_rows[[j]]<- subset_rows[[j]][1:100]}
    if(length(subset_rows[[j]])>0){
    if( "Gene" %in% colnames(dat)){
      plot(plot2+geom_text_repel(data=data.frame(dat[subset_rows[[j]],]),size=2,direction="x", aes(x=Mean, y=logFC, label=Gene))+theme(legend.position="none") )
    } else{
      plot(plot2+geom_text_repel(data=data.frame(dat[subset_rows[[j]],]),direction="x", aes(x=Mean, y=logFC,size=1, label=feature_identifier))+theme(legend.position="none") )
    }
    }
  }) }
}
grid.arrange(g_legend(plot1))
tmp<-dev.off();
}

