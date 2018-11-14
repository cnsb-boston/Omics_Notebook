#################################################
##
##           NOTEBOOK.RMD FUNCTIONS
##
#################################################

#-------------------------------------------------
#' Make Eset
#'
#' This function takes in a data frame that is a raw output from a proteomics search engine
#' (MaxQuant) and processes into an ExpressionSet object
#'
#' @param data A data.frame that contains sampes in rows, and intensity values and annotation data in columns
#' @param annotate An annotation data frame, including 
#' @param type Type of data to be processed
#' @param cutoff Fraction of samples that must be non-zero to include feature in analysis
#' @param data_format Format of data
#' @param uniprot_annotation Whether or not to query uniprot for annotation info
#' @param log_transform TRUE to log transform data
#'
#' @return an expression set object
#' 
#' @examples
#' eset <- makeEset(proteinGroups.phospho, annot, type='phos', cutoff=zero.percent);
#' 
#' @export

makeEset <- function(data, annotate, type, cutoff=0, log_transform=TRUE,
                     data_format, uniprot_annotation=TRUE) {
  
  if( grepl(".MQ.", data_format) ){
    if ( ("Potential.contaminant" %in% colnames(data)) ){
      if( "+" %in%  data[,"Potential.contaminant"] ) {
        data <- data[data[,"Potential.contaminant"]!='+',]; #remove potential contaminants
    } }
    if ("Reverse" %in% colnames(data)){
      if( "+" %in%  data[,"Reverse"] ) {
        data <- data[data[,"Reverse"]!='+',]; #remove reverse
    } }
  }
  
  annotate[,"SampleName"] <- make.names(annotate[,"SampleName"] ); #format sample names
  annotate[,type] <- make.names(annotate[,type] ); # data column headers
  
  if (data_format=='Protein.Groups..MQ.'){
    # pull out sample intensity values based on annotation
    data.matrix <- as.matrix(data[, annotate[ (annotate[,type]!="NA." & annotate[,"SampleName"]!="NA."),type] ]);
    # parse out protein namea
    data[,"Protein"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Majority.protein.IDs"])){ substr(x["Majority.protein.IDs"], 0, 
                                                        unlist(gregexpr(';', x["Majority.protein.IDs"]))[1]-1)
      } else {x["Majority.protein.IDs"]} } );
  }
  if( grepl("Sites..MQ.", data_format) ){
    if( "+" %in%  data[,"Diagnostic.peak"] ) { data <- data[ data[,"Diagnostic.peak"]!='+',]; } #remove diagnostic features
    data <- data[ data[,"Localization.prob"]>=0.70 ,]; #filter low probability features
    
    # pull out sample intensity values based on annotation and collapse multiple sites
    data1 <- data[, -grep(paste(annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),type],
                                collapse='|'), colnames(data)) ]; 
    data1 <- data1[rep(rownames(data1), 3),];
    for (i in 1:length(annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),type])){
      cols.keep <- colnames(data)[grep(annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),type][i],
                                       colnames(data) )];  
      data1 <- cbind(data1, melt(data[, cols.keep], measure.vars=cols.keep,
                                 variable.name=paste("Phospho.Site.", i-1, sep=''),
                                 value.name=annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),type][i] ));
    }
    data <- data1;
    data.matrix <- as.matrix(data[,annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),type]]);
    # parse out protein namea
    data[,"Protein"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Protein"])){ substr(x["Protein"], 0, 
                                           unlist(gregexpr(';', x["Protein"]))[1]-1)}
      else {x["Protein"]} } );
  }
  
  if( data_format=="Generic" ){
    data.matrix <- as.matrix(data[, annotate[ (annotate[,type]!="NA." & annotate[,"SampleName"]!="NA."),type] ]);

  }
  
  if( data_format=="Metabolites" ){ # Either OpenMS/Sirius/CSIFingerID output or XCMS Online output
    data.matrix <- as.matrix(data[, annotate[ (annotate[,type]!="NA." & annotate[,"SampleName"]!="NA."),type] ]);
    
    if( "mzmed" %in% colnames(data) ){ data[,"mz"] <- data[,"mzmed"] } # if xmcs output, make mz column
    
    if( "links" %in% colnames(data) ) { # if there are CSIFingerID guesses, grab first one and list in column
      data[,"CSIFingerID_Pubmed"] <- str_match( data[,"links"], "PubChem:\\((.*?)[ )]")[,2]
      data[,"CSIFingerID_HMDB"] <- str_match( data[,"links"], "HMDB:\\((.*?)\\)")[,2]
    }
    
  }
  
  # log2 Intensity Values
  if(log_transform){ data.matrix <- log2(data.matrix+1); }
  
  if ( "Protein" %in% colnames(data) ){
    # Make uniprot hyperlink
    data[, "Link"] <- paste("https://www.uniprot.org/uniprot/", data[,"Protein"], sep="")
    data[,"Uniprot"] <- paste("<a href='https://www.uniprot.org/uniprot/",data[,"Protein"],"'>",
                              data[,"Protein"],"</a>", sep="")
    # Add uniprot annotation
    if(uniprot_annotation == TRUE) { try({data <- cbind(data, getUniprotAnnotation(IDs=data[,"Protein"])) }, silent=TRUE) }
  }
  # Parse maxquant to get gene names
  if ( "Gene.names" %in% colnames(data) ){ 
    data[,"Gene"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Gene.names"])){ substr(x["Gene.names"], 0, unlist(gregexpr(';', x["Gene.names"]))[1]-1)
      } else {x["Gene.names"]} } );
  }
  # make feature identifiers/rownames
  if( "Gene" %in% colnames(data) ){
    
    if( "Protein" %in% colnames(data) ){
      data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"], data[,"Protein"], sep="_" ) );
    } else if( "Transcript" %in% colnames(data) ){
      data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"], data[,"Transcript"], sep="_" ) );
    } else {
      data[,"feature_identifier"] <- make.unique( data[,"Gene"] );
    }
    
  } else {
    data[,"feature_identifier"] <- make.unique( make.names(data[,1]) ); # take the first column as feature ID if nothing else is found.
    data[,"Gene"] <- data[,"feature_identifier"] ;
  }
  
  # Make column and row names for data matrix
  colnames(data.matrix) <- annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),"SampleName"]
  rownames(data.matrix) <- data[,"feature_identifier"]
  
  # make expression set object
  eset <- ExpressionSet(assayData=data.matrix)

  fData(eset) <- data[, -grep(paste(annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),type],
                                    collapse='|'), colnames(data)) ];
  pData(eset) <- cbind(annotate[ (annotate[,type]!="NA."& annotate[,"SampleName"]!="NA."),],
                       colnames(exprs(eset)));

  #if( !("Gene" %in% colnames(fData(eset))) ){ pData(eset)$Gene <- rownames(eset) }
  rownames(pData(eset)) <- colnames(data.matrix)
  
  # filter out features not detected in enough samples
  if( data_format!="Metabolites" ) { eset <- eset[ rowSums(exprs(eset)>0)>=cutoff*ncol(exprs(eset)), ]; }
  
  return(eset);
# Code for getting gene name when FASTA not properly configured:
# data[,"gene.name"] <- apply(data, 1, function(x) substr(x["Majority.protein.IDs"],
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[2]+1,
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[3]-3) )
}

#-------------------------------------------------
#' Get Uniprot Annotation
#'
#' This function takes in Uniprot ID's and uses the web API to add annotation,
#'
#' @param IDs A list of uniprot IDs
#'
#' @return a data frame with annotation inforamtion
#' 
#' @examples
#' 
#' @export

getUniprotAnnotation <- function(IDs){
  
  # Uniprot entries to fetch (and col names)
  uniprot_columns <- c("comment(FUNCTION)", "comment(SUBCELLULAR LOCATION)", "comment(DISEASE)",
                     "go(biological process)", "go(molecular function)", "go(cellular component)", "go-id", "database(Reactome)", "database(KEGG)",
                     "database(BioCyc)", "database(Ensembl)", "database(ChEMBL)", "database(IntAct)","database(STRING)")
  uniprot_col_names <- c("Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                       "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID",
                       "BioCyc_ID", "Ensembl_ID", "ChEMBL_ID", "IntAct_ID","STRING_ID")
  
  # create data frame for storing info
  annotUniprot <- data.frame(matrix(ncol=length(uniprot_col_names)+1))
  colnames(annotUniprot) <- c("ENTRY", uniprot_col_names)
  
  # List of uniprot IDs (remove duplicates for speed)
  IDs_unique <- IDs[!duplicated(IDs)] #get unique IDs
  
  # Query uniprot server 100 entries at a time
  for (i in 1: (length(IDs_unique) %/% 100)  ){
    info_url<-paste("https://www.uniprot.org/uniprot/?format=tab&columns=id,",gsub(" ","%20", paste(uniprot_columns, collapse=",")),
                    "&query=accession%3A",paste(IDs_unique[ (1 + 100*(i-1)) : (100*i) ], collapse="+OR+accession%3A"), sep="")
    if(url.exists(info_url)==TRUE){
      invisible(info_annot <- try(data.frame(read.delim(url(info_url),header=TRUE, stringsAsFactors=FALSE, quote=""))) )
      if (class(info_annot) != 'try-error') { if(dim(info_annot)[2]==length(uniprot_columns)+1){
        if(colnames(info_annot)[2]=="Function..CC."){ colnames(info_annot) <- colnames(annotUniprot)
                                                      annotUniprot <- rbind(annotUniprot, info_annot) }
      } }
    }
    Sys.sleep(2)
  }
  info_url<-paste("https://www.uniprot.org/uniprot/?format=tab&columns=id,",gsub(" ","%20", paste(uniprot_columns, collapse=",")),
                  "&query=accession%3A",paste(IDs_unique[ ((length(IDs) %/% 100)*100) : (((length(IDs_unique) %/% 100)*100) + (length(IDs_unique) %% 100))],
                                              collapse="+OR+accession%3A"), sep="")
  if(url.exists(info_url)==TRUE){
    invisible(info_annot <- try(data.frame(read.delim(url(info_url),header=TRUE, stringsAsFactors=FALSE, quote=""))) )
    if (class(info_annot) != 'try-error') {  if(dim(info_annot)[2]==length(uniprot_columns)+1){
      if(colnames(info_annot)[2]=="Function..CC."){ colnames(info_annot) <- colnames(annotUniprot)
                                                    annotUniprot <- rbind(annotUniprot, info_annot) }
    } }
  }
  annotUniprot <- annotUniprot[-1,] #remove first row (NAs)
  
  # make data frame corresponding to original ID list
  annotatedUniprot <- data.frame(matrix(ncol=length(uniprot_col_names)+1, nrow=length(IDs)))
  colnames(annotatedUniprot) <- colnames(annotUniprot)
  annotatedUniprot[,"ENTRY"] <- IDs
  for (r in 1:nrow(annotatedUniprot)){ if(sum(annotUniprot[,"ENTRY"]==annotatedUniprot[r,"ENTRY"])!=0) {
    annotatedUniprot[r,uniprot_col_names] <-annotUniprot[which(annotUniprot[,"ENTRY"]==annotatedUniprot[r,"ENTRY"])[1],uniprot_col_names]
  } }

  return (annotatedUniprot)
}

#-------------------------------------------------
#' Get HMDB Annotation
#'
#' This function takes in Uniprot ID's and uses the web API to add annotation,
#'
#' @param IDs A list of uniprot IDs
#'
#' @return a data frame with annotation inforamtion
#' 
#' @examples
#' 
#' @export

#getHMDBAnnotation <- function(IDs){



#################################################
##
##      NOTEBOOK_EXPL.RMD FUNCTIONS
##
#################################################

checkColor <- function(x){
  #return(all(x%in%colors() | grepl("^#(\\d|[a-f]){6,8}$",x,ignore.case=TRUE)))
  return(all(grepl("^#(\\d|[a-f]){6,8}$",x,ignore.case=TRUE)))
}


#-------------------------------------------------
#' Intensity Normalization
#'
#' This function takes proteomics results and performs normalization and QC plotting to inspect data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param annotate An annotation data frame
#' @param type Type of data to be processed
#' @param norm one of the supported normalization methods
#' @param outputpath output file path for plots
#'
#' @return an expression set object, with normalized data
#' 
#' @examples
#'   eset <- intensityNorm(eset, norm=norm_method)}
#' 
#' @export
intensityNorm <- function(eset, norm, type, outputpath=output_plots_path,
                          annotate=annot, cutoff=0, data_format ){
  
  col_palette <- rainbow(length(levels(as.factor(pData(eset)$Group))))
  annotCol <- col_palette[as.factor(pData(eset)$Group)]
  try({if("Color_Groups" %in% colnames(pData(eset))) {
    if( checkColor(pData(eset)[,"Color_Groups"]) ){
      annotCol <- pData(eset)[,"Color_Groups"]
    }
  } })

  if( data_format!="Metabolites" ) { eset <- eset[ rowSums(exprs(eset)>0)>=cutoff*ncol(exprs(eset)), ]; }
  
  eset_matrix <- exprs(eset)
  df <- data.frame(reshape2::melt(eset_matrix, id.vars = NULL));
  colnames(df) <- c("Feature","Sample", "Intensity");
  plot1 <- ggplot(df, aes(x = Sample, y = Intensity)) + geom_boxplot(fill=annotCol) +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) +  
    labs(title=paste(type, " log2 Intensity: Before Normalization", sep=''));
  plot2 <- ggplot(df, aes(x=Intensity, colour=Sample)) + geom_density()+ theme_bw() + 
    labs(title=paste(type, " log2 Intensity: Before Normalization", sep=''));
  
  if(norm=='quantile'){
    eset_matrix_norm <- as.matrix(normalize.quantiles(eset_matrix))
    dimnames(eset_matrix_norm) <- dimnames(eset_matrix)
  }
  else if (norm=='loess'){
    eset_matrix_norm <- as.matrix(normalizeCyclicLoess(eset_matrix, method='pairs'))
  }
  else {  eset_matrix_norm<-eset_matrix }
  
  if (norm != "none"){
    df <- data.frame(reshape2::melt(eset_matrix_norm, id.vars = NULL));
    colnames(df) <- c("Feature","Sample", "Intensity");
    plot3 <- ggplot(df, aes(x = Sample, y = Intensity)) + geom_boxplot(fill=annotCol) +
      theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
      labs(title=paste(type, " log2 Intensity: After ",
                       toupper(norm)," Normalization", sep=''));
    plot4 <- ggplot(df, aes(x=Intensity, colour=Sample)) + geom_density()+ theme_bw() + 
      labs(title=paste(type, " log2 Intensity: After ",
                       toupper(norm)," Normalization", sep=''));
  }
  
  output_filename<-file.path(outputpath, paste("boxplot_histogram_",type,".pdf", sep=''));
  pdf(output_filename);
  print(plot1);
  if(norm!="none") { print(plot3);}
  print(plot2);
  if(norm!="none") { print(plot4); }
  dev.off();
  
  if(norm!="none") { 
    grid.arrange(plot1+labs(title="Raw Data"), plot3+labs(title="Normalized"),
                 plot2+theme(legend.position="none")+labs(title="Raw Data"), 
                 plot4+theme(legend.position="none")+labs(title="Normalized"),
                 top=paste(type, ":", sep=""), ncol=4);
  } else {
    grid.arrange(plot1+labs(title="Intensity"),
                 plot2+theme(legend.position="none")+labs(title="Intensity"), 
                 top=paste(type, ":", sep=""), ncol=4);
    
  }
  
  output_filename<-file.path(outputpath, paste("MAplots_",type,".pdf", sep=''))
  pdf(output_filename)
  for (i in 1:ncol(eset)) {
    limma::plotMA(eset_matrix, array=i, main=paste(type, "MA Plot",i,sep=" "))
    if(norm!="none"){
      limma::plotMA(eset_matrix_norm, array=i, 
                  main=paste(type,toupper(norm),"Normalized MA Plot",i,sep=" "))
    }
  }
  dev.off();
  
  mn <- apply(eset_matrix_norm, 1, median);
  rle <- data.frame(sweep(eset_matrix_norm, MARGIN=1, STATS=mn, FUN='-'));
  output_filename <- file.path(outputpath, paste("RLE_",type,".pdf", sep=''));
  pdf(output_filename, width=11, height=8.5);
  par(mar=c(10,4,4,2)+0.1)
  boxplot(rle, main="RLE (Relative Log Expression)\nShould be centered on 0 (blue line)",
          names=colnames(eset_matrix_norm), las=2)
  lines(x=c(0,ncol(eset_matrix_norm)+1), y=rep(0,2), col="blue", lty=2)
  lines(x=c(0,ncol(eset_matrix_norm)+1), y=rep(0.1,2), col="red", lty=2)
  par(mar=c(5,4,4,2)+0.1)
  dev.off();
  
  eset_norm <- ExpressionSet(assayData=eset_matrix_norm);
  pData(eset_norm) <- pData(eset);
  rownames(pData(eset_norm)) <- colnames(eset_matrix_norm);
  fData(eset_norm) <- fData(eset);
  
  return (eset_norm);
}

#-------------------------------------------------
#' Draw PCA
#'
#' This function generates PCA plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return a pca plot
#' 
#' @examples
#'     if(isthereProteo) { plotProteo <- plot.pca(eset_global) }
#' 
#' @export
drawPCA <- function(eset, x_axis="PC1", y_axis="PC2", type, outputpath=output_plots_path, show_sample_names=TRUE) {
  
  PC_data <- prcomp(t(exprs(eset)))
  percent_variance <- summary(PC_data)$importance["Proportion of Variance",]*100
  
#   if ("Color_Groups" %in% colnames(pData(eset)) ){ 
#     colors_dots <- unique(pData(eset)[,"Color_Groups"])
#     names(colors_dots) <- unique(pData(eset)[,"Group"])
#   }
  
  pca_graph <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis], colour=pData(eset)$Group) )+ 
    geom_point() + theme_bw() + theme(legend.title=element_blank()) + 
    { if(show_sample_names==TRUE) geom_text_repel(aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis], label=rownames(PC_data$x)), colour="black") } +
#     { if ("Color_Groups" %in% colnames(pData(eset))) scale_color_manual(values=colors_dots)  } +
    labs(title=paste("PCA - ",type,sep=""),
       x=paste(x_axis, sprintf(" (%2.0f%%)", percent_variance[x_axis]), sep=""),
       y=paste(y_axis, sprintf(" (%2.0f%%)", percent_variance[y_axis]), sep=""))


  pca_loading_graph <- ggplot(data=as.data.frame(PC_data$rotation), aes(x=(PC_data$rotation[,x_axis]), y=PC_data$rotation[,y_axis])) +
    geom_text(label=rownames(PC_data$rotation), colour="black") + 
    labs(x=paste(x_axis, sep=""), y=paste(y_axis, sep=""),
         title=paste("PC Factor Loadings\n ",type, sep=""))

  output_filename <- file.path(outputpath, paste("PCAplots_",type,".pdf",sep=""))
  pdf(output_filename, width=4, height=3)
  print(pca_graph)
  print(pca_loading_graph)
  plot(PC_data, type = 'l', main=paste("PCs vs. Variance",type,sep=""))
  dev.off()

  return(pca_graph)
}

#-------------------------------------------------
#' Draw Venn Diagram
#'
#' This function generates venn diagrams
#' 
#' @param item_list a list of item lists
#' @param item_name a name for the items
#' @param outputpath output file path for plots
#' 
#' @examples
#' 
#' @export
drawVenn <- function(item_list, item_name, outputpath=output_plots_path){
  fill_col <- rainbow(length(item_list))
  futile.logger::flog.threshold(ERROR);
  venn <- venn.diagram( x=item_list, filename=NULL,lty="blank",# height=2000, width=2000, 
                        #cat.default.pos="outer", 
                        fill=fill_col,main=paste("Overlap: ", item_name, sep="") );
  output_filename <- file.path(outputpath, paste("VennDiagram_", item_name,".pdf", sep="") );
  pdf(output_filename, width=4, height=3);
  upset(fromList(item_list), order.by = "freq");
  grid.newpage();
  pushViewport(viewport(width=unit(0.6, "npc"), height=unit(0.8, "npc")));
  grid.draw(venn);
  tmp<-dev.off();  
}

drawXYCorr <- function(item_list, item_name, outputpath=output_plots_path, file_name=""){
  
  output_filename<-file.path(outputpath,paste("Correlation_Plots_",item_name,file_name,".pdf", sep=""))
  pdf(output_filename, width=4, height=3 )
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
        geom_point() + scale_color_viridis() + labs(x=names(item_list)[j], y=names(item_list)[k]) +
        labs(title=paste("Avg. ",item_name,": ",names(item_list)[j],"\n vs. ",names(item_list)[k],sep="" )) +
        theme_bw() + annotate("text",x=-Inf, y=Inf, label=label, color="black", vjust=1.5, hjust=-0.4 , parse=TRUE)
      print(plot)
    }
  }
  dev.off()
}


#-------------------------------------------------
#' Variation Plot and Filter
#'
#' This function generates plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return a matrix of most variable features by MAD
#' 
#' @examples
#' 
#' @export
variationPlot <- function(eset, type, outputpath=output_plots_path, 
                                percent_choice=0.1) {
  emat<-exprs(eset)
  MEAN <- apply(emat, 1, mean)
  STDEV <- apply(emat, 1, sd)
  MAD <- apply(emat, 1, mad)
  MED <- apply(emat, 1, median)
  top_hits <- names(MAD)[order(MAD, decreasing=TRUE)[1:(percent_choice*length(MAD))]]
  emat_top <- emat[top_hits[top_hits!=""],]
  
  output_filename<-file.path(outputpath, paste("variation_",type,".pdf", sep=""))
  pdf(output_filename)
  plot(MEAN, STDEV, pch=".", cex=4, main=paste("mean vs. stdev: ", type, sep=""))
  plot(MED, MAD, pch=19, cex=0.5, log="", main=paste("median vs MAD: ",type, sep=""))
  points(MED[top_hits], MAD[top_hits], pch=19, cex=0.5, col="red")
  legend("topright", pch=20, col=c("black", "red"), legend=c("all points", paste("filtered top ", percent_choice*100, "%", sep="")) )
  dev.off()
  
  output_filename<-file.path(outputpath, paste("correlation_",type,".pdf", sep=""))
  pdf(output_filename)
  corrplot::corrplot(cor(emat_top), type="upper", order="hclust", tl.col="black")
  dev.off()
  
  output_filename<-file.path(outputpath, paste("pairs_",type,".jpeg", sep=""))
  jpeg(output_filename)
  pairs(emat, pch=".")
  dev.off()
  return (top_hits[top_hits!=""])
}
#-------------------------------------------------
#' Draw the Heatmaps
#'
#' This function generates heatmap plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#' @param limmaSig is list of significant genes
#'
#' @return a matrix of most variable features by MAD
#' 
#' @examples
#' 
#' @export
drawHeatmaps <- function(eset, emat_top=FALSE, type, outputpath=output_plots_path,
                        mapcolor=map_color, subset=FALSE, k_clust=0, 
                        outputcontrastpath=output_contrast_path, limmaSig=FALSE, 
                        show_row_names=FALSE){
  # Annotation Colors
  annotLab <- data.frame(Group = factor(pData(eset)$Group, levels=unique(pData(eset)$Group)));
  annotCol <- list(Group = rainbow(length(levels(factor(pData(eset)$Group)))))
  try({ if("Color_Groups" %in% colnames(pData(eset))) {
    if( checkColor(pData(eset)[,"Color_Groups"]) ){
      annotCol <- list( Group=unique(pData(eset)[,"Color_Groups"])[1:length(levels(factor(pData(eset)$Group)))] )
    }
  } })
  names(annotCol$Group) <- levels(factor(annotLab$Group))
  ha_column <- HeatmapAnnotation(df=annotLab, col=annotCol)
  
  if(mapcolor=="viridis"){mapcolor <- viridis(11); maponeway <- viridis(11);
  } else {mapcolor <- (rev(brewer.pal(11, mapcolor)));
          maponeway <- rev(brewer.pal(9, "Blues")) }
  
  # Standard Heatmaps
  if(class(limmaSig)=="logical"){
    output_filename <- file.path(outputpath, paste("heatmaps_",type,".pdf", sep=''));
    pdf(output_filename);
    
    # all features, z-score
    emat_sel <- na.omit(t(scale(t(exprs(eset))))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2
    ht1 <- Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
                   row_names_gp=gpar(fontsize=4),
                   column_title=paste(type, ": All features, row z score", sep='') )
    draw(ht1)
    #Optional k clustering
    if(k_clust !=0){
      kclus <- kmeans(emat_sel, 3);
      split <- paste0("Cluster ", kclus$cluster)
      ht1 <-Heatmap(matrix=(emat_sel), col=mapcolor, name="", top_annotation=ha_column,cluster_columns=FALSE,
                    show_row_names=show_row_names,row_names_gp=gpar(fontsize=4), 
                    split=split,#km=k_clust,
                    column_title=paste(type, ": All features, row z score", sep='') )
      draw(ht1)
    }
    # Correlation
    #draw(Heatmap(matrix=cor(emat_sel), col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
    #             column_title=paste(type, ": Correlation, row z score", sep='') ))
  
    # Subset features only
    if( class(subset)!="logical" ){
      emat_sel <- exprs(eset)[subset,]
      emat_sel <- t(scale(t(emat_sel))) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      draw(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,cluster_columns=FALSE,
                   row_names_gp=gpar(fontsize=4),
                   column_title=paste(type, ": Select genes, row z score", sep='') ))
    }
  
    # variation filter, z score
    if ( class(emat_top) !="logical" ){
      emat_sel <- exprs(eset)[emat_top,]
      emat_sel <- na.omit(t(scale(t(emat_sel)))) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      draw(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
                   column_title=paste(type, ": Highest variation, row z score", sep='') ))
    }
  
    # all features, log2 intensity
    emat_sel <- na.omit(exprs(eset))
    draw(Heatmap(matrix=emat_sel, col=maponeway, name="", top_annotation=ha_column,show_row_names=FALSE,
                 column_title=paste(type, ": All features, log2 Intensity", sep='') ))
  
    tmp<-dev.off();
    return(ht1);
    
  } else if(class(limmaSig)!="logical"){
  # Differential Heatmap
    output_filename <- file.path(outputcontrastpath, paste("heatmaps_",type,".pdf", sep=''));
    pdf(output_filename);

    # limma differential expression, z score
    emat_sel <- exprs(eset[rownames(eset) %in% limmaSig,])
    emat_sel <- na.omit(t(scale(t(emat_sel)))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2 
    draw(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,
                 show_row_names=show_row_names,
                 row_names_gp=gpar(fontsize=4),
                 column_title=paste(type, ": Differential Expression, row z score", sep='') ))
    tmp<-dev.off();
  }
}

#-------------------------------------------------
#' Make Interactive Heatmap
#'
#' This function generates an interactive heatmap
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#' @mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#'
#' @return
#' 
#' @examples
#' 
#' @export
interactiveHeatmap <- function(eset, type, outputpath=output_plots_path,mapcolor=map_color ){
  
  annotCol <- c("red", "green", "blue", "yellow", "green", "purple",
                "brown", "black", "grey", "orange", "white", "light green", 
                "pink", "light blue" );
  sampleCols <- annotCol[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
  
  if(mapcolor=="viridis"){mapcolor <- viridis(11)
  } else {mapcolor <- colorRampPalette(rev(brewer.pal(11, mapcolor))) }
  
  #   emat.sel <- exprs(eset[rownames(eset) %in% rownames(limmaSig),])
  output_filename <- file.path(outputpath, paste("heatmap_all_",type,".html", sep=''));
  emat_sel <- exprs(eset)
  emat_sel <- na.omit(t(scale(t(emat_sel)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  suppressWarnings(invisible(capture.output(tmp<-heatmaply(emat_sel, scale='none', margins=c(100,100,40,20),
                                          col_side_colors=pData(eset)$Group,
                                          colors=mapcolor,file=output_filename ))) );
  rm(tmp)
}
#-------------------------------------------------
#' Save Files
#'
#' This function saves txt files for subsequent use
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path
#' @param workingdir specifies the working directory
#' @param saveTheEset specifies whether or not to save the eset as an RDS object
#' @param limmaRes topTable output
#' @param contrast_name name of the contrast
#' @param outputcontrastpath output file path
#'
#' @return
#' 
#' @examples
#' 
#' @export
saveFiles <- function(data, type, outputpath=output_files_path, workingdir=working_dir,
                      limmaRes=FALSE, contrast_name=FALSE, outputcontrastpath=output_contrast_path){
    
  output_filename <- file.path(outputpath, paste("Data_",type,".RDS", sep=""));             
  saveRDS(data, file=output_filename); 
  
  if ( "Gene" %in% colnames(fData(data[["eSet"]])) &  "Protein.names" %in% colnames(fData(data[["eSet"]])) ) {
    output_table <- cbind(fData(data[["eSet"]])[,c("Gene", "Protein.names")],exprs(data[["eSet"]]))
  } else {  output_table <- cbind( rownames(exprs(data[["eSet"]])), exprs(data[["eSet"]])) }
  output_filename <- file.path(outputpath,paste("Expression_matrix_",type,".txt", sep=''));
  write.table(x= output_table, file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
  
  if ( class(limmaRes) != "logical") {
    
    if("Gene"%in% colnames(limmaRes)){
      output_filename <- file.path(outputcontrastpath,paste("GSEA_",type, gsub("-","_",contrast_name),".rnk", sep=''));
      ranked <- cbind(limmaRes[,"Gene"],sign(limmaRes[,"logFC"]) * -log10(limmaRes[,"adj.P.Val"]))
      colnames(ranked)<-c("GeneName", "rank")
      colnames(ranked)<-c("GeneName", "rank")
      ranked <- ranked[order(as.numeric(ranked[,"rank"]), decreasing=TRUE),]
      ranked <- ranked[ranked[,"GeneName"]!="", ]
      ranked <- ranked[!duplicated(ranked[,"GeneName"]),]
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
    }
  }
}

#-------------------------------------------------
#' Write data to sheets
#'
#' This function outputs summary data to an excel sheet to share with collaborators
#' 
#' @param wb an openxlsx workbook object
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param data_format format of data to be processed
#' @param mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#' @param limmaFit a lmFit object
#' @param contrastnames specifies names of the contrasts/coefficients in the lmFit object
#'
#' @return wb with summary data formatted
#' 
#' @examples
#' 
#' @export
writeDataToSheets <- function(wb, eset, limmaFit=FALSE, data_format, mapcolor=map_color, type){
  
  # Make colors for sample names
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = gsub('.{2}$','', rainbow(length(levels(factor(pData(eset)$Group)))) ) )
  sampleCols <- annotCol$Group[1:length(levels(factor(pData(eset)$Group)))][factor(pData(eset)$Group)];
  mapcolor <- rev(brewer.pal(7, "RdYlBu"))
  
  if (class(limmaFit)!="logical"){
    # Create table for each coefficient
    DiffList <- vector("list", (ncol(limmaFit$coefficients)) )
    for(i in 1:length(DiffList) ) {
      DiffList[[i]] = topTable(limmaFit, adjust.method="BH", n=Inf, sort.by='p', coef=i)[,c("P.Value","adj.P.Val","logFC")]
    }
    # Match the row order to the first contrast
    eset <- eset[match(rownames(DiffList[[1]]),rownames(eset)),]
  } else { DiffList <- vector("list", 0 ) }
  
  # Create new sheet to add to the workbook
  stName <- as.character(type)
  addWorksheet(wb=wb, sheetName=stName)
  
  # Format Uniprot hyperlinks
  if ("Link" %in% colnames(fData(eset)) ){
    links <- fData(eset)$Link
    names(links) <- fData(eset)$Protein
    class(links) <- "hyperlink"
  }
  
  # Make row z-score values for "heatmap"
  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  
  # Create the data table
  formatted_table<-""
  if ( "feature_identifier" %in% colnames(fData(eset)) ) {  
    formatted_table <- data.frame( fData(eset)[,c("feature_identifier")],emat_sel )
    names <- c("Feature", colnames(eset)  )
  } else if ( "Gene" %in% colnames(fData(eset)) ) {  
    formatted_table <- data.frame( fData(eset)[,c("Gene")],emat_sel )
    names <- c("Gene", colnames(eset)  )
  } else { 
    formatted_table <- data.frame(  rownames(exprs(eset)),emat_sel )
    names <- c("Feature", colnames(eset)  )
  }
  col_widths <- c(20, rep(4, ncol(eset)));
  
  if ( "Gene" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Gene")] );
    names <- c(names, "Gene");
    col_widths<-c(col_widths, 16);
  }
  if ( "Protein.names" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Protein.names")] );
    names <- c(names, "Protein");
    col_widths<-c(col_widths, 50);
  }
  if ( "Protein" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Protein")] );
    names <- c(names, "Uniprot"); 
    col_widths<-c(col_widths, 16);
  }
  
  if ( "mz" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("mz")] );
    names <- c(names, "MZ"); 
    col_widths<-c(col_widths, 16);
  }
  
  col_index <- length(names)
  
  if(class(limmaFit)=="MArrayLM") { 
    for (i in 1:length(DiffList)) {
      DiffList[[i]] <- DiffList[[i]][match(rownames(DiffList[[1]]),rownames(DiffList[[i]])),];
      formatted_table <- data.frame(formatted_table, DiffList[[i]][,c("P.Value","adj.P.Val","logFC")] );
      names <- c(names, "P.Value","adj.P.Val","logFC");
      col_widths<-c(col_widths, 8,8,8);
    } 
  }
  
  if ( "Uniprot_Function" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Uniprot_Function")]); 
    names <- c(names,"Uniprot_Function"); col_widths<-c(col_widths, 50);
  }
  if ( "Uniprot_Cellular_Location" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Uniprot_Cellular_Location")] ); 
    names <- c(names,"Uniprot_Cellular_Location"); col_widths<-c(col_widths, 50);
  }
  if ( "Uniprot_Disease" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Uniprot_Disease")] );
    names <- c(names,"Uniprot_Disease"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_biological_process" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_biological_process")] ); 
    names <- c(names,"GO_biological_process"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_molecular_function" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_molecular_function")] ); 
    names <- c(names,"GO_molecular_function"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_cellular_component" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_cellular_component")] ); 
    names <- c(names,"GO_cellular_component"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_ID" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_ID")] ); 
    names <- c(names,"GO_ID"); col_widths<-c(col_widths, 8);
  }
  if ( "ReactomeID" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("ReactomeID")] ); 
    names <- c(names,"ReactomeID"); col_widths<-c(col_widths, 8);
  } 
  if ( "KEGG_ID" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("KEGG_ID")] ); 
    names <- c(names,"KEGG_ID"); col_widths<-c(col_widths, 8);
  }

  formatted_table<- data.frame(formatted_table, exprs(eset) );  
  names <- c(names,colnames(eset) );
  col_widths <- c(col_widths, rep(4, ncol(eset)));
  
  # Add phospho site probabilities and data
  try({ 
    if("Sites..MQ." %in% data_format) { 
      formatted_table <- data.frame(formatted_table, fData(eset)[,"Localization.prob"],
                                          fData(eset)[,grep("Probabilities", colnames(fData(eset)))],
                                          paste(fData(eset)[,"Amino.acid"],fData(eset)[,"Position"], sep=""),
                                          fData(eset)[,"Sequence.window"] ) 
      names <- c(names, "Localization Probability","Site Probabilities", "Amino Acid", "Peptide Sequence")
      col_widths <- c(col_widths, 16, 16, 16, 16)
    }
  })
  
  formatted_table <- data.frame(formatted_table, stringsAsFactors=FALSE)
  colnames(formatted_table) <-  make.unique(toupper(names));
  
  # Write data table to sheet
  writeDataTable(wb=wb, sheet=stName, x=formatted_table, xy=c("A",2), keepNA=FALSE, tableStyle="TableStyleLight1")
  if ("Link" %in% colnames(fData(eset)) ){ writeData(wb=wb, sheet=stName, x=links, xy=c( 1+(ncol(eset)+3 ) ,3)) }
  
  # Add heatmap color
  conditionalFormatting(wb=wb, sheet=stName, type="colourScale", cols=2:(1+ncol(eset)), rows=3:(2+nrow(eset)), style=mapcolor[c(1,4,7)])

  # Add color bars for fold change
  if(length(DiffList)>0) { i<-1; #for(i in 1:length(DiffList)){
    conditionalFormatting(wb=wb, sheet=stName, cols= ( 3+col_index+(3*(i-1))) , rows=3:(3+nrow(limmaFit)), type="databar", style=c("blue", "red"), showvalue=FALSE)
  } #}
  
  # Rotate text for sample names
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"),
             rows=2, cols=c+1)
  }
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"),
             rows=2, cols=(c + ( length(names) - ncol(eset) )  ) )
  }

  # Merge cells and add Intensity title
  mergeCells(wb=wb, sheet=stName, rows=1, cols=2:(1+ncol(eset)) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity, Row Z Score", xy=c(2,1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=2, stack=TRUE)
  
  # Add final intensity title
  mergeCells(wb=wb, sheet=stName, rows=1, cols=( length(names) : ( length(names) - ncol(eset) + 1 )) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity", xy=c(( length(names) - ncol(eset) + 1 ),1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=( length(names) - ncol(eset) + 1 ), stack=TRUE)
  
  # Add contrast titles
  if(length(DiffList)>0) { for(i in 1:length(DiffList)){
    startCol<- col_index + 1 + (3*(i-1))
    mergeCells(wb=wb, sheet=stName, rows=1, cols=startCol:(startCol+2) )
    writeData(wb=wb, sheet=stName, x=colnames(limmaFit$coefficients)[i], xy=c(startCol,1))
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=startCol, stack=TRUE)
  } }
  
  # Freeze columns/rows and bold column titles
  freezePane(wb=wb, sheet=stName, firstActiveRow=3, firstActiveCol=2)
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", fontColour='black'),
           rows=1:2, cols=1:length(names),gridExpand=TRUE, stack=TRUE)   
  
  # Don't treat gene names as dates
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(numFmt="TEXT"), rows=3:(2+nrow(eset)), cols=1, gridExpand=TRUE )
  
  # Set row heights and column widths
  setRowHeights(wb=wb, sheet=stName, rows=2, heights=100)
  setColWidths(wb=wb, sheet=stName, cols=1:ncol(formatted_table), widths= col_widths )
}




#################################################
##
##      NOTEBOOK_2_Diff.RMD FUNCTIONS
##
#################################################

#-------------------------------------------------
#'Draw Volcano Plots
#'
#' This function makes a volcano plot
#' 
#' @param dat is the result of topTable()
#' @param dt is the result of decideTests()
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#'  
#' 
#' @export
drawVolcano <- function(dat, type, outputpath=output_contrast_path){ 
  threshold_fc_high <- sort(dat$logFC, decreasing=TRUE)[8]
  threshold_fc_low <- sort(dat$logFC, decreasing=FALSE)[8]
  threshold_p_low <- sort(dat$P.Value, decreasing=FALSE)[10]
  
  dat$Significance <- "NS"
  dat$Significance[(dat$P.Value< 0.05 & sign(dat$logFC)>0) ] <- "Up"
  dat$Significance[(dat$P.Value< 0.05 & sign(dat$logFC)<0) ] <- "Down"
  dat$Signicance <- factor(dat$Significance, levels="NS", "Up", "Down")
  
  plot1 <- ggplot(data=data.frame(dat)) + geom_point(aes(x=logFC, y=-log10(P.Value),colour=(dat$Significance))) +
    scale_colour_manual(values=c(NS="grey", Up="Red", Down="blue"))   +
    labs(title=paste(type,"Volcano Plot", sep=" ")) + theme_bw() + 
    theme(legend.title=element_blank()) +
    geom_text_repel(data=data.frame(subset(dat, logFC>threshold_fc_high)), aes(x=logFC, y=-log10(P.Value), label=Gene)) + 
    geom_text_repel(data=data.frame(subset(dat, logFC<threshold_fc_low)), aes(x=logFC, y=-log10(P.Value), label=Gene)) + 
    geom_text_repel(data=data.frame(subset(dat, P.Value<threshold_p_low)), aes(x=logFC, y=-log10(P.Value), label=Gene))
  
  output_filename <- file.path(outputpath, paste("volcano_",type,".pdf", sep=''));
  pdf(output_filename);
  plot(plot1)
  tmp<-dev.off();
}

#-------------------------------------------------
#' Make Interactive Volcano
#'
#' This function makes a volcano plot
#' 
#' @param eset is an expression set object
#' @param dt is the result of decideTests()
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#'  
#' 
#' @export
interactiveVolcano <- function(eset, fit, limmaSig, dt, type, outputcontrastpath=output_contrast_path, col ){
  
  annotCol <- c("red", "green", "blue", "yellow", "green", "purple",
                "brown", "black", "grey", "orange", "white", "light green", 
                "pink", "light blue" );
  sampleCols <- annotCol[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
  
  #   if(!newcontrast){
  #     glMDSPlot(exprs(eset), groups=pData(eset)$Group, html=paste("MDS-Plot_",type,sep=""),
  #            folder="InteractivePlots", path=outputpath, launch=FALSE  );
  #   }
  #   glMDPlot(fit, groups=pData(eset)$Group, status=dt, counts=exprs(eset),
  #            anno=fit$genes[,c("Gene", "Protein", "Protein.names","Uniprot")],
  #            xlab="Normalized log2 Intensity",
  #            display.columns=c("Gene", "Protein", "Protein.names", "Uniprot"),
  #            side.main="Gene", side.ylab="Normalized Intensity",
  #            main=colnames(dt), sample.cols=sampleCols,
  #            html=paste("MD-Plot_",type,sep=""),
  #            folder="InteractivePlots", path=outputcontrastpath, launch=FALSE );
  annot_columns <- "feature_identifier";
  if("Gene" %in% colnames(fit$genes)) { annot_columns <- c(annot_columns, "Gene"); }
  if("Protein" %in% colnames(fit$genes)) { annot_columns <- c(annot_columns, "Protein"); }
  if("Protein.names" %in% colnames(fit$genes)) { annot_columns <- c(annot_columns, "Protein.names"); }
  if("Uniprot" %in% colnames(fit$genes)) { annot_columns <- c(annot_columns, "Uniprot"); }
  
  glXYPlot(x=fit$coef[,col], y=fit$lod[,col], xlab="logFC", status=dt[,col], counts=exprs(eset),
           groups=pData(eset)$Group, ylab="logOdds", #sample.cols=pData(eset)$Group,
           side.ylab="Normalized Intensity",
           anno=fit$genes[,annot_columns], side.main='feature_identifier',
           html=paste("Volcano-Plot_",type, sep=""),
           folder="InteractivePlots", path=outputcontrastpath, launch=FALSE  );
}

#################################################
##
##      NOTEBOOK_3_Enrch.RMD FUNCTIONS
##
#################################################

#-------------------------------------------------
#' Wrapper fnction for EnrichR
#'
#' This runs enrichment analysis with enrichr
#' 
#' @param genes, list of genes
#' @param type Type of data to be processed
#' @param search_data list of enrichr databases to search
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#' 
#' @export
runEnrichR <- function (genes, type, search_dat=search_databases, 
                        outputpath=output_contrast_path) {
  
  genes_up <- genes[genes[,"logFC"]>0,"Gene"]
  genes_down <- genes[genes[,"logFC"]<0,"Gene"]
  
  if(length(genes_up)>0){
    enriched_up <- enrichr(unique(genes_up), databases = search_dat);
    enriched_up <- bind_rows(enriched_up, .id="databases")
    enriched_up <- enriched_up[order(enriched_up[,"Adjusted.P.value"]),]
    output_filename <- file.path(outputpath,paste("enrichr_results_", type,"_up.txt", sep='')); 
    write.table(enriched_up, file=output_filename, sep="\t", row.names=FALSE);
  }

  if(length(genes_down)>0){
    enriched_down <- enrichr(unique(genes_down), databases = search_dat);
    enriched_down <- bind_rows(enriched_down, .id="databases")
    enriched_down <- enriched_down[order(enriched_down[,"Adjusted.P.value"]),]
    output_filename <- file.path(outputpath,paste("enrichr_results_", type,"_down.txt", sep='')); 
    write.table(enriched_down, file=output_filename, sep="\t", row.names=FALSE);
  }
}



