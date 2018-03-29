#-------------------------------------------------
#' Make Eset
#'
#' This function taes in a data frame that is a raw output from a proteomics search engine
#' (MaxQuant) and processes into an ExpressionSet object
#'
#' @param data A data.frame that contains sampes in rows, and intensity values and annotation data in columns
#' @param annotate An annotation data frame
#' @param type Type of data to be processed
#' @param cutoff Fraction of samples that must be non-zero to include feature in analysis
#' @param format Format of data
#'
#' @return an expression set object
#' 
#' @examples
#' eset.phospho <- makeEset(proteinGroups.phospho, annot, type='phos', cutoff=zero.percent);
#' 
#' @export
makeEset <- function(data, annotate, type='proteo', cutoff=zero.percent,
                     format="MaxQuant") {
  # remove contaminants and reverse hits
  data <- data[data[,"Potential.contaminant"]!='+',]; #remove potential contaminants
  data <- data[data[,"Reverse"]!='+',]; #remove reverse
  
  if (type=='proteo'){
    # pull out sample intensity values based on annotation
    annotate[,"Intensity"] <- make.names(annotate[,"Intensity"] );
    data.matrix <- as.matrix(data[, annotate[annotate[,"Intensity"]!="X","Intensity"] ]);
    # parse out protein namea
    data[,"Protein"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Majority.protein.IDs"])){ substr(x["Majority.protein.IDs"], 0, 
                                                        unlist(gregexpr(';', x["Majority.protein.IDs"]))[1]-1)
      } else {x["Majority.protein.IDs"]} } );
  }
  if(type=="phos"){
    data <- data[ data[,"Diagnostic.peak"]!='+',]; 
    #data <- data[ data[,"Localization.prob"]>=0.70 ,];
    
    # pull out sample intensity values based on annotation (added collapse multiple sites)
    annotate[,"Phospho.Intensity"] <- make.names(annotate[,"Phospho.Intensity"] );
    data1 <- data[, -grep(paste(annotate[annotate[,"Phospho.Intensity"]!="X","Phospho.Intensity"],
                                collapse='|'), colnames(data)) ]; 
    data1 <- data1[rep(rownames(data1), 3),];
    for (i in 1:length(annotate[annotate[,"Phospho.Intensity"]!="X","Phospho.Intensity"])){
      cols.keep <- colnames(data)[grep(annotate[annotate[,"Phospho.Intensity"]!="X","Phospho.Intensity"][i],
                                       colnames(data) )];  
      data1 <- cbind(data1, melt(data[, cols.keep], measure.vars=cols.keep,
                                 variable.name=paste("Phospho.Site.", i-1, sep=''),
                                 value.name=annotate[annotate[,"Phospho.Intensity"]!="X","Phospho.Intensity"][i] ));
    }
    data <- data1;
    data.matrix <- as.matrix(data[,annotate[annotate[,"Phospho.Intensity"]!="X","Phospho.Intensity"]]);
    # parse out protein namea
    data[,"Protein"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Protein"])){ substr(x["Protein"], 0, 
                                           unlist(gregexpr(';', x["Protein"]))[1]-1)}
      else {x["Protein"]} } );
  }
  
  # Make uniprot hyperlink
  data[,"Uniprot"] <- paste("<a href='https://www.uniprot.org/uniprot/",data[,"Protein"],"'>",
                            data[,"Protein"],"</a>", sep="")
  # Get gene names
  data[,"Gene"] <- apply(data, 1, function(x) {
    if(grepl(';', x["Gene.names"])){ substr(x["Gene.names"], 0, unlist(gregexpr(';', x["Gene.names"]))[1]-1)
    } else {x["Gene.names"]} } );
  rownames(data) <- make.unique(data[,"Gene"]); 
  
  # log 2 transform data
  data.matrix <- log2(data.matrix+1)
  
  # Make column and row names for data matrix
  if(type=="proteo"){
    colnames(data.matrix) <- annotate[annotate[,"Intensity"]!="X","SampleName"]
  } 
  if (type=="phos"){
    colnames(data.matrix) <- annotate[annotate[,"Phospho.Intensity"]!="X","SampleName"]
  }
  rownames(data.matrix) <- rownames(data)
  
  # make expression set object
  eset <- ExpressionSet(assayData=data.matrix)
  if (type=='proteo'){
    fData(eset) <- data[, -grep(paste(annotate[annotate[,"Intensity"]!="X","Intensity"],
                                      collapse='|'), colnames(data)) ];
    pData(eset) <- cbind(annotate[annotate[,"Intensity"]!="X",], colnames(exprs(eset)));
  }
  if(type=='phos'){
    fData(eset) <- data[, -grep(paste(annotate[annotate[,"Phospho.Intensity"]!="X","Phospho.Intensity"],
                                      collapse='|'), colnames(data)) ];
    pData(eset) <- cbind(annotate[annotate[,"Phospho.Intensity"]!="X",], colnames(exprs(eset)));
  }
  rownames(pData(eset)) <- colnames(data.matrix)
  
  # filter out features not detected in enough samples
  eset <- eset[ rowSums(exprs(eset)>0)>=cutoff*ncol(exprs(eset)), ];
  
  return(eset);
# Function for getting gene name when FASTA not properly configured
# data[,"gene.name"] <- apply(data, 1, function(x) substr(x["Majority.protein.IDs"],
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[2]+1,
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[3]-3) )
}






