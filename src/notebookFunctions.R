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

makeEset <- function(data, annotate, type='proteo', cutoff=0.05,
                     format="MaxQuant") {
  data <- data[data[,"Potential.contaminant"]!='+',]; #remove potential contaminants
  data <- data[data[,"Reverse"]!='+',]; #remove reverse
  annotate[,"SampleName"] <- make.names(annotate[,"SampleName"] ); #format sample names
  
  if (type=='proteo'){
    # pull out sample intensity values based on annotation
    annotate[,"Intensity"] <- make.names(annotate[,"Intensity"] );
    data.matrix <- as.matrix(data[, annotate[ (annotate[,"Intensity"]!="NA." & annotate[,"SampleName"]!="NA."),"Intensity"] ]);
    # parse out protein namea
    data[,"Protein"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Majority.protein.IDs"])){ substr(x["Majority.protein.IDs"], 0, 
                                                        unlist(gregexpr(';', x["Majority.protein.IDs"]))[1]-1)
      } else {x["Majority.protein.IDs"]} } );
  }
  if(type=="phos"){
    data <- data[ data[,"Diagnostic.peak"]!='+',]; #remove diagnostic features
    data <- data[ data[,"Localization.prob"]>=0.70 ,]; #filter low probability features
    
    # pull out sample intensity values based on annotation and collapse multiple sites
    annotate[,"Phospho.Intensity"] <- make.names(annotate[,"Phospho.Intensity"] );
    data1 <- data[, -grep(paste(annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Phospho.Intensity"],
                                collapse='|'), colnames(data)) ]; 
    data1 <- data1[rep(rownames(data1), 3),];
    for (i in 1:length(annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Phospho.Intensity"])){
      cols.keep <- colnames(data)[grep(annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Phospho.Intensity"][i],
                                       colnames(data) )];  
      data1 <- cbind(data1, melt(data[, cols.keep], measure.vars=cols.keep,
                                 variable.name=paste("Phospho.Site.", i-1, sep=''),
                                 value.name=annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Phospho.Intensity"][i] ));
    }
    data <- data1;
    data.matrix <- as.matrix(data[,annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Phospho.Intensity"]]);
    # parse out protein namea
    data[,"Protein"] <- apply(data, 1, function(x) {
      if(grepl(';', x["Protein"])){ substr(x["Protein"], 0, 
                                           unlist(gregexpr(';', x["Protein"]))[1]-1)}
      else {x["Protein"]} } );
  }
  
  # Make uniprot hyperlink
  data[, "Link"] <- paste("https://www.uniprot.org/uniprot/", data[,"Protein"], sep="")
  data[,"Uniprot"] <- paste("<a href='https://www.uniprot.org/uniprot/",data[,"Protein"],"'>",
                            data[,"Protein"],"</a>", sep="")
  # Get gene names
  data[,"Gene"] <- apply(data, 1, function(x) {
    if(grepl(';', x["Gene.names"])){ substr(x["Gene.names"], 0, unlist(gregexpr(';', x["Gene.names"]))[1]-1)
    } else {x["Gene.names"]} } );
  rownames(data) <- make.unique(data[,"Gene"]); 
  
  # Add uniprot annotation
  data <- cbind(data, addUniprotAnnotation(data[,"Protein"]))
  
  # log 2 transform data
  data.matrix <- log2(data.matrix+1)
  
  # Make column and row names for data matrix
  if(type=="proteo"){
    colnames(data.matrix) <- annotate[ (annotate[,"Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"SampleName"]
  } 
  if (type=="phos"){
    colnames(data.matrix) <- annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"SampleName"]
  }
  rownames(data.matrix) <- rownames(data)
  
  # make expression set object
  eset <- ExpressionSet(assayData=data.matrix)
  if (type=='proteo'){
    fData(eset) <- data[, -grep(paste(annotate[ (annotate[,"Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Intensity"],
                                      collapse='|'), colnames(data)) ];
    pData(eset) <- cbind(annotate[ (annotate[,"Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),], colnames(exprs(eset)));
  }
  if(type=='phos'){
    fData(eset) <- data[, -grep(paste(annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),"Phospho.Intensity"],
                                      collapse='|'), colnames(data)) ];
    pData(eset) <- cbind(annotate[ (annotate[,"Phospho.Intensity"]!="NA."& annotate[,"SampleName"]!="NA."),], colnames(exprs(eset)));
  }
  rownames(pData(eset)) <- colnames(data.matrix)
  
  # filter out features not detected in enough samples
  eset <- eset[ rowSums(exprs(eset)>0)>=cutoff*ncol(exprs(eset)), ];
  
  return(eset);
# Code for getting gene name when FASTA not properly configured:
# data[,"gene.name"] <- apply(data, 1, function(x) substr(x["Majority.protein.IDs"],
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[2]+1,
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[3]-3) )
}


#-------------------------------------------------
#' Add Uniprot Annotation
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

addUniprotAnnotation <- function(IDs){
  
  uniprot_columns <- c("comment(FUNCTION)", "comment(SUBCELLULAR LOCATION)", "comment(DISEASE)",
                     "go(biological process)", "go(molecular function)", "go(cellular component)", "go-id", "database(Reactome)", "database(KEGG)",
                     "database(BioCyc)", "database(Ensembl)", "database(ChEMBL)", "database(IntAct)","database(STRING)")
  uniprot_col_names <- c("Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                       "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID",
                       "BioCyc_ID", "Ensembl_ID", "ChEMBL_ID", "IntAct_ID","STRING_ID")

  annotUniprot <- data.frame(matrix(ncol=length(uniprot_col_names)+1))
  colnames(annotUniprot) <- c("ENTRY", uniprot_col_names)

  IDs_unique <- IDs[!duplicated(IDs)] #get unique IDs
  
  # Query uniprot server 100 entries at a time
  for (i in 1: (length(IDs_unique) %/% 100)  ){
    info_url<-paste("http://www.uniprot.org/uniprot/?format=tab&columns=id,",gsub(" ","%20", paste(uniprot_columns, collapse=",")),
                    "&query=accession%3A",paste(IDs_unique[ (1 + 100*(i-1)) : (100*i) ], collapse="+OR+accession%3A"), sep="")
    if(url.exists(info_url)==TRUE){
      invisible(info_annot <- try(data.frame(read.delim(url(info_url),header=TRUE, stringsAsFactors=FALSE))) )
      if (class(info_annot) != 'try-error') { if(dim(info_annot)[2]==length(uniprot_columns)+1){
        if(colnames(info_annot)[2]=="Function..CC."){ colnames(info_annot) <- colnames(annotUniprot)
                                                      annotUniprot <- rbind(annotUniprot, info_annot) }
      } }
    }
    Sys.sleep(2)
  }
  info_url<-paste("http://www.uniprot.org/uniprot/?format=tab&columns=id,",gsub(" ","%20", paste(uniprot_columns, collapse=",")),
                  "&query=accession%3A",paste(IDs_unique[ ((length(IDs) %/% 100)*100) : (((length(IDs_unique) %/% 100)*100) + (length(IDs_unique) %% 100))],
                                              collapse="+OR+accession%3A"), sep="")
  if(url.exists(info_url)==TRUE){
    invisible(info_annot <- try(data.frame(read.delim(url(info_url),header=TRUE, stringsAsFactors=FALSE))) )
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