##-------------------------------------------------
#' Make Eset
#'
#' This function takes in a data frame that is a raw output from a proteomics search engine
#' (MaxQuant) or other Omics data and processes into an ExpressionSet object. 
#'
#' @param data A data.frame that contains sampes in rows, and intensity values and annotation data in columns
#' @param annotate An annotation data frame that specifies columns in data that are samples. See Omics_Notebook example for format.
#' @param type Type of data to be processed or a name for the Omics Data.
#' @param cutoff Fraction of samples that must be non-zero to include feature in analysis (filter for low detected features).
#' @param data_format Format of data. See Omics_Notebook example for options. 
#' @param uniprot_annotation Whether or not to query uniprot for annotation data. (Requires web access and time consuming).
#' @param log_transform TRUE to log transform data
#'
#' @return an expression set object
#' 
#' @examples
#' eset <- makeEset(proteinGroups, annot, type='phos');
#' 
#' @export

makeEset <- function(data, annotate, type, log_transform=TRUE,
                     data_format, uniprot_annotation=TRUE, outputpath=output_plots_path) {
  # For MaxQuant formats, remove potential contaminants and reverse hits
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
  # Get and format Sample Names from annotation
  annotate[,"SampleName"] <- make.names(annotate[,"SampleName"] ); #format sample names
  annotate[,type] <- make.names(annotate[,type] ); # data column headers
  samp_cols = annotate[ (annotate[,type]!="NA." & annotate[,"SampleName"]!="NA."),type]
  missing_cols = setdiff(samp_cols, colnames(data))
  if(!grepl("Sites..MQ.", data_format) && length(missing_cols)>0){
    stop(paste0("Annotation file requests columns not present in data file (",type,", / ",data_format,"): ", paste0(missing_cols, collapse=", ")))
  }

  split_semi = function(data, dest, src){
    if(src %in% colnames(data)){
      data[,dest] = sub(';.*', '', data[,src])
    } else if(!(dest %in% colnames(data))) {
      stop(paste0(data_format, "input file requires one of the following columns:", dest, src))
    }

    data
  }

  if (data_format=='Protein.Groups..MQ.'){
    # pull out sample intensity values based on annotation
    data.matrix <- as.matrix(data[, samp_cols]);
    # parse out protein names
    data = split_semi(data, "Protein", "Majority.protein.IDs")
  }
  if (data_format=='Peptides..MQ.'){
    # pull out sample intensity values based on annotation
    data.matrix <- as.matrix(data[, samp_cols]);
    # parse out protein names
    data = split_semi(data, "Protein", "Proteins")
  }
  # Get data for MaxQuant Sites format
  if( grepl("Sites..MQ.", data_format) ){
    if( "+" %in%  data[,"Diagnostic.peak"] ) { data <- data[ data[,"Diagnostic.peak"]!='+',]; } #remove diagnostic features
    data <- data[ data[,"Localization.prob"]>=0.70 ,]; #filter low probability features
    
    # pull out sample intensity values based on annotation and collapse multiple sites
    data1 <- data[, !(colnames(data) %in% samp_cols)];
    data1 <- data1[rep(rownames(data1), 3),];
    data2 = do.call("cbind",lapply(1:length(samp_cols), FUN=function(i){
      cols.keep <- grep(samp_cols[i], colnames(data), value=T);  
      reshape2::melt(data[, cols.keep], measure.vars=cols.keep,
                     variable.name=paste("Phospho.Site.", i-1, sep=''),
                     value.name=samp_cols[i]);
    }))
    data <- cbind(data1, data2)
    data.matrix <- as.matrix(data[, samp_cols]);
    # parse out protein names
    data = split_semi(data, "Protein", "Protein")
  }
 
  sample_column_names <- colnames(data[, which(colnames(data) %in% samp_cols) ])
  if( data_format=="Generic" ){
    data.matrix <- as.matrix(data[, sample_column_names ]);
  }
  
  if( grepl("Metabolites", data_format) ){ # Either OpenMS output or XCMS Online output
    data.matrix <- as.matrix(data[, samp_cols]);
    
    if( "mzmed" %in% colnames(data) ){ data[,"mz"] <- data[,"mzmed"] } # if xmcs output, make mz column
    if( "mz_cf" %in% colnames(data) ){ data[,"mz"] <- data[,"mz_cf"] } # if openms output, make mz column
    
    if( "rtmed" %in% colnames(data) ){ data[,"rt"] <- data[,"rtmed"] } # if xmcs output, make mz column
    if( "rt_cf" %in% colnames(data) ){ data[,"rt"] <- data[,"rt_cf"] } # if openms output, make mz column
    
    if( "links" %in% colnames(data) ) { # if there are CSIFingerID guesses, grab first one and list in column
      data[,"CSIFingerID_Pubmed"] <- str_match( data[,"links"], "PubChem:\\((.*?)[ )]")[,2]
      data[,"CSIFingerID_HMDB"] <- str_match( data[,"links"], "HMDB:\\((.*?)\\)")[,2]
      data[,"identifier"] <- data[,"CSIFingerID_HMDB"]
    }
    
  }
  
  # data.matrix <- as.numeric(as.character(data.matrix) )
  # data.matrix[is.na(data.matrix)] <- 0
  
  try({ drawRange(matrix=data.matrix, file_name = type, outputpath = outputpath ) })
  
  # log2 Intensity Values
  if(log_transform){ data.matrix <- log2(data.matrix+1); }
  
  if ( "Protein" %in% colnames(data) ){
    # Add uniprot annotation
    uprot = data$Protein
    if(uniprot_annotation == TRUE) { try({
      data <- cbind(data, getUniprotAnnotation(IDs=data[,"Protein"], genes=!("Gene.names" %in% colnames(data))))
      uprot = data$ENTRY
    }, silent=FALSE)}
    # Make uniprot hyperlink
    data[, "Link"] <- paste("https://www.uniprot.org/uniprot/", uprot, sep="")
    data[,"Uniprot"] <- paste("<a href='https://www.uniprot.org/uniprot/",uprot,"'>",
                              uprot,"</a>", sep="")
  }
  # Parse maxquant to get gene names
  if ( "Gene.names" %in% colnames(data) ){ 
    data = split_semi(data, "Gene", "Gene.names")
  } else if("Fasta.headers" %in% colnames(data)){ # For species that have Uniprot but not MQ gene annotations
    uni=grepl("^(sp|tr)\\|",data[,"Fasta.headers"])
    if(any(uni)){
      data[,"Gene"]=""
      data[uni,"Gene"]=sub("^..\\|.*\\|([^_]*).*","\\1",data[uni,"Fasta.headers"])
    }
  }
  # make feature identifiers/rownames
  if( "Gene" %in% colnames(data) ){
    if( "Protein" %in% colnames(data) ){
      if(("Amino.acid" %in% colnames(data))&("Position"%in%colnames(data))){
        data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"],"_", data[,"Protein"],"_",
                                                          data[,"Amino.acid"],"", data[,"Position"], sep="" ) );
      } else {
        data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"], data[,"Protein"], sep="_" ) );
      }
    } else if( "Transcript" %in% colnames(data) ){
      data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"], data[,"Transcript"], sep="_" ) );
    } else {
      data[,"feature_identifier"] <- make.unique( data[,"Gene"] );
    }
  } else {
    data[,"feature_identifier"] <- make.unique( make.names(data[,1]) ); # take the first column as feature ID if nothing else is found.
    #data[,"Gene"] <- data[,"feature_identifier"] ;
  }
  
  data.matrix[is.na(data.matrix)] <- 0
  
  # Make column and row names for data matrix
  colnames(data.matrix) <- annotate[which(annotate[,type] %in% sample_column_names), "SampleName"]
  rownames(data.matrix) <- data[,"feature_identifier"]
  
  # make expression set object
  eset <- ExpressionSet(assayData=data.matrix)

  fData(eset) <- data[, which(!(colnames(data) %in% sample_column_names)) ];
  pData(eset) <- cbind(annotate[ which(annotate[,type] %in% sample_column_names),],
                       colnames(exprs(eset)) );

  #if( !("Gene" %in% colnames(fData(eset))) ){ pData(eset)$Gene <- rownames(eset) }
  rownames(pData(eset)) <- colnames(data.matrix)
  
  return(eset);
# Code for getting gene name when FASTA not properly configured:
# data[,"gene.name"] <- apply(data, 1, function(x) substr(x["Majority.protein.IDs"],
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[2]+1,
#                                                         unlist(gregexpr('=', x["Majority.protein.IDs"]))[3]-3) )
}




