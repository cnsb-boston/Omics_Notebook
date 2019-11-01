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
  for (i in 1: (length(IDs_unique) %/% 100)  ){ try({
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
  }) }
  try({
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
  })
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


