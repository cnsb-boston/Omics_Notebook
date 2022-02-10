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

  colnames_annotUniprot <- c("ENTRY", uniprot_col_names)
  
  # List of uniprot IDs (remove duplicates for speed)
  IDs_unique <- IDs[!duplicated(IDs)] #get unique IDs

  id_groups=split(IDs_unique,factor((1:length(IDs_unique)) %/% 100))
  
  # Query uniprot server 100 entries at a time
  annotUniprot <- do.call("rbind",lapply(id_groups,FUN=function(q_ids){ try({
    ret <- NULL
    info_url<-paste0("https://www.uniprot.org/uniprot/?format=tab&columns=id,",gsub(" ","%20", paste(uniprot_columns, collapse=",")),
                    "&query=accession%3A",paste(q_ids, collapse="+OR+accession%3A"))
    if(RCurl::url.exists(info_url)==TRUE){
      invisible(info_annot <- try(data.frame(read.delim(url(info_url),header=TRUE, stringsAsFactors=FALSE, quote=""))) )
      if (class(info_annot) != 'try-error' &&
          dim(info_annot)[2]==length(uniprot_columns)+1 &&
          colnames(info_annot)[2]=="Function..CC."){
        colnames(info_annot) <- colnames_annotUniprot
        ret <- info_annot
      }
    }
    Sys.sleep(2)
    ret
  }) 
  }))
 
  # make data frame corresponding to original ID list
  annotatedUniprot <- data.frame(matrix(ncol=length(uniprot_col_names)+1, nrow=length(IDs)))
  colnames(annotatedUniprot) <- colnames_annotUniprot
  annotatedUniprot[,"ENTRY"] <- IDs
  for (r in 1:nrow(annotatedUniprot)){ if(sum(annotUniprot[,"ENTRY"]==annotatedUniprot[r,"ENTRY"])!=0) {
    annotatedUniprot[r,uniprot_col_names] <-annotUniprot[which(annotUniprot[,"ENTRY"]==annotatedUniprot[r,"ENTRY"])[1],uniprot_col_names]
  } }

  return (annotatedUniprot)
}


