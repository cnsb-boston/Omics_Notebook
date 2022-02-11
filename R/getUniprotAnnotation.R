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
  uniprot_columns <- c("id","comment(FUNCTION)", "comment(SUBCELLULAR LOCATION)", "comment(DISEASE)",
                     "go(biological process)", "go(molecular function)", "go(cellular component)", "go-id", "database(Reactome)", "database(KEGG)",
                     "database(BioCyc)", "database(Ensembl)", "database(ChEMBL)", "database(IntAct)","database(STRING)")
  uniprot_col_names <- c("ENTRY","Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                       "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID",
                       "BioCyc_ID", "Ensembl_ID", "ChEMBL_ID", "IntAct_ID","STRING_ID")

  colnames_annotUniprot <- c("ENTRY", uniprot_col_names)
  
  # List of uniprot IDs (remove duplicates for speed)
  IDs_unique <- IDs[!duplicated(IDs)] #get unique IDs

  id_groups=split(IDs_unique,factor((1:length(IDs_unique)) %/% 5000))
  
  # Query uniprot server in batches
  annotUniprot <- do.call("rbind",lapply(id_groups,FUN=function(q_ids){ try({
    ret <- NULL
    res = httr::POST("https://www.uniprot.org/uploadlists/",
                     body=list(query=paste(q_ids,collapse=","), format="tab", from="ACC", to="ACC",
                               columns=paste0(uniprot_columns,collapse=",")))
    info_annot <- try(data.frame(read.delim(textConnection(httr::content(res)),
                                            header=TRUE, stringsAsFactors=FALSE, quote="")))
    if (class(info_annot) != 'try-error' &&
        ncol(info_annot)==length(uniprot_columns)+1 &&
        colnames(info_annot)[2]=="Function..CC."){
      info_annot=info_annot[,-ncol(info_annot)] # yourlist...
      colnames(info_annot) <- uniprot_col_names
      ret <- info_annot
    }
    Sys.sleep(2)
    ret
  }) 
  }))
 
  annotatedUniprot <- data.frame(ENTRY=IDs)
  annotatedUniprot <- merge(annotatedUniprot, annotUniprot, by="ENTRY", all.x=T)

  return (annotatedUniprot)
}


