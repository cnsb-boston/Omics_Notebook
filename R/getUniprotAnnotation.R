batchUniprotQuery = function(IDs, from="ACC", cols, col_names){
  id_groups=split(IDs, factor((1:length(IDs)) %/% 5000))
  
  # Query uniprot server in batches
  annotUniprot <- do.call("rbind",lapply(id_groups,FUN=function(q_ids){ try({
    ret <- NULL
    res = httr::POST("https://www.uniprot.org/uploadlists/",
                     body=list(query=paste(q_ids,collapse=","), format="tab", from=from, to="ACC",
                               columns=paste0(cols,collapse=",")))
    info_annot <- try(data.frame(read.delim(textConnection(httr::content(res)),
                                            header=TRUE, stringsAsFactors=FALSE, quote="")))
    if (class(info_annot) != 'try-error' &&
        ncol(info_annot)==length(cols)+1 &&
        colnames(info_annot)[2]=="Function..CC."){
      #info_annot=info_annot[,-ncol(info_annot)] # yourlist...
      colnames(info_annot) <- c(col_names,"OrigID")
      ret <- info_annot
    }
    Sys.sleep(2)
    ret
  }) 
  }))

  annotUniprot
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
getUniprotAnnotation <- function(IDs, genes=F){
  
  # Uniprot entries to fetch (and col names)
  uniprot_columns <- c("id","comment(FUNCTION)", "comment(SUBCELLULAR LOCATION)", "comment(DISEASE)",
                     "go(biological process)", "go(molecular function)", "go(cellular component)", "go-id", "database(Reactome)", "database(KEGG)",
                     "database(BioCyc)", "database(Ensembl)", "database(ChEMBL)", "database(IntAct)","database(STRING)")
  uniprot_col_names <- c("ENTRY","Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                       "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID",
                       "BioCyc_ID", "Ensembl_ID", "ChEMBL_ID", "IntAct_ID","STRING_ID")

  if(genes){
    uniprot_columns = c(uniprot_columns, "genes")
    uniprot_col_names= c(uniprot_col_names, "Gene.names")
  }
  
  colnames_annotUniprot <- c("ENTRY", uniprot_col_names)

  # List of uniprot IDs (remove duplicates for speed)
  IDs_unique <- IDs[!duplicated(IDs)] #get unique IDs

  uniprot_regex="^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
  uniparc_regex="^UPI[0-9A-Z]*$"
  idfilter = lapply(c(uniprot_regex,uniparc_regex),FUN=function(r) grepl(r, IDs_unique))
  idtype = rep("",length(IDs_unique))
  idtype[idfilter[[1]]] = "ACC"
  idtype[idfilter[[2]]] = "UPARC"

  unknown_types = idtype==""
  idtype = idtype[!unknown_types]
  IDs_unique = IDs_unique[!unknown_types]

  idgroups = split(IDs_unique, as.factor(idtype))

  annotUniprot = do.call("rbind",lapply(1:length(idgroups),FUN=function(gi){
    batchUniprotQuery(idgroups[[gi]],
               from = names(idgroups[gi]),
               cols = uniprot_columns,
               col_names = uniprot_col_names)
  }))

  if(genes){
    annotUniprot$Gene.names = gsub(" ",";",annotUniprot$Gene.names)
  }

  annotatedUniprot <- data.frame(matrix(ncol=length(uniprot_col_names), nrow=length(IDs)))
  colnames(annotatedUniprot) <- uniprot_col_names 
  annotatedUniprot[,"OrigID"] <- IDs
  xi=match(IDs, annotUniprot$OrigID)
  annotatedUniprot[,-ncol(annotatedUniprot)] = annotUniprot[xi,-ncol(annotUniprot)]
  annotatedUniprot[is.na(annotatedUniprot[,1]),1] = IDs[is.na(annotatedUniprot[,1])]
  annotatedUniprot[is.na(annotatedUniprot)] = ""

  return (annotatedUniprot)
}


