get_res = function(u, cols){
  for(i in 1:10){
    res = httr::GET(u, httr::accept_json())
    js = httr::content(res)$jobStatus
    if(length(js)>0 && js=="RUNNING"){
      Sys.sleep(2)
    } else {
      break
    }
  }

  pagesize = 500
  #total = httr::headers(res)$`x-total-results`
  next_u = paste0(sub("status","uniprotkb/results",u),"?size=500&format=tsv&fields=", cols)
  ret = NULL
  while(TRUE){
    res = httr::GET(next_u, httr::accept_json())
    #print(next_u)
    hdr = httr::headers(res)
    map <- try(data.frame(read.delim(textConnection(httr::content(res)),
                                     header=TRUE, stringsAsFactors=FALSE, quote="")))
    ret <- rbind(ret, map)
    next_u = sub(".*<(.*)>.*","\\1",hdr$link)
    if(length(next_u) == 0)
      break
  }

  do.call("rbind", lapply(split(ret, as.factor(ret$From)), function(f) f[1,]))
}

batchUniprotIDmap = function(IDs, from="UniProtKB_AC-ID", cols){
  args =  list(from=from, to="UniProtKB")
  if(from=="Gene_Name")
    args$taxId = sub(".*\\(([0-9]*)\\).*","\\1", species)
  # Query uniprot server in batches
  id_groups=split(IDs, factor((1:length(IDs)) %/% 5000))
  annotUniprot <- do.call("rbind",lapply(id_groups,FUN=function(q_ids){ try({
    ret <- NULL
    res = httr::POST("https://rest.uniprot.org/idmapping/run",
                     body=c(list(ids=paste(q_ids,collapse=",")), args),
                     httr::accept_json())

    res_id = httr::content(res)$jobId
    print(httr::content(res))
    if(is.null(res_id))
      return(NULL)
    get_res(paste0("https://rest.uniprot.org/idmapping/status/", res_id), paste0(cols,collapse=","))
  })}))

  annotUniprot
}

blank_results <- function(IDs, col_names, noNA=T){
  ret <- data.frame(matrix(ncol=length(col_names), nrow=length(IDs)))
  colnames(ret) = col_names
  ret$OrigID=IDs
  if(noNA)
    ret[is.na(ret)] = ""
  ret
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
  uniprot_columns <- c("accession","cc_function", "cc_subcellular_location", "cc_disease",
                     "go_p", "go_f", "go_c", "go_id")
  uniprot_col_names <- c("ENTRY","Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                       "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID")

  if(genes){
    uniprot_columns = c(uniprot_columns, "gene_names")
    uniprot_col_names= c(uniprot_col_names, "Gene.names")
  }
  
  colnames_annotUniprot <- c("OrigID", uniprot_col_names)

  uniprot_regex="^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
  uniparc_regex="^UPI[0-9A-Z]*$"
  gn_regex="^gene="
  idfilter = lapply(c(uniprot_regex,uniparc_regex,gn_regex),FUN=function(r) grepl(r, IDs))
  idtype = rep("",length(IDs))
  idtype[idfilter[[1]]] = "UniProtKB_AC-ID"
  idtype[idfilter[[2]]] = "UniParc"
  idtype[idfilter[[3]]] = "Gene_Name"

  gene_types = idtype=="Gene_Name"
  IDs[gene_types] = sub(gn_regex, "", IDs[gene_types])

  unknown_types = idtype==""
  if(all(unknown_types))
    return(blank_results(IDs, colnames_annotUniprot))
  search_inds = !unknown_types & !duplicated(IDs)
  idtype = idtype[search_inds]
  IDs_unique = IDs[search_inds]

  idgroups = split(IDs_unique, as.factor(idtype))

  annotUniprot = do.call("rbind",lapply(1:length(idgroups),FUN=function(gi){
    batchUniprotIDmap(idgroups[[gi]],
                      from = names(idgroups[gi]),
                      cols = uniprot_columns)
  }))

  if(is.null(annotUniprot))
    return(blank_results(IDs, colnames_annotUniprot))
  colnames(annotUniprot) = c("OrigID", uniprot_col_names)
  if(genes){
    annotUniprot$Gene.names = gsub(" ",";",annotUniprot$Gene.names)
  }

  annotatedUniprot <- blank_results(IDs, colnames_annotUniprot, noNA=F)
  xi=match(IDs, annotUniprot$OrigID)
  origid = which(colnames(annotatedUniprot) == "OrigID")
  annotatedUniprot = annotUniprot[xi,]
  annotatedUniprot[is.na(annotatedUniprot[,origid]),origid] = IDs[is.na(annotatedUniprot[,origid])]
  annotatedUniprot[is.na(annotatedUniprot)] = ""

  return (annotatedUniprot)
}
