#' Print Omics List
#'
#' Print Omics List
#'
#' S3 method for printing OmicsNotebook state
#' 
#' @param object the OmicsNotebook analysis state
#' 
#' @name ON-print
NULL

#' @rdname ON-print
#' @export 
print.ON = function(object){
  cat("An OmicsNotebook (ON) analysis object.\n")
}

#' @rdname ON-print
#' @export 
summary.ON = function(object){
  list(NumGroups=length(unique(object$annot$Group)),
       NumSamples=length(unique(object$annot$SampleName)),
       NumDataFiles=which(colnames(object$annot)=="SampleName.1")-3, # not exactly safe, but oh well
       Calls=object$calls)
}
