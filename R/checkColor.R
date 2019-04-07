#-------------------------------------------------
#' Check color
#'
#' Internal function to check that colors are valid
#' 
#' @param x vector of potential hex colors
#'
#' @return TRUE or FALSE
#' 
#' @examples
#'  
#' 
#' @export
checkColor <- function(x){
  #return(all(x%in%colors() | grepl("^#(\\d|[a-f]){6,8}$",x,ignore.case=TRUE)))
  return(all(grepl("^#(\\d|[a-f]){6,8}$",x,ignore.case=TRUE)))
}

