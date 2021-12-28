.modify_stop_propagation <- function(x){
  if(!is.null(x$children[[1]]))
    x$children[[1]]$attribs$onclick = "event.stopPropagation()"
  x
}


#' dround
#'
#' Trim to a certain number of digits (equivalent to `format(...,digits=digits)`, except that the output is numeric)-
#'
#' @param x A vector of numeric values
#' @param digits The number of digits to keep
#' @param roundGreaterThan1 Whether to trim also numbers greater than 1 (default TRUE)
#'
#' @return A numeric vector of the same length as `x`
#' @export
#'
#' @examples
#' dround( c(0.00002345, 554356, 12.56) )
dround <- function(x, digits=3, roundGreaterThan1=TRUE){
  if(is.matrix(x) || is.data.frame(x)){
    for(i in 1:ncol(x)){
      if(is.numeric(x[,i])){
        tryCatch(x[,i] <- dround(x[,i], digits, roundGreaterThan1), error=function(e) warning(e))
      }
    }
    return(x)
  }
  if(roundGreaterThan1){
    w <- 1:length(x)
  }else{
    w <- which(abs(x)<1)
  }
  if(length(w)==0) return(x)
  e <- ceiling(-log10(abs(x[w])))
  x[w] <- round(10^e*x[w],digits-1)/10^e
  x
}


.getWordsFromString <- function(ss){
  ss <- lapply(strsplit(gsub(" ", ",", ss, fixed=TRUE), ",", fixed=TRUE),
               FUN=function(x) x[which(x != "")])
  if(length(ss)==1) return(ss[[1]])
  ss
}

#' grepGene
#'
#' Finds genes names in a vector or Summarized Experiment object
#'
#' @param x A vector of name composites (e.g. ensemblID.symbol), or a
#' SummarizedExperiment object with such row names
#' @param g A set of genes names to select
#' @param ignore.case Logical; whether to ignore case when matchiing
#'
#' @return The subsetted elements of `x` matching `g`
#'
#' @export
#' @examples
#' x <- c("ENSG00000000003.TSPAN6", "ENSG00000000005.TNMD",
#'        "ENSG00000000419.DPM1")
#' grepGene(x, c("TSPAN6", "TNMD"))
grepGene <- function(x, g, ignore.case=TRUE){
  if(!is.character(x)){
    g <- grepGene(row.names(x), g)
    return(x[g,drop=FALSE])
  }
  if(all(g %in% x)) return(g)
  g <- paste0("^",g,"\\.|^",g,"$|\\.",g,"$", collapse="|")
  g <- grep(g, x, value=TRUE, ignore.case=ignore.case)
  return(unique(g))
}
