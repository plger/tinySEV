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


.homogenizeDEA <- function(x, keepTop=TRUE){
  x <- as.data.frame(x)
  colnames(x) <- gsub("log2FoldChange|log2Fold|log2FC|log2\\(fold_change\\)|log2\\.fold_change\\.",
                      "logFC", colnames(x))
  
  abf <- colnames(df)[which(colnames(df) %in% c("meanExpr", "AveExpr", 
                                                "baseMean", "logCPM"))]
  if (length(abf) == 1) {
    x$meanExpr <- df[, abf]
    if (abf == "baseMean") 
      x$meanExpr <- log(x$meanExpr + 1)
  }else if(all(c("value_1","value_2") %in% colnames(x))){ # cufflinks
    x$meanExpr <- log(1+x$value_1+x$value_2)
  }
  colnames(x) <- gsub("P\\.Value|pvalue|p_value|pval", "PValue", colnames(x))
  colnames(x) <- gsub("padj|adj\\.P\\.Val|q_value|qval", "FDR", colnames(x))
  if (!("FDR" %in% colnames(x))) 
    x$FDR <- p.adjust(x$PValue, method = "fdr")
  f <- grep("^logFC$",colnames(x),value=TRUE)
  if(length(f)==0) f <- grep("logFC",colnames(x),value=TRUE)
  if(length(f)==0) warning("No logFC found.")
  if(length(f)>1){
    message("Using ",f[1])
    x[["logFC"]] <- x[[f[1]]]
  }
  x$FDR[is.na(x$FDR)] <- 1
  x <- x[!is.na(x$logFC),]
  if(!is.null(x$PValue)){
    x <- x[!is.na(x$PValue),]
    x[order(x$PValue),]
  }else{
    x <- x[order(x$FDR),]
  }
  x
}