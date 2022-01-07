#' Get residuals
#'
#' This function returns the residual of the underlying linear mixed models
#'
#' @param object An RMASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @return A list of residuals per variable
#' 
#' @examples
#' load("PE.Rdata")
#' residuals(model, variable = c("IL-6", "PlGF"))
#' 
#' @export
residuals.ALASCA <- function(object, variable = NA){
  if(any(is.na(variable))){
    return(lapply(object$regr.model, residuals))
  }else{
   varList <- names(object$regr.model)
   resList <- lapply(seq_along(object$regr.model), function(x) if(varList[x] %in% variable){residuals(object$regr.model[[x]])})
   names(resList) <- names(object$regr.model)
   resList[sapply(resList, is.null)] <- NULL
   return(resList)
  }
}

#' Plot residuals
#'
#' This function plots the residuals of the underlying linear mixed models
#'
#' @param object An RMASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @param plottitle If `TRUE` (default), include variable name as title
#' @return A list of ggplot2 objects per variable
#' 
#' @export
plotresiduals <- function(object, variable = NA, plottitle = TRUE, myTheme = ggplot2::theme_classic()){
   if(any(is.na(variable))){
      resList <- lapply(object$regr.model, residuals)
   }else{
      varList <- names(object$regr.model)
      resList <- lapply(seq_along(object$regr.model), function(x) if(varList[x] %in% variable){residuals(object$regr.model[[x]])})
      names(resList) <- names(object$regr.model)
      resList[sapply(resList, is.null)] <- NULL
   }
   lapply(seq_along(resList), function(x){
      g <- ggplot2::ggplot(data.frame(Residuals = resList[[x]]),
         ggplot2::aes(sample = Residuals)) + 
         ggplot2::stat_qq() + ggplot2::stat_qq_line()
      if(plottitle) g <- g + ggplot2::labs(title = names(resList)[x])
      g + myTheme
   })
}