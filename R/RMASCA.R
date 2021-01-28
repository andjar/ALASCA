#' Get an RMASCA object
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param df Data Frame to be analyzed
#' @param formula Regression model
#' @param separateTimeAndGroup Logical: should time and group effect be separated?
#' @param pAdjustMethod Method for correcting p values for multiple testing
#' @return An RMASCA object
#' @export
RMASCA <- function(df, formula, separateTimeAndGroup = TRUE, pAdjustMethod = "BH"){
  object <- list(df = df, formula = formula, separateTimeAndGroup = separateTimeAndGroup, pAdjustMethod = pAdjustMethod)
  class(object) <- "RMASCA"

  object <- sanitizeObject(object)
  object <- getLMECoefficients(object)
  object <- removeCovars(object)
  object <- separateLMECoefficients(object)
  object <- getEffectMatrix(object)
  object <- doPCA(object)
  object <- cleanPCA(object)

  # gt_2 <- getLoadingsPlot(scores = score_time, comp = "PC1", pca = pca_time)
  # gt_4 <- getLoadingsPlot(score_time, comp = "PC2", pca = pca_time)
  # gt_1 <- getScorePlot(PC_time, comp = "PC1", pca = pca_time, separateTimeAndGroup = separateTimeAndGroup)
  # gt_3 <- getScorePlot(PC_time, comp = "PC2", pca = pca_time, separateTimeAndGroup = separateTimeAndGroup)
  # g_time <- ggarrange(gt_1, gt_2, gt_3, gt_4)
  #
  # if(separateTimeAndGroup){
  #   gg_2 <- getLoadingsPlot(score_group, comp = "PC1", pca = pca_group)
  #   gg_4 <- getLoadingsPlot(score_group, comp = "PC2", pca = pca_group)
  #   gg_1 <- getScorePlot(PC_group, comp = "PC1", pca = pca_group, separateTimeAndGroup = separateTimeAndGroup)
  #   gg_3 <- getScorePlot(PC_group, comp = "PC2", pca = pca_group, separateTimeAndGroup = separateTimeAndGroup)
  #   g_group <- ggarrange(gg_1, gg_2, gg_3, gg_4)
  # }
  #
  # if(separateTimeAndGroup){
  #   ggg <- ggarrange(g_time, g_group)
  #   return(ggg)
  # }else{
  #   return(g_time)
  # }
  return(object)
}

#' Sanitize an RMASCA object
#'
#' This function checks that the input to an RMASCA object is as expected
#'
#' @param object An RMASCA object to be sanitized
#' @return An RMASCA object
sanitizeObject <- function(object){
  # Check that the input is as expected
  if(!("time" %in% colnames(object$df))){
    stop("The dataframe must contain a column names 'time'")
  }
  if(!("group" %in% colnames(object$df))){
    stop("The dataframe must contain a column names 'group'")
  }
  if(!("variable" %in% colnames(object$df))){
    stop("The dataframe must contain a column names 'variable'")
  }

  formulaTerms <- colnames(attr(terms.formula(object$formula),"factors"))
  if(!any(grepl("\\|",formulaTerms))){
    stop("The model must contain at least one random effect")
  }

  object$df$time <- factor(object$df$time)
  object$df$group <- factor(object$df$group)
  object$df$variable <- factor(object$df$variable)

  object$hasGroupTerm <- ifelse(any(formulaTerms == "group"), TRUE, FALSE)
  object$hasInteractionTerm <- ifelse(any(formulaTerms == "group:time" | formulaTerms == "time:group"), TRUE, FALSE)
  object$covars <- formulaTerms[!(formulaTerms %in% c("time","group","group:time","time:group"))]

  return(object)
}
