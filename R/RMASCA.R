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
#' @param validate Logical. If `TRUE`, give estimates for robustness
#' @param participantColumn String. Name of the column containing participant identification
#' @param minimizeObject Logical. If `TRUE`, remove unnecessary clutter, optimize for validation
#' @param nValFold Partitions when validating
#' @param nValRuns number of validation runs
#' @param validationMethod among  `loo` (leave-one-out, default)
#' @param validationObject Don't worry about me:)
#' @param validationParticipants Don't worry about me:)
#' @return An RMASCA object
#' @export
RMASCA <- function(df,
                   formula,
                   separateTimeAndGroup = TRUE,
                   pAdjustMethod = "BH",
                   participantColumn = FALSE,
                   validate = FALSE,
                   minimizeObject = FALSE,
                   nValFold = 7,
                   nValRuns = 50,
                   validationMethod = "loo",
                   validationObject = NA,
                   validationParticipants = NA){
  if(!is.na(validationObject[1])){
    object <- list(df = validationObject$df, #inherit
                formula = validationObject$formula, #inherit
                separateTimeAndGroup = validationObject$separateTimeAndGroup, #inherit
                pAdjustMethod = validationObject$pAdjustMethod, #inherit
                participantColumn = validationObject$participantColumn, #inherit
                validate = validate, #overwrite
                minimizeObject = minimizeObject, #overwrite
                nValFold = validationObject$nValFold, #inherit
                nValRuns = validationObject$nValRuns, #inherit
                validationMethod = validationObject$validationMethod, #inherit
                validationObject = validationObject, #overwrite
                validationParticipants = validationParticipants #overwrite
    )
  }else{
    object <- list(df = df,
                   formula = formula,
                   separateTimeAndGroup = separateTimeAndGroup,
                   pAdjustMethod = pAdjustMethod,
                   participantColumn = participantColumn,
                   validate = validate,
                   minimizeObject = minimizeObject,
                   nValFold = nValFold,
                   nValRuns = nValRuns,
                   validationMethod = validationMethod,
                   validationObject = validationObject,
                   validationParticipants = validationParticipants
    )
    class(object) <- "RMASCA"
  }

  object <- sanitizeObject(object)
  object <- getLMECoefficients(object)
  object <- removeCovars(object)
  object <- separateLMECoefficients(object)
  object <- getEffectMatrix(object)
  if(object$minimizeObject){
    # To save space, we remove unnecessary embedded data
    object <- removeEmbedded(object)
  }
  object <- doPCA(object)
  object <- cleanPCA(object)
  if(object$validate){
    object <- validate(object)
  }

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

  if(object$minimizeObject){
    object$df <- object$df[object$validationParticipants,]
    partColumn <- which(colnames(object$df) == object$participantColumn)
    object$df[,partColumn] <- factor(object$df[,partColumn])
  }
  object$df$time <- factor(object$df$time)
  object$df$group <- factor(object$df$group)
  object$df$variable <- factor(object$df$variable)

  object$hasGroupTerm <- ifelse(any(formulaTerms == "group"), TRUE, FALSE)
  # if(!object$hasGroupTerm){
  #   object$separateTimeAndGroup <- FALSE
  # }
  object$hasInteractionTerm <- ifelse(any(formulaTerms == "group:time" | formulaTerms == "time:group"), TRUE, FALSE)
  object$covars <- formulaTerms[!(formulaTerms %in% c("time","group","group:time","time:group"))]

  return(object)
}

#' Remove df from objectt
#'
#' This function checks that the input to an RMASCA object is as expected
#'
#' @param object An RMASCA object
#' @return An RMASCA object
removeEmbedded <- function(object){
  object$df <- NULL
  object$validationObject <- NULL
  return(object)
}
