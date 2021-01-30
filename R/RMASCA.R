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
#' @param pAdjustMethod Method for correcting p values for multiple testing, see p.adjust.methods
#' @param validate Logical. If `TRUE`, give estimates for robustness
#' @param participantColumn String. Name of the column containing participant identification
#' @param minimizeObject Logical. If `TRUE`, remove unnecessary clutter, optimize for validation
#' @param scale If `TRUE` (default), each variable is scaled to unit SD and zero mean
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
                   pAdjustMethod = NA,
                   participantColumn = FALSE,
                   validate = FALSE,
                   scale = TRUE,
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
                scale = validationObject$scale, #inherit
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
                   scale = scale,
                   minimizeObject = minimizeObject,
                   nValFold = nValFold,
                   nValRuns = nValRuns,
                   validationMethod = validationMethod,
                   validationObject = validationObject,
                   validationParticipants = validationParticipants
    )
    class(object) <- "RMASCA"
  }

  start.time <- Sys.time()
  object <- sanitizeObject(object)
  end.time <- Sys.time()
  cat("Time 1: ", end.time - start.time, "\n")
  start.time <- Sys.time()
  object <- getLMECoefficients(object)
  end.time <- Sys.time()
  cat("Time 2: ", end.time - start.time, "\n")
  start.time <- Sys.time()
  object <- removeCovars(object)
  end.time <- Sys.time()
  cat("Time 3: ", end.time - start.time, "\n")
  start.time <- Sys.time()
  object <- separateLMECoefficients(object)
  end.time <- Sys.time()
  cat("Time 4: ", end.time - start.time, "\n")
  start.time <- Sys.time()
  object <- getEffectMatrix(object)
  end.time <- Sys.time()
  cat("Time 5: ", end.time - start.time, "\n")
  start.time <- Sys.time()
  object <- doPCA(object)
  end.time <- Sys.time()
  cat("Time 6: ", end.time - start.time, "\n")
  start.time <- Sys.time()
  object <- cleanPCA(object)
  end.time <- Sys.time()
  cat("Time 7: ", end.time - start.time, "\n")
  if(object$minimizeObject){
    # To save space, we remove unnecessary embedded data
    object <- removeEmbedded(object)
  }
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

  if(any(is.na(object$df$time) | is.na(object$df$group))){
    cat("\n\n!!! -> Oh dear, at least on of your rows is missing either time or group. I have removed it/them for now, but you should check if this is a problem...\n\n")
    object$df <- object$df[!is.na(object$df$time) & !is.na(object$df$group),]
  }
  if(!object$minimizeObject & object$scale){
    cat("Scaling data...\n")
    for(i in unique(object$df$variable)){
      valColumn <- which(as.character(object$formula)[2] == colnames(object$df))
      object$df[object$df$variable == i,valColumn] <- (object$df[object$df$variable == i,valColumn] - mean(object$df[object$df$variable == i,valColumn]))/sd(object$df[object$df$variable == i,valColumn])
    }
  }

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
  object$lmer.models <- NULL
  object$LMM.coefficients <- NULL
  object$effect.matrix <- NULL
  return(object)
}
