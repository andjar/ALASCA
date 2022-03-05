#' Organize the ALASCA model construction
#'
#' This function builds the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
buildModel <- function(object) {
  if (!object$minimizeObject) {
    # This is not a validation run
    cat("Calculating ", object$method, " coefficients...\n")
  }
  
  if (object$reduceDimensions) {
    if (object$doDebug) currentTs <- Sys.time()
    object <- reduce_dimensions(object)
    if (object$doDebug) cat("* reduce_dimensions:", Sys.time() - currentTs, "s\n")
  }

  if (object$doDebug) currentTs <- Sys.time()
  object <- runRegression(object)
  if (object$doDebug) cat("* runRegression:", Sys.time() - currentTs, "s\n")
  if (!object$useRfast) {
    # With Rfast, we've already got the coefficients
    object <- getRegressionCoefficients(object)
  }

  if (!object$minimizeObject) {
    # This is not a validation run
    cat("Finished calculating regression coefficients!\n")
  }
  
  if (object$doDebug) currentTs <- Sys.time()
  object <- removeCovars(object)
  if (object$doDebug) cat("* removeCovars:", Sys.time() - currentTs, "s\n")
  if (object$doDebug) currentTs <- Sys.time()
  object <- separateLMECoefficients(object)
  if (object$doDebug) cat("* separateLMECoefficients:", Sys.time() - currentTs, "s\n")
  if (object$doDebug) currentTs <- Sys.time()
  object <- getEffectMatrix(object)
  if (object$doDebug) cat("* getEffectMatrix:", Sys.time() - currentTs, "s\n")
  if (object$doDebug) currentTs <- Sys.time()
  object <- doPCA(object)
  if (object$doDebug) cat("* doPCA:", Sys.time() - currentTs, "s\n")
  if (object$doDebug) currentTs <- Sys.time()
  object <- cleanPCA(object)
  if (object$doDebug) cat("* cleanPCA:", Sys.time() - currentTs, "s\n")
  if (object$doDebug) currentTs <- Sys.time()
  object <- cleanALASCA(object)
  if (object$doDebug) cat("* cleanALASCA:", Sys.time() - currentTs, "s\n")
  if (object$method %in% c("LM", "LMM")) {
    if (object$minimizeObject) {
      if (object$validateRegression) {
        object <- getRegressionPredictions(object)
      }
    } else {
      if (object$validateRegression) {
        object <- getRegressionPredictions(object)
      }
    }
  }

  return(object)
}

#' Run regressions
#'
#' This function runs the underlying regression models
#'
#' @param object An ALASCA object
#' @return An ALASCA object
runRegression <- function(object) {
  if (object$useRfast & object$method %in% c("LMM")) {
    # start.time <- Sys.time()
    if (any(is.na(object$df[, value]))) {
      stop("Rfast does NOT like NA's! Check your scaling function or value column.")
    }
    object$RegressionCoefficients <- rbindlist(
      lapply(object$variablelist, function(x) {
        df <- object$df[variable == x]
        modmat <- model.matrix(object$newformula, data = df)
        if (object$forceEqualBaseline) {
          modmat <- modmat[, !grepl(paste0("time", levels(object$df$time)[1]), colnames(modmat))]
        }
        data.frame(
          estimate = Rfast::rint.reg(
            y = df[, value],
            x = modmat[, 2:ncol(modmat)],
            id = as.numeric(factor(df[, ID])),
            ranef = FALSE
          )$be,
          pvalue = NA,
          covar = as.character(x),
          variable = colnames(modmat)
        )
      })
    )
    # end.time <- Sys.time()
    # cat("\n\n",end.time - start.time,"\n")
    return(object)
  } else if (!object$useRfast & object$method %in% c("LM")) {
    object$regr.model <- lapply(object$variablelist, function(x) {
      modmat <- model.matrix(object$formula, data = object$df[variable == x])
      modmat <- modmat[, -1] # Remove intercept
      if (object$forceEqualBaseline) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0("time", unique(object$df$time)[1]), colnames(modmat))]
      }
      environment(object$newformula) <- environment()
      regr.model <- lm(object$newformula, data = object$df[variable == x])
      attr(regr.model, "name") <- x
      regr.model
    })
  } else if (object$useRfast & object$method %in% c("LM")) {
    # start.time <- Sys.time()
    if (any(is.na(object$df[, value]))) {
      stop("Rfast does NOT like NA's! Check your scaling function or value column.")
    }
    object$RegressionCoefficients <- rbindlist(
      lapply(object$variablelist, function(x) {
        df <- object$df[variable == x]
        modmat <- model.matrix(object$formula, data = object$df[variable == x])
        if (object$forceEqualBaseline) {
          # Remove interaction between group and first time point
          modmat <- modmat[, !grepl(paste0("time", unique(object$df$time)[1]), colnames(modmat))]
        }
        data.frame(
          estimate = Rfast::lmfit(
            y = df[, value],
            x = modmat
          )$be,
          pvalue = NA,
          covar = as.character(x),
          variable = colnames(modmat)
        )
      })
    )
    # end.time <- Sys.time()
    # cat("\n\n",end.time - start.time,"\n")
    return(object)
  } else if (object$method %in% c("LMM")) {
    object$regr.model <- lapply(object$variablelist, function(x) {
      modmat <- model.matrix(object$formula, data = object$df[variable == x])
      modmat <- modmat[, -1] # Remove intercept
      if (object$forceEqualBaseline) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0("time", levels(object$df$time)[1]), colnames(modmat), fixed = TRUE)]
      }
      # modmat <- modmat[,ncol(modmat):1]
      environment(object$newformula) <- environment()
      regr.model <- lmerTest::lmer(object$newformula, data = object$df[variable == x])
      attr(regr.model, "name") <- x
      regr.model
    })
  } else if (object$method == "KM") {
    object$regr.model <- lapply(object$variablelist, function(x) {
      regr.model <- survival::coxph(
        formula(paste("survival::Surv(value, ifelse(belowLowerLimit, 0, 1)) ~", as.character(object$formula)[3])),
        data = subset(object$df, variable == x)
      )
      attr(regr.model, "name") <- x
      regr.model
    })
  } else if (object$method == "KMM") {
    object$regr.model <- lapply(object$variablelist, function(x) {
      regr.model <- coxme::coxme(
        formula(paste("survival::Surv(value, ifelse(belowLowerLimit, 0, 1)) ~", as.character(object$formula)[3], "+ (1|ID)")),
        data = subset(object$df, variable == x)
      )
      attr(regr.model, "name") <- x
      regr.model
    })
  }
  names(object$regr.model) <- object$variablelist
  return(object)
}

#' Get regression coefficients
#'
#' This function extract the regression coefficients for the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
getRegressionCoefficients <- function(object) {

  fdf <- data.table::rbindlist(
    lapply(object$regr.model, function(y) {
      if (object$method %in% c("LM", "Lim")) {
        tmp_ef <- coef(y)
        a <- as.data.frame(summary(y)[["coefficients"]][, c(1, 4)])
      } else {
        tmp_ef <- lme4::fixef(y)
        a <- as.data.frame(summary(y)[["coefficients"]][, c(1, 5)])
      }
      a$covar <- as.character(attr(y, "name"))
      a$variable <- rownames(a)
      rownames(a) <- NULL
      a
    })
  )
  fdf$variable <- gsub("modmat", "", fdf$variable, fixed = TRUE)
  colnames(fdf) <- c("estimate", "pvalue", "covar", "variable")
  object$RegressionCoefficients <- fdf
  
  if (!is.na(object$pAdjustMethod)) {
    cat("Adjusting p values...\n")
    object$RegressionCoefficients[, pvalue_unadj := pvalue]
    object$RegressionCoefficients[, pvalue := p.adjust(pvalue, method = object$pAdjustMethod), by = variable]
  }
  
  return(object)
}

#' Remove unwanted covariables
#'
#' This function removes coefficients that we do not want in our PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
removeCovars <- function(object) {
  object$CovarCoefficients <- data.frame()
  for (i in unique(object$covars)) {
    object$CovarCoefficients <- rbind(object$CovarCoefficients, subset(object$RegressionCoefficients, substr(variable, 1, nchar(i)) == i))
    object$RegressionCoefficients <- subset(object$RegressionCoefficients, substr(variable, 1, nchar(i)) != i)
  }
  if (object$reduceDimensions) {
    object$CovarCoefficients <- rbindlist(lapply(unique(object$CovarCoefficients$variable), function(v){
      ref <- object$CovarCoefficients[variable == v]
      data.frame(
        variable = v,
        covar = rownames(object$Limm$loadings),
        pvalue = NA,
        estimate = rowSums(object$Limm$loadings*ref$estimate[match(colnames(object$Limm$loadings), ref$covar)][col(object$Limm$loadings)])
      )
    }))
  }
  return(object)
}

#' Separate time and group effects
#'
#' This function separates time and group variables if separateTimeAndGroup = TRUE
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
separateLMECoefficients <- function(object) {
  object$RegressionCoefficients$comp <- "TIME"
  if (object$separateTimeAndGroup) {
    object$RegressionCoefficients$comp[!(object$RegressionCoefficients$variable == "(Intercept)" |
      (substr(object$RegressionCoefficients$variable, 1, 4) == "time" & !grepl(":", object$RegressionCoefficients$variable, fixed = "TRUE")))] <- "GROUP"
  }
  return(object)
}

#' Get effect matrix
#'
#' This function calculates the effect matrix
#'
#' @param object An ALASCA object
#' @return An ALASCA object
getEffectMatrix <- function(object) {
  if (!object$minimizeObject) {
    cat("Calculating effect matrix\n")
  }
  parts <- object$df[variable == object$variablelist[1]]
  # parts <- object$df[!duplicated(cbind(object$df$ID, object$df$time))]
  Dmatrix <- model.matrix(object$formula, data = object$df[variable == object$variablelist[1]])
  # Dmatrix <- Dmatrix[,ncol(Dmatrix):1]

  if (object$separateTimeAndGroup) {
    BmatrixTime <- object$RegressionCoefficients[object$RegressionCoefficients$comp == "TIME", c("covar", "estimate", "variable")]
    BmatrixTime <- reshape2::dcast(BmatrixTime, formula = variable ~ covar, value.var = "estimate")
    selectDColumnsTime <- colnames(Dmatrix) %in% BmatrixTime$variable
    rowOrder <- c()
    for (i in seq_len(nrow(BmatrixTime))) {
      rowOrder[i] <- which(BmatrixTime$variable == colnames(Dmatrix[, selectDColumnsTime])[i])
    }
    BmatrixTime <- BmatrixTime[rowOrder, ]
    if (any(colnames(Dmatrix[, selectDColumnsTime]) != BmatrixTime$variable)) {
      stop("Column mismatch for time in getEffectMatrix")
    }
    AmatrixTime <- as.data.frame(as.matrix(Dmatrix[, selectDColumnsTime]) %*% as.matrix(BmatrixTime[, -1]))
    AmatrixTime$comp <- "TIME"

    BmatrixGroup <- object$RegressionCoefficients[object$RegressionCoefficients$comp == "GROUP", c("covar", "estimate", "variable")]
    BmatrixGroup <- reshape2::dcast(BmatrixGroup, formula = variable ~ covar, value.var = "estimate")
    selectDColumnsGroup <- colnames(Dmatrix) %in% BmatrixGroup$variable
    rowOrder <- c()
    for (i in seq_len(nrow(BmatrixGroup))) {
      rowOrder[i] <- which(BmatrixGroup$variable == colnames(Dmatrix[, selectDColumnsGroup])[i])
    }
    BmatrixGroup <- BmatrixGroup[rowOrder, ]
    if (any(colnames(Dmatrix[, selectDColumnsGroup]) != BmatrixGroup$variable)) {
      stop("Column mismatch for group in getEffectMatrix")
    }
    AmatrixGroup <- as.data.frame(as.matrix(Dmatrix[, selectDColumnsGroup]) %*% as.matrix(BmatrixGroup[, -1]))
    AmatrixGroup$comp <- "GROUP"
    object$effect.matrix <- rbind(AmatrixTime, AmatrixGroup)
  } else {
    BmatrixTime <- object$RegressionCoefficients[object$RegressionCoefficients$comp == "TIME", c("covar", "estimate", "variable")]
    BmatrixTime <- reshape2::dcast(BmatrixTime, formula = variable ~ covar, value.var = "estimate")

    selectDColumnsTime <- colnames(Dmatrix) %in% BmatrixTime$variable
    rowOrder <- c()
    for (i in seq_len(nrow(BmatrixTime))) {
      rowOrder[i] <- which(BmatrixTime$variable == colnames(Dmatrix[, selectDColumnsTime])[i])
    }
    BmatrixTime <- BmatrixTime[rowOrder, ]
    if (any(colnames(Dmatrix[, selectDColumnsTime]) != BmatrixTime$variable)) {
      stop("Column mismatch for time in getEffectMatrix")
    }
    AmatrixTime <- as.data.frame(as.matrix(Dmatrix[, selectDColumnsTime]) %*% as.matrix(BmatrixTime[, -1]))
    AmatrixTime$comp <- "TIME"
    object$effect.matrix <- AmatrixTime
  }

  object$parts$time <- parts$time
  if (object$keepTerms != "") {
    keepTerms <- c("group", object$keepTerms)
    object$parts$group <- apply(parts[, ..keepTerms], 1, paste, collapse = " - ")
  } else {
    object$parts$group <- parts$group
  }
  if (!object$minimizeObject) {
    cat("Finished calculating effect matrix!\n")
  }
  return(object)
}
