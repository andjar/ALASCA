#' Organize the ALASCA model construction
#'
#' This function builds the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
buildModel <- function(object) {
  if (!object$minimize_object) {
    # This is not a validation run
    object <- add_to_log(object, message = paste0("Calculating ", object$method, " coefficients"), print = TRUE)
  }
  
  if (object$reduce_dimensions) {
    if (object$do_debug) currentTs <- Sys.time()
    object <- reduce_dimensions(object)
    if (object$do_debug) cat("* reduce_dimensions:", Sys.time() - currentTs, "s\n")
  }

  if (object$do_debug) currentTs <- Sys.time()
  object <- run_regression(object)
  if (object$do_debug) cat("* runRegression:", Sys.time() - currentTs, "s\n")
  if (!object$use_Rfast) {
    # With Rfast, we've already got the coefficients
    object <- get_regression_coefficients(object)
  }

  if (!object$minimize_object) {
    # This is not a validation run
    object <- add_to_log(object, message = "Finished calculating regression coefficients!", print = TRUE)
  }
  
  if (object$do_debug) currentTs <- Sys.time()
  object <- remove_covars(object)
  if (object$do_debug) cat("* removeCovars:", Sys.time() - currentTs, "s\n")
  if (object$do_debug) currentTs <- Sys.time()
  object <- separate_regression_coefficients(object)
  if (object$do_debug) cat("* separateLMECoefficients:", Sys.time() - currentTs, "s\n")
  if (object$do_debug) currentTs <- Sys.time()
  object <- get_effect_matrix(object)
  if (object$do_debug) cat("* getEffectMatrix:", Sys.time() - currentTs, "s\n")
  if (object$do_debug) currentTs <- Sys.time()
  object <- do_pca(object)
  if (object$do_debug) cat("* doPCA:", Sys.time() - currentTs, "s\n")
  if (object$do_debug) currentTs <- Sys.time()
  object <- clean_pca(object)
  if (object$do_debug) cat("* clean_pca:", Sys.time() - currentTs, "s\n")
  if (object$do_debug) currentTs <- Sys.time()
  object <- clean_alasca(object)
  if (object$do_debug) cat("* clean_alasca:", Sys.time() - currentTs, "s\n")
  if (object$method %in% c("LM", "LMM")) {
    if (object$minimize_object) {
      if (object$validate_regression) {
        object <- get_regression_predictions(object)
      }
    } else {
      if (object$validate_regression) {
        object <- get_regression_predictions(object)
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
run_regression <- function(object) {
  df_by_variable <- split(object$df, object$df$variable)
  if (object$use_Rfast && object$method %in% c("LMM")) {
    # start.time <- Sys.time()
    if (any(is.na(object$df[, value]))) {
      add_to_log(object, message = "Rfast does NOT like NA's! Check your scaling function or value column.", level = "STOP")
    }
    object$regression_coefficients <- rbindlist(
      lapply(object$variablelist, function(x) {
        modmat <- model.matrix(object$new_formula, data = df_by_variable[[x]])
        if (object$equal_baseline) {
          modmat <- modmat[, !grepl(paste0("time", object$timelist[1]), colnames(modmat))]
        }
        data.frame(
          estimate = Rfast::rint.reg(
            y = df_by_variable[[x]][, value],
            x = modmat[, 2:ncol(modmat)],
            id = as.numeric(factor(df_by_variable[[x]][, ID])),
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
  } else if (!object$use_Rfast & object$method %in% c("LM")) {
    object$regression_model <- lapply(object$variablelist, function(x) {
      modmat <- model.matrix(object$formula, data = df_by_variable[[x]])
      modmat <- modmat[, -1] # Remove intercept
      if (object$equal_baseline) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0("time", object$timelist[1]), colnames(modmat))]
      }
      environment(object$new_formula) <- environment()
      regression_model <- lm(object$new_formula, data = df_by_variable[[x]])
      attr(regression_model, "name") <- x
      regression_model
    })
  } else if (object$use_Rfast & object$method %in% c("LM")) {
    # start.time <- Sys.time()
    if (any(is.na(object$df[, value]))) {
      add_to_log(object, message = "Rfast does NOT like NA's! Check your scaling function or value column.", level = "STOP")
    }
    object$regression_coefficients <- rbindlist(
      lapply(object$variablelist, function(x) {
        modmat <- model.matrix(object$formula, data = df_by_variable[[x]])
        if (object$equal_baseline) {
          # Remove interaction between group and first time point
          modmat <- modmat[, !grepl(paste0("time", object$timelist[1]), colnames(modmat))]
        }
        data.frame(
          estimate = Rfast::lmfit(
            y = df_by_variable[[x]][, value],
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
    object$regression_model <- lapply(object$variablelist, function(x) {
      modmat <- model.matrix(object$formula, data = df_by_variable[[x]])
      modmat <- modmat[, -1] # Remove intercept
      if (object$equal_baseline) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0("time", object$timelist[1]), colnames(modmat), fixed = TRUE)]
      }
      # modmat <- modmat[,ncol(modmat):1]
      environment(object$new_formula) <- environment()
      regression_model <- lmerTest::lmer(object$new_formula, data = df_by_variable[[x]])
      attr(regression_model, "name") <- x
      regression_model
    })
  }
  names(object$regression_model) <- object$variablelist
  return(object)
}

#' Get regression coefficients
#'
#' This function extract the regression coefficients for the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_regression_coefficients <- function(object) {

  fdf <- data.table::rbindlist(
    lapply(object$regression_model, function(y) {
      if (object$method %in% c("LM")) {
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
  object$regression_coefficients <- fdf
  
  if (!is.na(object$p_adjust_method)) {
    object <- add_to_log(object, message = "Adjusting p values", print = TRUE)
    object$regression_coefficients[, pvalue_unadj := pvalue]
    object$regression_coefficients[, pvalue := p.adjust(pvalue, method = object$p_adjust_method), by = variable]
  }
  
  return(object)
}

#' Remove unwanted covariables
#'
#' This function removes coefficients that we do not want in our PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
remove_covars <- function(object) {
  object$covar_coefficients <- data.frame()
  for (i in unique(object$covars)) {
    object$covar_coefficients <- rbind(object$covar_coefficients, subset(object$regression_coefficients, substr(variable, 1, nchar(i)) == i))
    object$regression_coefficients <- subset(object$regression_coefficients, substr(variable, 1, nchar(i)) != i)
  }
  if (object$reduce_dimensions) {
    object$covar_coefficients <- rbindlist(lapply(unique(object$covar_coefficients$variable), function(v){
      ref <- object$covar_coefficients[variable == v]
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
#' This function separates time and group variables if separate_time_and_group = TRUE
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
separate_regression_coefficients <- function(object) {
  object$regression_coefficients$comp <- "TIME"
  if (object$separate_time_and_group) {
    object$regression_coefficients$comp[!(object$regression_coefficients$variable == "(Intercept)" |
      (substr(object$regression_coefficients$variable, 1, 4) == "time" & !grepl(":", object$regression_coefficients$variable, fixed = "TRUE")))] <- "GROUP"
  }
  return(object)
}

#' Get effect matrix
#'
#' This function calculates the effect matrix
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_effect_matrix <- function(object) {
  if (!object$minimize_object) {
    object <- add_to_log(object, message = "Calculating effect matrix", print = TRUE)
  }
  parts <- object$df[variable == object$variablelist[1]]
  # parts <- object$df[!duplicated(cbind(object$df$ID, object$df$time))]
  Dmatrix <- model.matrix(object$formula, data = object$df[variable == object$variablelist[1]])
  # Dmatrix <- Dmatrix[,ncol(Dmatrix):1]

  if (object$separate_time_and_group) {
    BmatrixTime <- object$regression_coefficients[object$regression_coefficients$comp == "TIME", c("covar", "estimate", "variable")]
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

    BmatrixGroup <- object$regression_coefficients[object$regression_coefficients$comp == "GROUP", c("covar", "estimate", "variable")]
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
    object$effect_matrix <- rbind(AmatrixTime, AmatrixGroup)
  } else {
    BmatrixTime <- object$regression_coefficients[object$regression_coefficients$comp == "TIME", c("covar", "estimate", "variable")]
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
    object$effect_matrix <- AmatrixTime
  }

  object$parts$time <- parts$time
  if (object$keep_terms != "") {
    keep_terms <- c("group", object$keep_terms)
    object$parts$group <- apply(parts[, ..keep_terms], 1, paste, collapse = " - ")
  } else {
    object$parts$group <- parts$group
  }
  if (!object$minimize_object) {
    object <- add_to_log(object, message = "Finished calculating effect matrix!", print = TRUE)
  }
  return(object)
}
