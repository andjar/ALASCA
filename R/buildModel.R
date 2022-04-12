#' Organize the ALASCA model construction
#'
#' This function builds the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
build_model <- function() {
  if (!self$minimize_self) {
    # This is not a validation run
    log4r::debug(self$log, "Starting to build model")
    log4r::info(self$log, paste0("Calculating ", self$method, " coefficients"))
  }
  
  if (self$reduce_dimensions) {
    log4r::debug(self$log, "Starting to reduce dimensions")
    self$reduce_dimensions()
    log4r::debug(self$log, "Finished to reduce dimensions")
  }

  if (self$do_debug) currentTs <- Sys.time()
  self$run_regression()
  log4r::debug(self$log, "Starting to calculate regression coefficients")
  if (!self$use_Rfast) {
    # With Rfast, we've already got the coefficients
    self$get_regression_coefficients()
  }

  if (!self$minimize_self) {
    # This is not a validation run
    log4r::info(self$log, "Finished calculating regression coefficients!")
  }
  
  if (self$do_debug) currentTs <- Sys.time()
  self$remove_covars()
  if (self$do_debug) cat("* removeCovars:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$separate_regression_coefficients()
  if (self$do_debug) cat("* separateLMECoefficients:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$get_effect_matrix()
  if (self$do_debug) cat("* getEffectMatrix:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$do_pca()
  if (self$do_debug) cat("* doPCA:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$clean_pca()
  if (self$do_debug) cat("* clean_pca:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$clean_alasca()
  if (self$do_debug) cat("* clean_alasca:", Sys.time() - currentTs, "s\n")
  if (self$method %in% c("LM", "LMM")) {
    if (self$minimize_self) {
      if (self$validate_regression) {
        self$get_regression_predictions()
      }
    } else {
      if (self$validate_regression) {
        self$get_regression_predictions()
      }
    }
  }

  #invivisble(self)
}

#' Run regressions
#'
#' This function runs the underlying regression models
#'
#' @param object An ALASCA object
#' @return An ALASCA object
run_regression <- function() {

  #df_by_variable <- split(self$df, self$df$variable)
  rows_by_variable <- lapply(self[["variablelist"]], function(x) self[["df"]][, .I[variable == x]] )
  names(rows_by_variable) <- self[["variablelist"]]
  
  if (self$use_Rfast && self$method %in% c("LMM")) {
    # start.time <- Sys.time()
    if (any(is.na(self[["df"]][, value]))) {
      log4r::error("Rfast does NOT like NA's! Check your scaling function or value column.")
      stop()
    }
    if (!self[["minimize_self"]]) {
      self[["modmat"]] <- model.matrix(self[["new_formula"]], data = self$df)
      if (self[["equal_baseline"]]) {
        self[["modmat"]] <- self[["modmat"]][, !grepl(paste0("time", self[["timelist"]][1]), colnames(self[["modmat"]]))]
      }
      self[["cnames_modmat"]] <- colnames(self[["modmat"]])
    }
    self$regression_coefficients <- rbindlist(
      lapply(self[["variablelist"]], function(x) {
        list(
          estimate = Rfast::rint.reg(
            y = self[["df"]][ rows_by_variable[[x]], value],
            x = self[["modmat"]][rows_by_variable[[x]], -1],
            id = as.numeric(factor(self[["df"]][ rows_by_variable[[x]], ID])),
            ranef = FALSE
          )$be,
          pvalue = NA,
          covar = as.character(x),
          variable = self[["cnames_modmat"]]
        )
      })
    )
    # end.time <- Sys.time()
    # cat("\n\n",end.time - start.time,"\n")
    #invisible(self)
  } else if (!self$use_Rfast & self$method %in% c("LM")) {
    self$regression_model <- lapply(self[["variablelist"]], function(x) {
      modmat <- model.matrix(self[["formula"]], data = self$df[ rows_by_variable[[x]] ])
      modmat <- modmat[, -1] # Remove intercept
      if (self[["equal_baseline"]]) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0("time", self[["timelist"]][1]), colnames(modmat))]
      }
      environment(self$new_formula) <- environment()
      regression_model <- lm(self[["new_formula"]], data = self$df[ rows_by_variable[[x]] ])
      attr(regression_model, "name") <- x
      regression_model
    })
  } else if (self$use_Rfast & self$method %in% c("LM")) {
    # start.time <- Sys.time()
    if (any(is.na(self$df[, value]))) {
      log4r::error(self$log, "Rfast does NOT like NA's! Check your scaling function or value column.")
      stop()
    }
    if (!self$minimize_self) {
      self$modmat <- model.matrix(self[["formula"]], data = self$df[ rows_by_variable[[x]] ])
      if (self[["equal_baseline"]]) {
        # Remove interaction between group and first time point
        self[["modmat"]] <- self[["modmat"]][, !grepl(paste0("time", self[["timelist"]][1]), colnames(self[["modmat"]]))]
      }
      self[["cnames_modmat"]] <- colnames(self[["modmat"]])
    }
    self$regression_coefficients <- rbindlist(
      lapply(self[["variablelist"]], function(x) {
        list(
          estimate = Rfast::lmfit(
            y = self[["df"]][ rows_by_variable[[x]], value],
            x = self[["modmat"]][rows_by_variable[[x]], -1]
          )$be,
          pvalue = NA,
          covar = as.character(x),
          variable = self[["cnames_modmat"]]
        )
      })
    )
    # end.time <- Sys.time()
    # cat("\n\n",end.time - start.time,"\n")
    #invisible(self)
  } else if (self$method %in% c("LMM")) {
    self$regression_model <- lapply(self[["variablelist"]], function(x) {
      modmat <- model.matrix(self[["formula"]], data = self$df[ rows_by_variable[[x]] ])
      modmat <- modmat[, -1] # Remove intercept
      if (self[["equal_baseline"]]) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0("time", self[["timelist"]][1]), colnames(modmat), fixed = TRUE)]
      }
      # modmat <- modmat[,ncol(modmat):1]
      environment(self$new_formula) <- environment()
      regression_model <- lmerTest::lmer(self[["new_formula"]], data = self$df[ rows_by_variable[[x]] ])
      attr(regression_model, "name") <- x
      regression_model
    })
  }
  names(self$regression_model) <- self$variablelist
  #invisible(self)
}

#' Get regression coefficients
#'
#' This function extract the regression coefficients for the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_regression_coefficients <- function() {

  fdf <- data.table::rbindlist(
    lapply(self$regression_model, function(y) {
      if (self$method %in% c("LM")) {
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
  self$regression_coefficients <- fdf
  
  if (!is.na(self$p_adjust_method)) {
    log4r::info(self$log, "Adjusting p values")
    self$regression_coefficients[, pvalue_unadj := pvalue]
    self$regression_coefficients[, pvalue := p.adjust(pvalue, method = self$p_adjust_method), by = variable]
  }
  
  #invisible(self)
}

#' Remove unwanted covariables
#'
#' This function removes coefficients that we do not want in our PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
remove_covars <- function() {
  self$covar_coefficients <- data.frame()
  for (i in unique(self$covars)) {
    self$covar_coefficients <- rbind(self$covar_coefficients, subset(self$regression_coefficients, substr(variable, 1, nchar(i)) == i))
    self$regression_coefficients <- subset(self$regression_coefficients, substr(variable, 1, nchar(i)) != i)
  }
  if (self$reduce_dimensions) {
    self$covar_coefficients <- rbindlist(lapply(unique(self$covar_coefficients$variable), function(v){
      ref <- self$covar_coefficients[variable == v]
      list(
        variable = v,
        covar = rownames(self$Limm$loadings),
        pvalue = NA,
        estimate = rowSums(self$Limm$loadings*ref$estimate[match(colnames(self$Limm$loadings), ref$covar)][col(self$Limm$loadings)])
      )
    }))
  }
  #invisible(self)
}

#' Separate time and group effects
#'
#' This function separates time and group variables if separate_time_and_group = TRUE
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
separate_regression_coefficients <- function() {
  self$regression_coefficients$comp <- "TIME"
  if (self$separate_time_and_group) {
    self$regression_coefficients$comp[!(self$regression_coefficients$variable == "(Intercept)" |
      (substr(self$regression_coefficients$variable, 1, 4) == "time" & !grepl(":", self$regression_coefficients$variable, fixed = "TRUE")))] <- "GROUP"
  }
  #invisible(self)
}

#' Get effect matrix
#'
#' This function calculates the effect matrix
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_effect_matrix <- function() {
  if (!self$minimize_self) {
    log4r::info(self$log, "Calculating effect matrix")
  }
  parts <- self$df[variable == self$variablelist[1]]
  # parts <- self$df[!duplicated(cbind(self$df$ID, self$df$time))]
  Dmatrix <- model.matrix(self$formula, data = self$df[variable == self$variablelist[1]])
  # Dmatrix <- Dmatrix[,ncol(Dmatrix):1]

  if (self$separate_time_and_group) {
    BmatrixTime <- self$regression_coefficients[self$regression_coefficients$comp == "TIME", c("covar", "estimate", "variable")]
    BmatrixTime <- dcast(BmatrixTime, formula = variable ~ covar, value.var = "estimate")
    selectDColumnsTime <- colnames(Dmatrix) %in% BmatrixTime$variable
    rowOrder <- c()
    for (i in seq_len(nrow(BmatrixTime))) {
      rowOrder[i] <- which(BmatrixTime$variable == colnames(Dmatrix[, selectDColumnsTime])[i])
    }
    BmatrixTime <- BmatrixTime[rowOrder, ]
    if (any(colnames(Dmatrix[, selectDColumnsTime]) != BmatrixTime$variable)) {
      log4r::error(self$log, "Column mismatch for time in getEffectMatrix")
      stop()
    }
    AmatrixTime <- as.data.frame(as.matrix(Dmatrix[, selectDColumnsTime]) %*% as.matrix(BmatrixTime[, -1]))
    AmatrixTime$comp <- "TIME"

    BmatrixGroup <- self$regression_coefficients[self$regression_coefficients$comp == "GROUP", c("covar", "estimate", "variable")]
    BmatrixGroup <- dcast(BmatrixGroup, formula = variable ~ covar, value.var = "estimate")
    selectDColumnsGroup <- colnames(Dmatrix) %in% BmatrixGroup$variable
    rowOrder <- c()
    for (i in seq_len(nrow(BmatrixGroup))) {
      rowOrder[i] <- which(BmatrixGroup$variable == colnames(Dmatrix[, selectDColumnsGroup])[i])
    }
    BmatrixGroup <- BmatrixGroup[rowOrder, ]
    if (any(colnames(Dmatrix[, selectDColumnsGroup]) != BmatrixGroup$variable)) {
      log4r::error(self$log, "Column mismatch for group in getEffectMatrix")
      stop()
    }
    AmatrixGroup <- as.data.frame(as.matrix(Dmatrix[, selectDColumnsGroup]) %*% as.matrix(BmatrixGroup[, -1]))
    AmatrixGroup$comp <- "GROUP"
    self$effect_matrix <- rbind(AmatrixTime, AmatrixGroup)
  } else {
    BmatrixTime <- self$regression_coefficients[self$regression_coefficients$comp == "TIME", c("covar", "estimate", "variable")]
    BmatrixTime <- dcast(BmatrixTime, formula = variable ~ covar, value.var = "estimate")

    selectDColumnsTime <- colnames(Dmatrix) %in% BmatrixTime$variable
    rowOrder <- c()
    for (i in seq_len(nrow(BmatrixTime))) {
      rowOrder[i] <- which(BmatrixTime$variable == colnames(Dmatrix[, selectDColumnsTime])[i])
    }
    BmatrixTime <- BmatrixTime[rowOrder, ]
    if (any(colnames(Dmatrix[, selectDColumnsTime]) != BmatrixTime$variable)) {
      log4r::error(self$log, "Column mismatch for time in getEffectMatrix")
      stop()
    }
    AmatrixTime <- as.data.frame(as.matrix(Dmatrix[, selectDColumnsTime]) %*% as.matrix(BmatrixTime[, -1]))
    AmatrixTime$comp <- "TIME"
    self$effect_matrix <- AmatrixTime
  }

  self$parts$time <- parts$time
  if (self$keep_terms != "") {
    keep_terms <- c("group", self$keep_terms)
    self$parts$group <- apply(parts[, ..keep_terms], 1, paste, collapse = " - ")
  } else {
    self$parts$group <- parts$group
  }
  if (!self$minimize_self) {
    log4r::info(self$log, "Finished calculating effect matrix!")
  }
  #invisible(self)
}
