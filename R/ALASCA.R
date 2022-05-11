#' Create an ALASCA model
#'
#' `ALASCA` initializes an ALASCA model and returns an ALASCA object.
#'
#' This function builds your ALASCA model. It needs, at least, a data frame and a formula. The effect matrices can be specified with `effects`, e.g., `c("time", "time+group+time:group", "group+time:group")`.
#'
#' @param df Data frame to be analyzed
#' @param formula Regression formula
#' @param effects Vector with effect terms. If `NULL`, ALASCA will guess (default)
#' @param scale_function Either a custom function or string to define scaling function: `sdall`, `sdref`, `sdt1`, `sdreft1`
#' @param separate_effects Boolean. Separate effects?
#' @param equal_baseline Set to `TRUE` to remove interaction between effects
#' @param validate Boolean. If `TRUE`, give estimates for robustness
#' @param n_validation_runs number of validation runs
#' @param validation_method Choose between `bootstrap` (default) and "jack-knife"
#' @param stratification_column The column to stratify participants by
#' @param save Save models and plots automatically (default: `FALSE`)
#' @param filename File name to save model and plots (when `save = TRUE`)
#' @param use_Rfast Boolean. Defaults to `TRUE`
#' @param p_adjust_method Method for correcting p values for multiple testing, see p.adjust.methods
#' @param participant_column String. Name of the column containing participant identification
#' @param n_validation_folds Partitions when validating
#' 
#' @return An ALASCA object
#'
#' @export
ALASCA <- function(df,
                   formula,
                   effects = NULL,
                   ...) {
  
  object <- AlascaModel$new(df, formula, effects, ...)
  
  # Validate the model ----
  if (object$validate) {
    object$log("Starting validation", level = "DEBUG")
    object$do_validate()
    object$log("Completing validation", level = "DEBUG")
  } else {
    if (object$do_debug) currentTs <- Sys.time()
    object$clean_alasca()
    if (object$do_debug) cat("* clean_alasca:", Sys.time() - currentTs, "s\n")
  }
  
  object$run_time <- Sys.time() - object$init_time
  
  # Save the model
  if (object$save) {
    object$save_model()
  }
  
  if (object$save_to_disk) {
    DBI::dbDisconnect(object$db_con)
  }
  
  object$log("==== ALASCA has finished ====")
  object$finished <- TRUE
  object$log("To visualize the model, try `plot(<object>, effect = 1, component = 1, type = 'effect')`")
  return(object)
}



#' Check if response variables are missing
#' 
#' Called from AlascaDataset
#' 
#' @return An ALASCA object
find_missing_response_variables <- function() {
  if(self$df[, uniqueN(variable), by = .(ID, time)][, uniqueN(V1)] > 1) {
    if (self$model$ignore_missing) {
      self$model$log("Response variables missing for some samples! Continue with caution!", level = "WARN")
    } else {
      self$model$log("Response variables missing for some samples! To ignore this, use `ignore_missing = TRUE`", level = "ERROR")
      stop()
    }
  }
}

#' Check if predictor variables are missing
#' 
#' Called from AlascaDataset
#'
#' @return void
find_missing_predictor_variables <- function() {
  if (any(is.na(self$df))) {
    if (self$model$ignore_missing_covars) {
      self$model$log("Predictor variables missing for some samples! Continue with caution!", level = "WARN")
    } else {
      self$model$log("Predictor variables missing for some samples! To ignore this, use `ignore_missing_covars = TRUE`", level = "ERROR")
      stop()
    }
  }
}

#' Modify formula
#'
#' The regression formula must be adapted to the regression method
#'
#' @return void
get_regression_formula <- function() {
  if (self$model$method %in% c("LMM")) {
    if (self$ID != "ID") {
      
      # The user has specified a column
      self$model$df[, ID := get(self$ID)]
      self$model$participant_column <- "ID"
      self$replace(old_term = self$ID, new_term = "ID")
      
    } else if (self$ID == "ID") {
      
      self$model$log("Using ID for participants!", level = "DEBUG")
      
    } else if (length(self$random_terms) > 1) {
      
      self$model$log("Multiple random effects, couldn't determine participant-id. Please specify `participant_column`", level = "ERROR")
      stop()
      
    } else {
      
      self$model$log("Failed to find ID column", level = "ERROR")
      stop()
      
    }
    
    if (self$model$use_Rfast) {
      # Using Rfast
      self$regression_formula <- self$formula_wo_random
    } else {
      # Using lme4
      self$regression_formula <- formula(paste("value ~ modmat+", paste(self$random_terms, collapse = "+")))
    }
  }else if (self$model$method %in% c("LM")) {
    
    self$regression_formula <- as.formula("value ~ modmat")
    
  } else {
    
    self$model$log("Sorry, an error occurred! Please check your model", level = "ERROR")
    stop()
    
  }
  #invisible(self)
}

#' Determine whether to use LM or LMM
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
set_method <- function() {
  if (is.null(self$model$method)) {
    # Find which default method to use
    
    if (self$has_random()) {
      self$model$method <- "LMM"
      self$model$log("Will use linear mixed models!")
      if (!self$compatible_with_Rfast) {
        self$model$log("Cannot use Rfast in this case. Use lme4 with `use_Rfast = FALSE` instead!", level = "ERROR")
        stop()
      }
    } else {
      self$model$method <- "LM"
      self$model$log("Will use linear models!")
    }
  } else {
    # The user has specified a method to use
    if (self$model$method == "LMM") {
      if (!self$has_random()) {
        self$model$log("The model must contain at least one random effect. Are you sure you wanted linear mixed models?", level = "ERROR")
        stop()
      }
    } else if (self$model$method == "LM") {
      if (self$has_random()) {
        self$model$log("The model contains at least one random effect. Are you sure you wanted linear models?", level = "ERROR")
        stop()
      }
    } else {
      self$model$log("You entered an undefined method. Use `LMM` or `LM`!", level = "ERROR")
      stop()
    }
  }
  
  #invisible(self)
}

#' Remove df from object
#'
#' This function removes unnecessary data
#'
#' @param object An ALASCA object
#' @return An ALASCA object
remove_embedded_data <- function() {
  self$partID <- self$df$ID
  self$bootPartID <- self$df$originalIDbeforeBootstrap
  self$df <- NULL
  self$parts <- NULL
  self$stratification_vector <- NULL
  self$parts_with_variable <- NULL
  self$validation_self <- NULL
  self$effect_list$pca <- NULL
  
  attr(self$new_formula, ".Environment") <- NULL
  attr(self$formula, ".Environment") <- NULL
  #invisible(self)
}

#' Get a scaling function
#'
#' Return scaling function
#'
#' @param scale_function_string String to define scaing function: `sdall`, `sdref`, `sdt1`, `sdreft1`
#' @param scale_function.center Boolean. Mean centering
#' @return A scaling function
get_default_scaling_function <- function() {
  
  if (length(self$effect_terms) < 2) {
    if (self$scale_function %in% c("sdref", "sdreft1"))
    self$log(paste0("The scaling `", self$scale_function, "` has been replaced by `sdt1` as there is only one effect term. This corresponds to the column `", self$effect_terms[[1]], "`"), level = "WARN")
    self$scale_function <- "sdt1"
  }
  
  scale_function_string <- self$scale_function
  scale_function.center <- self$scale_function.center
  
  if (scale_function_string == "sdall") {
    if (scale_function.center) {
      self$scale_function <- function(df) {
        # Scale by the SD of all rows
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all rows
        df[, value := as.double(value)][, value := value / sd(value), by = variable]
      }
    }
  } else if (scale_function_string == "sdref") {
    if (scale_function.center) {
      self$scale_function <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value[get(self$effect_terms[[2]]) == self$get_ref(self$effect_terms[[2]])]), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := value / sd(value[get(self$effect_terms[[2]]) == self$get_ref(self$effect_terms[[2]])]), by = variable]
      }
    }
  } else if (scale_function_string == "sdt1") {
    if (scale_function.center) {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[get(self$effect_terms[[1]]) == self$get_ref(self$effect_terms[[1]])]), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := value / sd(value[get(self$effect_terms[[1]]) == self$get_ref(self$effect_terms[[1]])]), by = variable]
      }
    }
  } else if (scale_function_string == "sdreft1") {
    if (scale_function.center) {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[get(self$effect_terms[[1]]) == self$get_ref(self$effect_terms[[1]]) & get(self$effect_terms[[2]]) == self$get_ref(self$effect_terms[[2]])]), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := value / sd(value[get(self$effect_terms[[1]]) == self$get_ref(self$effect_terms[[1]]) & get(self$effect_terms[[2]]) == self$get_ref(self$effect_terms[[2]])]), by = variable]
      }
    }
  } else {
    self$log("Unknown scaling method. Please use of one the following: `none`, `sdall`, `sdref`, `sdreft1`, `sdt1`", level = "ERROR")
    stop()
  }
  #invisible(scale_function)
}

#' Get a scaling function
#'
#' Return scaling function
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_scaling_function <- function() {
  if (is.function(self$scale_function)) {
    # The user provided a custom function
    if (!self$minimize_object) {
      self$log("Scaling data with custom function...")
    }
  } else if (self$scale_function == "none") {
    # The user do not want to scale
    if (!self$minimize_object) {
      self$log("Not scaling data...", level = "WARN")
    }
    self$scale_function <- identity()
  } else if (is.character(self$scale_function)) {
    # Use a deafult scaling
    if (!self$minimize_object) {
      self$log(paste("Scaling data with",self$scale_function,"..."))
    }
    self$get_default_scaling_function()
  } else {
    self$log("Unknown scaling function", level = "ERROR")
    stop()
  }
  #invisible(self)
}

#' Flip an ALASCA object
#'
#' Changes the sign of loadings and scores
#'
#' @param x An ALASCA object
#' @return The rotated object
#' 
#' @export
flip <- function(x, ...) {
  UseMethod("flip")
}

#' Flip an ALASCA object
#'
#' Changes the sign of loadings and scores
#'
#' @param object An ALASCA object
#' @param component Component(s) to be flipped, use `0` or `NULL` to flip all (default)
#' @param effect Specify effect(s) to be flipped, use `0` or `NULL` to flip all (default)
#' @return The rotated object
#' 
#' @export
flip.AlascaModel <- function(object, effect = 0, component = 0) {
  object$flip(effect_i = effect, component = component)
  invisible(object)
}

#' Summary
#'
#' Print info about the ALASCA object
#'
#' @param object An ALASCA object
#' 
#' @export
summary.AlascaModel <- function(object) {
  cat("---- The log file:\n")
  cat(paste0(readLines(object$log_file), collapse="\n"))
}

#' Summary
#'
#' Gives some general information
#'
#' @param object
#' @return An ALASCA object
get_pca_function <- function () {
  if (is.function(self$pca_function)) {
    self$function.pca <- self$pca_function
  } else if (is.character(self$pca_function)) {
    if (self$pca_function == "prcomp") {
      self$function.pca <- function(df, center = TRUE) {
        prcomp(df, scale = FALSE, center = center)
      }
    } else if (self$pca_function == "irlba") {
      self$function.pca <- function(df, center = TRUE) {
        k <- irlba::prcomp_irlba(df, scale = FALSE, center = center, n = floor(0.9*min(dim(df))))
        rownames(k$rotation) <- colnames(df)
        k
      }
    } else if (self$pca_function == "princomp") {
      self$function.pca <- function(df, center = TRUE) {
        k <- princomp(df)
        l <- k$loadings
        s <- k$scores
        out <- list(
          rotation = data.frame(matrix(as.numeric(l), attributes(l)$dim, dimnames=attributes(l)$dimnames)),
          x = data.frame(matrix(as.numeric(s), attributes(s)$dim, dimnames=attributes(s)$dimnames)),
          sdev = k$sdev
        )
        colnames(out$x) <- gsub("Comp.", "PC", colnames(out$x), fixed = TRUE)
        colnames(out$rotation) <- gsub("Comp.", "PC", colnames(out$rotation), fixed = TRUE)
        return(out)
      }
    } else {
      self$log("Unknown PCA function", level = "ERROR")
      stop()
    }
  } else {
    self$log("Unknown PCA function", level = "ERROR")
    stop()
  }
  #invisible(self)
}

#' Organize the ALASCA model construction
#'
#' This function builds the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
build_model <- function() {
  if (!self$minimize_object) {
    # This is not a validation run
    self$log("Starting to build model", level = "DEBUG")
    self$log(paste0("Calculating ", self$method, " coefficients"))
  }
  
  if (self$reduce_dimensions) {
    self$log("Starting to reduce dimensions", level = "DEBUG")
    self$do_reduce_dimensions()
    self$log("-> Finished to reduce dimensions", level = "DEBUG")
  }
  
  if (self$do_debug) currentTs <- Sys.time()
  self$run_regression()
  self$log("Starting to calculate regression coefficients", level = "DEBUG")
  if (!self$use_Rfast) {
    # With Rfast, we've already got the coefficients
    self$get_regression_coefficients()
  }
  
  if (!self$minimize_object) {
    # This is not a validation run
    self$log("-> Finished calculating regression coefficients!", level = "DEBUG")
  }
  
  self$remove_covars()
  self$get_effect_matrix()
  self$do_pca()
  self$clean_pca()
  
  if (self$method %in% c("LM", "LMM")) {
    if (self$minimize_object) {
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
  self$log("Find rows by variable", level = "DEBUG")
  
  # https://stackoverflow.com/questions/71946874/what-is-the-most-efficient-method-for-finding-row-indices-by-group-in-a-data-tab
  if (nrow(self[["df"]]) > 1000) {
    tmp <- self[["df"]][, .(idx = .(.I)), variable]
    rows_by_variable <- tmp$idx
    names(rows_by_variable) <- tmp$variable
  } else {
    rows_by_variable <- lapply(self$get_levels("variable"), function(x) data.table(I = which(self[["df"]][["variable"]] == x)) )
    names(rows_by_variable) <- self$get_levels("variable")
  }
  
  
  
  #names(rows_by_variable) <- self$get_levels("variable")
  
  if (self$use_Rfast && self$method %in% c("LMM")) {
    # start.time <- Sys.time()
    if (any(is.na(self[["df"]][, value]))) {
      self$log("Rfast does NOT like NA's! Check your scaling function or value column.", level = "ERROR")
      stop()
    }
    self$log("Make model matrix", level = "DEBUG")
    if (!self[["minimize_object"]] || self$reduce_dimensions) {
      self[["modmat"]] <- model.matrix(self[["formula"]][["regression_formula"]], data = self[["df"]])
      if (self[["equal_baseline"]]) {
        self[["modmat"]] <- self[["modmat"]][, !grepl(paste0(self$effect_terms[[1]], self$get_ref(self$effect_terms[[1]])), colnames(self[["modmat"]]))]
      }
      self[["cnames_modmat"]] <- colnames(self[["modmat"]])
    }
    self$log("-> Finished model matrix", level = "DEBUG")
    
    self$log("Starting regression", level = "DEBUG")
    # https://stackoverflow.com/questions/61013078/fastest-way-to-convert-a-matrix-to-a-data-table
    self$regression_coefficients <- setDT(as.data.frame(
      vapply(self$get_levels("variable"), function(x) {
        Rfast::rint.reg(
          y = self[["df"]][ rows_by_variable[[x]], value],
          x = self[["modmat"]][rows_by_variable[[x]], -1],
          id = as.numeric(factor(self[["df"]][ rows_by_variable[[x]], ID])),
          ranef = FALSE
        )$be
      }, FUN.VALUE = numeric(ncol(self[["modmat"]])))
    ))[]
    self$log("-> Finished regression", level = "DEBUG")
    colnames(self$regression_coefficients) <- self$get_levels("variable")
    self$regression_coefficients[, variable := self[["cnames_modmat"]]]
    self$regression_coefficients <- melt(self$regression_coefficients, id.vars = "variable", variable.name = "covar", variable.factor = FALSE, value.name = "estimate")
    self$regression_coefficients[, pvalue := NA]
    
    # end.time <- Sys.time()
    # cat("\n\n",end.time - start.time,"\n")
    #invisible(self)
  } else if (!self$use_Rfast && self$method %in% c("LM")) {
    
    #' We need to modify the model matrix. Therefore, self[["formula"]][["formula"]] contains the original terms, whereas
    #' self[["formula"]][["regression_formula"]] has replaced the terms with the modified model matrix
    
    self$regression_model <- lapply(self$get_levels("variable"), function(x) {
      modmat <- model.matrix(self[["formula"]][["formula"]], data = self$df[ rows_by_variable[[x]] ])
      #self[["modmat"]] <- self[["modmat"]][, -1] # Remove intercept
      if (self[["equal_baseline"]]) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0(self$effect_terms[[1]], self$get_ref(self$effect_terms[[1]])), colnames(modmat))]
      }
      self[["cnames_modmat"]] <- colnames(modmat)
      environment(self[["formula"]][["regression_formula"]]) <- environment()
      regression_model <- lm(self[["formula"]][["regression_formula"]], data = self$df[ rows_by_variable[[x]] ])
      attr(regression_model, "name") <- x
      regression_model
    })
    names(self$regression_model) <- self$get_levels("variable")
  } else if (self$use_Rfast && self$method %in% c("LM")) {
    # start.time <- Sys.time()
    if (any(is.na(self$df[, value]))) {
      self$log("Rfast does NOT like NA's! Check your scaling function or value column.", level = "ERROR")
      stop()
    }
    
    if (!self$minimize_object || self$reduce_dimensions) {
      self[["modmat"]] <- model.matrix(self[["formula"]][["formula"]], data = self$df)
      if (self[["equal_baseline"]]) {
        # Remove interaction between group and first time point
        self[["modmat"]] <- self[["modmat"]][, !grepl(paste0(self$effect_terms[[1]], self$get_ref(self$effect_terms[[1]])), colnames(self[["modmat"]]))]
      }
      self[["cnames_modmat"]] <- colnames(self[["modmat"]])
    }
    
    self$regression_coefficients <- setDT(as.data.frame(
      vapply(self$get_levels("variable"), function(x) {
        Rfast::lmfit(
          y = self[["df"]][ rows_by_variable[[x]], value],
          x = self[["modmat"]][rows_by_variable[[x]], ]
        )$be
      }, FUN.VALUE = numeric(ncol(self[["modmat"]])))
    ))[]
    self$log("-> Finished regression", level = "DEBUG")
    colnames(self$regression_coefficients) <- self$get_levels("variable")
    self$regression_coefficients[, variable := self[["cnames_modmat"]]]
    self$regression_coefficients <- melt(self$regression_coefficients, id.vars = "variable", variable.name = "covar", variable.factor = FALSE, value.name = "estimate")
    self$regression_coefficients[, pvalue := NA]

    # end.time <- Sys.time()
    # cat("\n\n",end.time - start.time,"\n")
    #invisible(self)
  } else if (!self$use_Rfast && self$method %in% c("LMM")) {
    
    #' We need to modify the model matrix. Therefore, self[["formula"]][["formula"]] contains the original terms, whereas
    #' self[["formula"]][["regression_formula"]] has replaced the terms with the modified model matrix
    
    self$regression_model <- lapply(self$get_levels("variable"), function(x) {
      modmat <- model.matrix(self[["formula"]][["formula"]], data = self$df[ rows_by_variable[[x]] ])
      #odmat <- modmat[, -1] # Remove intercept
      if (self[["equal_baseline"]]) {
        # Remove interaction between group and first time point
        modmat <- modmat[, !grepl(paste0(self$effect_terms[[1]], self$get_ref(self$effect_terms[[1]])), colnames(modmat), fixed = TRUE)]
      }
      self[["cnames_modmat"]] <- colnames(modmat)
      # modmat <- modmat[,ncol(modmat):1]
      environment(self[["formula"]][["regression_formula"]]) <- environment()
      regression_model <- lmerTest::lmer(self[["formula"]][["regression_formula"]], data = self$df[ rows_by_variable[[x]] ])
      attr(regression_model, "name") <- x
      regression_model
    })
    names(self$regression_model) <- self$get_levels("variable")
  }
  
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
  
  if (!is.null(self$p_adjust_method)) {
    self$log("Adjusting p values")
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
  if (self$formula$has_covars()) {
    
    self$covar_coefficients <- rbindlist(lapply(unique(self$formula$covars), function(x) {
      to_list <- self$regression_coefficients[substr(variable, 1, nchar(x)) == x]
      self$regression_coefficients <- self$regression_coefficients[substr(variable, 1, nchar(x)) != x]
      to_list
    }))
    
    if (self$reduce_dimensions) {
      self$covar_coefficients <- rbindlist(lapply(unique(self$covar_coefficients$variable), function(v){
        ref <- self$covar_coefficients[variable == v]
        list(
          variable = v,
          covar = self$reduced_df$loading$covars,
          pvalue = NA,
          estimate = as.matrix(self$reduced_df$loading[, -"covars"]) %*% as.matrix(ref$estimate[match(colnames(self$reduced_df$loading[,-"covars"]), ref$covar)])
        )
      }))
    }
    
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
  if (!self$minimize_object) {
    self$log("Calculating effect matrix", level = "DEBUG")
  }
  #stop()
  reg_coefs <- dcast(self$regression_coefficients, variable~covar, value.var = "estimate")
  rownames(reg_coefs) <- reg_coefs$variable
  self$effect_list$effect_matrix <- lapply(self$set_design_matrices(), function(mm) {
    mm[self$df[variable == self$get_ref("variable"), which = TRUE ], ] %*% as.matrix(reg_coefs[colnames(mm), -1])
  }
  )
  
  if (!self$minimize_object) {
    self$log("-> Finished calculating effect matrix!", level = "DEBUG")
  }
  #invisible(self)
}

#' Save ALASCA model
#'
#' @param object An ALASCA model
#' @return void
#' @export
save.AlascaModel <- function(object) {
  object$save_model()
}

#' Perform PCA
#'
#' This function performs PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
do_pca <- function() {
  
  self$effect_list$pca <- lapply(
    self$effect_list$effect_matrix,
    function(x) {
      k <- prcomp(x, scale = FALSE, center = !self$scale_function.center)
      
      if (ncol(k$x) > self$max_PC) {
        if(!self$minimize_object) self$log(paste("Keeping", self$max_PC, "of", ncol(k$x), "components. Change `max_PC` if necessary."), level = "WARN")
        k$x <- k$x[, seq_len(self$max_PC)]
        k$rotation <- k$rotation[, seq_len(self$max_PC)]
      }
      k
    }
  )
  #invisible(self)
}

#' Perform PCA for Limm-PCA
#'
#' This function performs PCA before "the real PCA"
#'
#' @param object An ALASCA object
#' @return An ALASCA object
do_reduce_dimensions <- function(){
  
  if (!self$minimize_object) {
    self$log("Reducing the number of dimensions with PCA")
  }
  
  wide_data <- dcast(data = self$df, paste(paste(self$formula$all_formula_terms, collapse = " + "), "~ variable"), value.var = "value")
  
  temp_pca_values <- self$function.pca(
    wide_data[, .SD, .SDcols = -self$formula$all_formula_terms],
    center = !self$scale_function.center
  )
  
  self$reduced_df[["explanatory_power"]] <- temp_pca_values$sdev^2 / sum(temp_pca_values$sdev^2)
  
  # Remove surplus columns
  if (is.null(self$reduced_df[["nComps"]])) {
    self$reduced_df[["nComps"]] <- which(cumsum(self$reduced_df[["explanatory_power"]]) >= self$reduced_df[["limit"]])[1]
    if (!self$minimize_object) {
      self$log(paste("Keeping",
                     self$reduced_df[["nComps"]],
                     "components from initial PCA, explaining",
                     round(100*cumsum(self$reduced_df[["explanatory_power"]])[self$reduced_df[["nComps"]]], 2),
                     "% of variation. The limit can be changed with `reduce_dimensions.limit`"))
    }
  }
  
  if(ncol(temp_pca_values$rotation) > self$reduced_df[["nComps"]]){
    temp_pca_values$rotation <- temp_pca_values$rotation[, seq_len(self$reduced_df[["nComps"]])]
    temp_pca_values$x <- temp_pca_values$x[, seq_len(self$reduced_df[["nComps"]])]
  }
  
  temp_pca_values$rotation <- data.table(temp_pca_values$rotation, keep.rownames = "covars")
  setkey(temp_pca_values$rotation, covars)
  
  if (is.null(self$df_raw$reduced_df$loading)) {
    
    # Make a copy of the main model loadings to use as reference
    self$df_raw$reduced_df$loading <- temp_pca_values$rotation
    
    
  } else {
    
    # Check if the pca model needs reflection to better fit the main model
    
    for (i in seq(2, ncol(temp_pca_values$rotation))) {
      V1 <- sum((temp_pca_values$rotation[self$df_raw$reduced_df$loading$covars, ..i] - self$df_raw$reduced_df$loading[, ..i])^2)
      V2 <- sum((-temp_pca_values$rotation[self$df_raw$reduced_df$loading$covars, ..i] - self$df_raw$reduced_df$loading[, ..i])^2)
      if(V2 < V1){
        temp_pca_values$rotation[, i] <- -temp_pca_values$rotation[, ..i]
        temp_pca_values$x[, i-1] <- -temp_pca_values$x[, i-1]
      }
    }
    
  }
  
  self$reduced_df$loading <- temp_pca_values$rotation
  self$reduced_df$score <- temp_pca_values$x
  self$reduced_df$df <- self$df
  self$df <- melt(data = cbind(wide_data[, .SD, .SDcols = self$formula$all_formula_terms], self$reduced_df$score),
                  id.vars = self$formula$all_formula_terms, variable.factor = FALSE)
  if (is.null(self$splot$loading_group_column)) {
    self$df <- merge(self$df,
                     unique(self$reduced_df$df[, .SD, .SDcols = self$formula$all_terms]),
                     by = self$formula$all_formula_terms,
                     all = TRUE)
  } else {
    self$df <- merge(self$df,
                     unique(self$reduced_df$df[, .SD, .SDcols = self$formula$all_terms[
                          self$formula$all_terms != self$splot$loading_group_column]
                        ]
                      ),
      by = self$formula$all_formula_terms,
      all = TRUE)
    if (any(is.na(self$df$value))) {
      self$log("Something went wrong. Check your data set for duplicated values", level = "ERROR")
      stop()
    }
  }
  
  #if (!self$minimize_object) {
  #  self$raw_data$modmat
  #}
  self$reduced_df[["variables"]] <- unique(self$df$variable)
  self$stratification_vector <- self$df[, get(self$stratification_column)]
  if (!self$minimize_object) {
    self$log("-> Finished the reduction of dimensions!")
  }
  
  #invisible(self)
}

#' Clean the PCA data
#'
#' This function makes the pca output more useful
#'
clean_pca <- function() {
  self$log("Starting clean_pca", level = "DEBUG")
  
  self$ALASCA$score <- lapply(seq_along(self$effect_list$pca), function(i){
    unique(
      cbind(
        self$df[variable == self$get_ref("variable"), .SD, .SDcols = self$effect_list$terms[[i]]],
        self$effect_list$pca[[i]]$x
      )
    )
  })
  
  self$ALASCA$loading <- lapply(seq_along(self$effect_list$pca), function(i){
    data.table(self$effect_list$pca[[i]]$rotation, keep.rownames = "covars")
  })
  
  self$ALASCA$explained <- lapply(self$effect_list$pca, function(x){
    x$sdev^2 / sum(x$sdev^2)
  })
  
  self$ALASCA$significant_PCs <- lapply(self$ALASCA$explained, function(x){
    unique(c(1, 2, which(x >= self$explanatory_limit)))
  })
  
  if(self$reduce_dimensions){
    # Loadings must be back-transformed
    
    for (i in seq_along(self$effect_list$pca)) {
      setkey(self$ALASCA$loading[[i]], covars)
      self$ALASCA$loading[[i]] <- data.table(
        covars = self$reduced_df$loading[, covars],
        as.matrix(self$reduced_df$loading[, -1]) %*% as.matrix(self$ALASCA$loading[[i]][colnames(self$reduced_df$loading[, -1]), -1])
      )
    }
    
    setkey(self$regression_coefficients, covar)
    
    self$regression_coefficients <- rbindlist(
      lapply(unique(self$regression_coefficients$variable), function(x){
        list(
        variable = x,
        pvalue = NA,
        covar = self$reduced_df$loading[, covars],
        estimate = as.matrix(self$reduced_df$loading[, -1]) %*% as.matrix(
          self$regression_coefficients[variable == x][colnames(self$reduced_df$loading[, -1])][, estimate]
          )
          )
      })
    )
    
    self$reduced_df$df <- NULL
    self$reduced_df$loading <- NULL
    self$reduced_df$score <- NULL
    self$reduced_df$nComps <- NULL
    
  }
  
  # Ensure that the highest loading has positive sign
  for (i in seq_along(self$ALASCA$loading)) {
    setkeyv(self$ALASCA$loading[[i]], cols = "covars")
    setkeyv(self$ALASCA$score[[i]], cols = self$effect_list$terms[[i]])
    
    # Get columns where the sign needs to be changed
    nVar <- max.col(abs(t(self$ALASCA$loading[[i]][, -1])))
    sVar <- sign(diag(as.matrix(self$ALASCA$loading[[i]][nVar, -1])))
    cols_to_change <- colnames(self$ALASCA$loading[[i]])[1+which(sVar < 0)]
    
    # Change column signs
    self$ALASCA$loading[[i]][, (cols_to_change) := -.SD, .SDcols = cols_to_change]
    self$ALASCA$score[[i]][, (cols_to_change) := -.SD, .SDcols = cols_to_change]
  }
  
  self$log("Completed clean_pca", level = "DEBUG")
  
  #invisible(self)
}

#' Clean the ALASCA data
#'
#' This function makes the pca output more useful
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
clean_alasca <- function() {
  
  for (i in seq_along(self$ALASCA$loading)) {
    self$ALASCA$score[[i]] <- melt(self$ALASCA$score[[i]],
                                   id.vars = self$effect_list$terms[[i]],
                                   variable.name = "PC",
                                   value.name = "score",
                                   variable.factor = FALSE)
    self$ALASCA$score[[i]][, PC := as.integer(substr(PC, 3, nchar(PC))), ]
    setkeyv(self$ALASCA$score[[i]], cols = self$effect_list$terms[[i]])
    
    self$ALASCA$loading[[i]] <- melt(self$ALASCA$loading[[i]],
                                     id.vars = "covars",
                                     variable.name = "PC",
                                     value.name = "loading",
                                     variable.factor = FALSE)
    self$ALASCA$loading[[i]][,   PC := as.integer(substr(PC, 3, nchar(PC))), ]
    setkeyv(self$ALASCA$loading[[i]], cols = "covars")
  }
  
  #invisible(self)
}

#' Validate the ALASCA model LMM
#'
#' This function performs leave-one-out robustness testing of your ALASCA model. If you didn't specify the number of runs `n_validation_runs` when initializing the model (see \code{\link{ALASCA}}), you can do it by running for example `model$n_validation_runs <- 100` prior to calling `validate`. Your dataset is divided into `n_validation_folds` partitions, keeping group proportions, and one of these are left out. `n_validation_folds` is set the same way as  `n_validation_runs`.
#'
#' @param object An ALASCA object
#' @param participant_column The name of the column containing participant identifier. Needed if not set during initialization of the model.
#' @param validate_regression Whether to validate regression models
#' @return An ALASCA object
#'
#' @examples
#' load("PE.Rdata")
#' model$n_validation_runs <- 10
#' model.val <- validate(model, participant_column = "ID")
#' @export
do_validate <- function() {
  if (self$validate) {
    # stop("The object has already been validated")
  }
  
  self$log("Starting validation")
  
  start_time_all <- Sys.time()
  self$get_validation_ids()
  
  temp_object <- lapply(seq_len(self$n_validation_runs), FUN = function(ii) {
    self$log(paste0("- Run ", ii, " of ", self$n_validation_runs))
    start.time.this <- Sys.time()
    
    # Make resampled model
    temp_object <- self$prepare_validation_run(runN = ii)
    if (nrow(self$df) > 100000) gc()
    
    # Rotate new loadings/scores to the original model
    if (self$optimize_score) {
      temp_object$rotate_matrix_optimize_score(target = self)
    } else {
      temp_object$rotate_matrix(target = self)
    }
    
    temp_object$clean_alasca()
    
    if (self$save_to_disk) {
      temp_object$save_validation(ii)
    }
    temp_object$effect_list <- NULL
    
    time_all <- difftime(Sys.time(), start_time_all, units = c("secs")) / ii
    self$log(paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((self$n_validation_runs - ii) * time_all, 2), " seconds"))
    temp_object
  })
  
  self$clean_alasca()
  
  self$log("Calculating percentiles for score and loading")
  self$get_validation_percentiles(objectlist = temp_object)
  
  #invisible(object)
}

.procrustes <- function(loadings, target) {
  s <- t(loadings) %*% target
  w1 <- s %*% t(s)
  v1 <- t(s) %*% s
  w <- eigen(w1)$vectors
  ew <- diag(eigen(w1)$values)
  v <- eigen(v1)$vectors
  ev <- diag(eigen(v1)$values)
  o <- t(w) %*% s %*% v
  k <- diag(((diag(o)) / abs(diag(o))), nrow = nrow(o), ncol = nrow(o))
  ww <- w %*% k
  out <- list()
  out$t1 <- ww %*% t(v) # Rotation matrix
  out$procrust <- loadings %*% out$t1 # Rotated loadings
  return(out)
}

#' Extract percentiles
#'
#' This function extract percentiles during validation
#'
#' @param objectlist List of ALASCA objects
get_validation_percentiles <- function(objectlist) {
  
  self$get_validation_percentiles_loading(objectlist)
  self$get_validation_percentiles_score(objectlist)
  if (self$validate_regression) {
    self$get_validation_percentiles_regression(objectlist)
  }
  if (self$formula$has_covars()) {
    self$get_validation_percentiles_covars(objectlist)
  }
  
  #invisible(self)
}

#' Extract percentiles for regressions
#'
#' This function extract percentiles for validation of regression
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_regression <- function(objectlist) {
  
  if (self$save_to_disk) {
    res <- DBI::dbSendQuery(self$db_con, paste0("SELECT prediction.", paste0(self$effect_terms, collapse = ", prediction."), ", prediction.pred, variables.variable AS variable FROM prediction
                            LEFT JOIN variables
                            ON variables.id = prediction.variable"))
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df <- rbindlist(lapply(objectlist, function(x) x$model_prediction))
  }
  
  df <- df[, as.list(quantile(pred, probs = self$limitsCI, type = self$validation_quantile_method)),
           by = c(self$effect_terms, "variable")]
  colnames(df) <- c(self$effect_terms, "variable", "low", "high")
  
  self$model_prediction <- merge(self$model_prediction, df)
  #invisible(self)
}

#' Extract percentiles for covariates
#'
#' This function extract percentiles for validation of covariates
#'
#' @inheritParams get_validation_percentiles
get_validation_percentiles_covars <- function(objectlist) {
  
  if (self$save_to_disk) {
    res <- DBI::dbSendQuery(self$db_con, paste0("SELECT covars.variable, covars.estimate, variables.variable AS covar FROM covars
                                                LEFT JOIN variables
                                                ON variables.id = covars.covar"))
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df <- rbindlist(lapply(objectlist, function(x) x$get_covars()))
  }
  df <- df[, as.list(quantile(estimate, probs = self$limitsCI, type = self$validation_quantile_method)),
           by = .(covar, variable)]
  colnames(df) <- c("covar", "variable", "low", "high")
  
  self$covar_coefficients <- merge(self$covar_coefficients, df)
  #invisible(self)
}

#' Extract percentiles for loading
#'
#' This function extract percentiles during validation of loadings
#'
#' @inheritParams get_validation_percentiles
get_validation_percentiles_loading <- function(objectlist) {
  
  for (effect_i in seq_along(self$ALASCA$loading)){
    df <- self$get_validation_loadings(objectlist, effect_i = effect_i)
    tmp <- df[, as.list(
      quantile(loading, probs = self$limitsCI, type = self$validation_quantile_method)
    ), by = .(PC, covars)]
    colnames(tmp) <- c("PC", "covars", "low", "high")
    self$ALASCA$loading[[effect_i]] <- merge(self$ALASCA$loading[[effect_i]], tmp,  all.x = TRUE, by = c("PC", "covars"))
  }
  #invisible(self)
}


get_validation_scores <- function(objectlist = NULL, effect_i = 1) {
  if (self$save_to_disk) {
    res <- DBI::dbSendQuery(self$db_con,
                            paste0("SELECT * FROM score WHERE effect = ", effect_i, " AND PC IN(", paste(self$get_PCs(effect_i), collapse = ", "), ")"))
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    if (is.null(self$effect_list$saved_scores)) {
      self$log("Saving validation scores", level = "DEBUG")
      for (effect_k in seq_along(self$effect_list$expr)) {
        df <- rbindlist(
          lapply(seq_along(objectlist), function(x) {
            data.table(
              model = x,
              objectlist[[x]]$ALASCA$score[[effect_k]][PC %in% self$get_PCs(effect_k), ]
            )
          }), fill = TRUE)
        fname <- paste0(self$filepath, self$filename, "_validation_scores_effect_", effect_k, ".fst")
        fst::write_fst(df, fname, compress = self$compress_validation)
        self$effect_list$saved_scores[[effect_k]] <- fname
      }
    } 
    df <- fst::read_fst(self$effect_list$saved_scores[[effect_i]], as.data.table = TRUE)
  }
  return(df)
}

get_validation_loadings <- function(objectlist = NULL, effect_i = 1) {
  if (self$save_to_disk) {
    res <- DBI::dbSendQuery(self$db_con,
                            paste0("SELECT loading.PC, loading.loading, variables.variable AS covars FROM loading
                              LEFT JOIN variables
                              ON variables.id = loading.covars
                              WHERE effect = ", effect_i, " AND PC IN(", paste(self$get_PCs(effect_i), collapse = ", "), ")"))
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    
    if (is.null(self$effect_list$saved_loadings)) {
      self$log("Saving validation loadings", level = "DEBUG")
      for (effect_k in seq_along(self$effect_list$expr)) {
        df <- rbindlist(
          lapply(seq_along(objectlist), function(x) {
            data.table(
              model = x,
              objectlist[[x]]$ALASCA$loading[[effect_k]][PC %in% self$get_PCs(effect_k), ]
            )
          }),
          fill = TRUE
        )
        fname <- paste0(self$filepath, self$filename, "_validation_loadings_effect_", effect_k, ".fst")
        fst::write_fst(df, fname, compress = self$compress_validation)
        self$effect_list$saved_loadings[[effect_k]] <- fname
      }
    }
    
    df <- fst::read_fst(self$effect_list$saved_loadings[[effect_i]], as.data.table = TRUE)
  }
  return(df)
}

#' Extract percentiles for scores
#'
#' This function extract percentiles during validation of scores
#'
#' @inheritParams get_validation_percentiles
get_validation_percentiles_score <- function(objectlist) {
  
  for (effect_i in seq_along(self$ALASCA$score)) {
    
    all_scores <- self$get_validation_scores(objectlist = objectlist, effect_i = effect_i)
    
    tmp <- all_scores[, as.list(
      quantile(score, probs = self$limitsCI, type = self$validation_quantile_method)
    ), by = c("PC", self$effect_list$terms[[effect_i]])]
    colnames(tmp) <- c("PC", self$effect_list$terms[[effect_i]], "low", "high")
    
    self$ALASCA$score[[effect_i]] <- merge(self$ALASCA$score[[effect_i]], tmp, all.x = TRUE, by = c("PC",  self$effect_list$terms[[effect_i]]))
  }
  
  #invisible(self)
}

#' Validate underlying regression models
#'
#' This function calcuates predictions from each regression model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_regression_predictions <- function() {
  if (!self$minimize_object) {
    # This is not a validation run
    self$log("Calculating predictions from regression models", level = "DEBUG")
  }
  
  regCoeffAll <- dcast(data = self[["regression_coefficients"]], variable~covar, value.var = "estimate")
  rownames(regCoeffAll) <- regCoeffAll$variable
  regModel <- unique(model.matrix(self$formula$formula_wo_random, data = self$df))
  if (self$equal_baseline) {
    regModel <- regModel[, !grepl(paste0(self$formula$all_terms[[1]], self$get_ref(self$formula$all_terms[[1]])), colnames(regModel), fixed = TRUE), drop = FALSE]
  }
  regModel <- regModel[, grepl(paste0(self$effect_terms, collapse = "|"), colnames(regModel)), drop = FALSE]
  regModel <- unique(regModel)
  
  self$model_prediction <- melt(
    cbind(
      as.matrix(regModel) %*% as.matrix(regCoeffAll[colnames(regModel), -1]),
      self$df[as.numeric(rownames(regModel)), .SD, .SDcols = c(self$effect_terms)]
    ),
    id.vars = c(self$effect_terms), value.name = "pred"
  )
  
  if (!self$minimize_object) {
    # This is not a validation run
    self$log("-> Finished calculating predictions from regression models!", level = "DEBUG")
  }
  #invisible(self)
}

get_validation_ids <- function() {
  if (is.null(self$validation_ids)) {
    self$log("Generating random validation sample", level = "DEBUG")
    
    original_IDs <- unique(self$df_raw$df[, .SD, .SDcols = c(self$formula$ID, self$stratification_column)])
    colnames(original_IDs) <- c("ID", "group")
    
    if (self$validation_method == "bootstrap") {
      
      tmp <- lapply(unique(self$stratification_vector), function(strat_group) {
        IDs_to_choose_from <- original_IDs[group == strat_group, ID]
        t(
          vapply(
            seq_len(self$n_validation_runs),
            function(x) sample(IDs_to_choose_from, replace = TRUE),
            FUN.VALUE = integer(length(sample(IDs_to_choose_from, replace = TRUE)))
            )
          )
      })
      
      self$validation_ids <- do.call(cbind, tmp)
      
    } else {
      
      tmp <- lapply(unique(self$stratification_vector), function(strat_group) {
        IDs_to_choose_from <- original_IDs[group == strat_group, ID]
        n_to_choose <- floor(length(IDs_to_choose_from) - length(IDs_to_choose_from)/self$n_validation_folds)
        if (n_to_choose <= self$n_validation_folds) {
          self$log(
              paste0("The stratification group ", strat_group, " has only ", length(IDs_to_choose_from), " members. Choosing ", n_to_choose, " of them. Consider adjusting `n_validation_folds`"),
              level = "WARN"
            )
        }
        t(
          vapply(
            seq_len(self$n_validation_runs),
            function(x) sample(IDs_to_choose_from, size = n_to_choose),
            FUN.VALUE = integer(n_to_choose)
          )
        )
      })
      
      self$validation_ids <- do.call(cbind, tmp)
      
    }
    if (self$save_validation_ids) {
      self$log("Saving validation IDs", level = "DEBUG")
      fwrite(as.data.table(self$validation_ids), file = paste0(self$filepath, "validation_IDs.csv"), col.names = FALSE)
    }
  }
}

#' Make a single validation run
#'
#' This function ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
prepare_validation_run <- function(runN) {
  if (self$validation_method %in% c("loo", "jack-knife", "jackknife")) {
    # Use leave-one-out validation
    
    bootobject <- self$clone()
    bootobject$my_df_rows <- unlist(lapply(self$validation_ids[runN, ], function(x) self$df_raw$rows_by_ID[[as.character(x)]]))
    bootobject[["df"]] <- NULL
    bootobject[["effect_list"]][["model_matrix"]] <- NULL
    bootobject[["df"]] <- self$scale_function(self$df_raw$df[bootobject$my_df_rows])
    bootobject[["modmat"]] <- self$modmat[bootobject$my_df_rows,]
    bootobject$update()
    return(bootobject)

  } else if (self$validation_method == "bootstrap") {
    # Use bootstrap validation
    # When using bootstrapping, we resample participants with replacement
    
    self$log("Starting preparation of bootstrap sample", level = "DEBUG")
    
    # Create bootstrap object without data
    bootobject <- self$clone()
    bootobject[["df"]] <- NULL
    bootobject[["modmat"]] <- NULL
    bootobject[["effect_list"]][["model_matrix"]] <- NULL
    bootobject$get_bootstrap_data(df_raw = self$df_raw,
                                  participants_in_bootstrap = self$validation_ids[runN, ],
                                  modmat = self$modmat)
    
    if (self$validation_assign_new_ids) {
      bootobject$df[, ID := uniqueIDforBootstrap] # Replace ID
    }
    
    self$log("-> Completed preparation of bootstrap sample", level = "DEBUG")
    
    bootobject$update()
    bootobject$validation$original_ids <- self$validation_ids[runN, ]
    return(bootobject)
  }
}

#' Make bootstrap data set
#'
#' Get data frame
#'
#' @param object An ALASCA object
#' @return A data frame
get_bootstrap_data <- function(df_raw, participants_in_bootstrap, modmat) {
  
  selected_rows <- rbindlist(
    lapply(seq_along(participants_in_bootstrap), function(participant){
      list(
        new_id = participant,
        row_nr = df_raw$rows_by_ID[[as.character(participants_in_bootstrap[participant])]]
      )
    })
  )
  
  self[["my_df_rows"]] <- selected_rows$row_nr
  self$df <- copy(self$scale_function(df_raw$df[self$my_df_rows,]))
  self$log(paste("Length df:", nrow(self[["df"]])), level = "DEBUG")
  self$log(paste("First 10 rows:", paste0(selected_rows$row_nr[1:10], collapse = " - ")), level = "DEBUG")
  self$df[, uniqueIDforBootstrap := selected_rows$new_id]
  self$df[, originalIDbeforeBootstrap := get(self$formula$ID)]
  if (self$reduce_dimensions) {
    self$modmat <- NULL
  } else {
    self$modmat <- modmat[selected_rows$row_nr,]
    self$log(paste("Length mm:", nrow(self[["modmat"]])), level = "DEBUG")
  }
  
  #invisible(self)
}



#' Get screeplot
#'
#' This function returns a screeplot for an ALASCA model showing what proportion of the variance each component of the model explains
#'
#' @param object An ALASCA object
#' @param component The highest principal component to plot
#'
#' @return An ggplot2 object
#'
#' @export
screeplot.AlascaModel <- function(object, component = 10, ...) {
  object$splot$call_plot(type = "scree", ...)
}

#' Get screeplot
#' 
#' @inheritParams screeplot.AlascaModel
#' 
#' @return An ggplot2 object
#' 
#' @export
screeplot <- function(object, component = 10, ...) {
  UseMethod("screeplot")
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param inheritParams plot_parts
#' @return A ggplot2 object
#'
#' @export
plot_participants <- function(...){
  plot_parts(...)
}

#' Plot participants
#'
#' This function returns the scores for an ALASCA model
#'
#' @param object An ALASCA object or a data frame. If a data frame, you need to specify the column names for participant and value. This also applies if you have not specified the participant column in the ALASCA model before.
#' @param variable List of variable names to print. If `NA`, return all (default).
#' @param participant_column Specify the column with participant identifier. Not necessary if you have already provided it to the ALASCA object
#' @param valueColumn Specify column with values (y axis). Not necessary to provide if you are plotting an ALASCA object.
#' @param timeColumn Specify column with times (x axis). Defaults to `time`.
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param addSmooth. Specify which geom_smooth model you want to apply, eg. `lm`, `glm`, `gam`, `loess` (default). Set to `NA` to remove.
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#' 
#' @export
plot_parts <- function(object,
                       variables = NA,
                       participant_column = FALSE,
                       x_label = NA,
                       group_label = NA,
                       valueColumn = FALSE,
                       timeColumn = "time",
                       addSmooth = "loess",
                       filename = NA,
                       filetype = object$plot.filetype,
                       figunit = object$plot.figunit,
                       plot.ylabel = "value",
                       as.list = FALSE,
                       figsize = object$plot.figsize,
                       my_theme = object$plot.my_theme) {
  
  if (!is.na(filename)) object$filename <- filename
  
  if(as.list){
    if (is.data.frame(object)) {
      df <- object
      if (any(participant_column == FALSE) | any(valueColumn == FALSE)) {
        object$log("You need to specify participant and value columns", level = "ERROR")
        stop()
      } else {
        participant_column <- participant_column
        valueColumn <- valueColumn
      }
      plotFunction <- function(df, timeColumn, valueColumn, participant_column, xi, addSmooth, x_label, group_label, my_theme) {
        g <- ggplot2::ggplot(subset(df, variable == xi), ggplot2::aes_string(x = timeColumn, y = valueColumn, color = "group", group = participant_column)) +
          ggplot2::geom_point(alpha = 0.7) +
          ggplot2::geom_line(alpha = 0.3) +
          ggplot2::scale_color_manual(values = get_plot_palette(list(df = df))) +
          ggplot2::scale_fill_manual(values = get_plot_palette(list(df = df))) +
          my_theme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = x_label, y = xi, color = group_label, fill = group_label)
        if (!any(is.na(addSmooth))) {
          g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group, fill = group), se = TRUE)
        }
        return(g)
      }
    } else if (is(object, "ALASCA")) {
      df <- object$df_raw
      valueColumn <- as.character(object$formula)[2]
      if (any(participant_column == FALSE)) {
        if (any(object$participant_column == FALSE)) {
          object$log("You need to specify participant column", level = "ERROR")
          stop()
        } else {
          participant_column <- object$participant_column
        }
      }
      
      if (is.na(x_label)) x_label <- object$plot.x_label
      if (is.na(group_label)) group_label <- object$plot.group_label
      
      plotFunction <- function(df, timeColumn, valueColumn, participant_column, xi, addSmooth, x_label, group_label, my_theme) {
        g <- ggplot2::ggplot(subset(df, variable == xi), ggplot2::aes_string(x = timeColumn, y = valueColumn, color = "group", group = participant_column)) +
          ggplot2::geom_point(alpha = 0.7) +
          ggplot2::geom_line(alpha = 0.3) +
          ggplot2::scale_color_manual(values = get_plot_palette(object)) +
          ggplot2::scale_fill_manual(values = get_plot_palette(object)) +
          my_theme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = x_label, y = xi, color = group_label, fill = group_label)
        if (!any(is.na(addSmooth))) {
          g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group, fill = group), se = TRUE)
        }
        return(g)
      }
    } else {
      object$log("Wrong input object: must be a ALASCA model or a data frame", level = "ERROR")
      stop()
    }
    
    if (any(is.na(variables))) variables <- unique(df$variable)
    
    g <- lapply(variables, function(xi) {
      plotFunction(df, timeColumn, valueColumn, participant_column, xi, addSmooth, x_label = x_label, group_label = group_label, my_theme = my_theme)
    })
    names(g) <- variables
    if (is(object, "ALASCA")) {
      if (object$save) {
        for (i in seq_along(g)) {
          saveALASCAPlot(object = object, g = g[[i]], filetype = filetype, figsize = figsize, figunit = figunit, suffix = names(g)[i])
        }
      }
    }
    return(g)
  } else {
    if (any(is.na(variables))) {
      variables <- object$variablelist
    }
    df <- object$df_raw
    if (is.na(x_label)) {
      x_label <- object$plot.x_label
    }
    if (is.na(group_label)) {
      group_label <- object$plot.group_label
    }
    g <- ggplot2::ggplot(subset(df, variable %in% variables), ggplot2::aes(x = time, y = value, color = group, group = ID)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_line(alpha = 0.3) +
      ggplot2::scale_color_manual(values = get_plot_palette(object)) +
      ggplot2::scale_fill_manual(values = get_plot_palette(object)) +
      my_theme +
      ggplot2::facet_wrap(~variable, scales = "free_y") +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = x_label, y = plot.ylabel, color = group_label, fill = group_label)
    if (!any(is.na(addSmooth))) {
      g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group, fill = group), se = TRUE)
    }
    
    if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
    
    return(g)
  }
}

plot_prediction <- function() {
  
  if (is.null(self$variable)) {
    self$model$log(paste0("Selecting the ",round(self$n_limit/2)," variables with highest/lowest loading on `", self$model$effect_list$expr[[self$effect_i[[1]]]], "` (PC",self$component[[1]],"). Use `variable` to specify variables to plot"))
    variables_to_plot <- self$model$get_loadings(effect_i = self$effect_i[[1]], component = self$component[[1]], n_limit = round(self$n_limit/2))[[1]]
  } else if (self$variable == 0) {
    variables_to_plot <- self$model$get_loadings(effect_i = self$effect_i[[1]], component = self$component[[1]])
  }
  
  effect_terms <- self$model$effect_list$terms[[self$effect_i[[1]]]]
  data_to_plot <- self$model$get_predictions()
  data_to_plot <- data_to_plot[variable %in% variables_to_plot$covars]
  
  data_to_plot[, variable := factor(variable, levels = variables_to_plot[order(loading, decreasing = TRUE), covars])]
  colnames(data_to_plot)[colnames(data_to_plot) == effect_terms[[1]]] <- "x_data"
  
  if (length(self$model$effect_terms) == 1) {
    # Use some reference
    data_to_plot[, group_data := self$model$get_ref(self$model$get_plot_group)]
  } else  {
    if (self$make_group_column) {
      data_to_plot[, group_data := do.call(paste, c(.SD, sep = "-")), .SDcols = self$model$effect_terms[-1]]
    } else {
      data_to_plot[, group_data := get(self$model$get_plot_group)]
    }
  }
  
  
  if(self$model$validate) {
    g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "x_data",
                                                           y = "pred",
                                                           group = "group_data",
                                                           color = "group_data",
                                                           linetype = "group_data",
                                                           ymin = "low",
                                                           ymax = "high")) +
      ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = self$dodgewidth)) +
      ggplot2::geom_line(position = ggplot2::position_dodge(width = self$dodgewidth))
    if (self$ribbon) {
      g <- g + ggplot2::geom_ribbon(ggplot2::aes_string(fill = "group_data"),
                                    alpha = .1,
                                    position = ggplot2::position_dodge(width = self$dodgewidth), color = NA
      ) + ggplot2::scale_fill_manual(values = self$get_plot_palette()) + ggplot2::labs(fill = self$group_label)
    }
  } else {
    g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "x_data",
                                                           y = "pred",
                                                           group = "group_data",
                                                           color = "group_data",
                                                           linetype = "group_data")) +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  }
  g <- g + ggplot2::scale_color_manual(values = self$get_plot_palette()) +
    ggplot2::scale_linetype_manual(values = self$get_plot_linetypes()) +
    ggplot2::labs(color = self$group_label,
                  linetype = self$group_label,
                  x = self$x_label,
                  y = "Std. value") +
    ggplot2::facet_wrap(~variable, scales = "free_y") + 
    self$my_theme + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = self$x_angle, vjust = self$x_v_just, hjust = self$x_h_just))
}

#' Plot marginal means
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param inheritParams plot_validation
#' @return A ggplot2 object.
#'
#' @export
plot_mm <- function(...){
  plot_prediction(...)
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param inheritParams plot_validation
#' @return A ggplot2 object.
#'
#' @export
plot_val <- function(...){
  plot_validation(...)
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param inheritParams plot_covars
#' @return A ggplot2 objects\.
#'
#' @export
plot_covar <- function(...){
  plot_covars(...)
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param inheritParams plot_covars
#' @return A ggplot2 objects\.
#'
#' @export
plot_covariates <- function(...){
  plot_covars(...)
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param inheritParams plot_covars
#' @return A ggplot2 objects\.
#'
#' @export
plot_covariate <- function(...){
  plot_covars(...)
}

plot_effect <- function() {
  
  effect_i <- self$effect_i
  component <- self$component
  
  if (length(effect_i) > 1 || length(component) > 1) {
    g_list <- list()
    ii <- 1
    # 
    # There are several cases to consider:
    # - No variable group
    #   - n_col = 1: Share group legend from last plot
    #   - n_col = 2: Share group legend from last plot
    # - Defined variable group
    #   - n_col = 1: Share combined legend from last plot
    #   - n_col = 2: Split group legend and variable legend
    
    for (effect_k in rev(effect_i)) {
      for (component_k in rev(component)) {
        
        gs <- self$plot_effect_score(effect_i = effect_k, component = component_k)
        gl <- self$plot_effect_loading(effect_i = effect_k, component = component_k)
        
        if (ii == 1 && (self$n_col == 1 || is.null(self$loading_group_column))) {
          # The first iteration is, in fact, the last, so add legend here
          if (is.null(self$loading_group_column)) {
            if (self$n_col == 2) {
              gs_legend <- gs # use this later
              g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "none", align = "h")
            } else {
              g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h")
            }
            
          } else {
            g <- ggpubr::ggarrange(gs, gl, align = "h")
          }
        } else if (ii <= 2 && self$n_col == 2 && !is.null(self$loading_group_column)) {
          if (ii == 1) {
            # Should hold the variable legend
            g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h", legend.grob = ggpubr::get_legend(gl))
          } else {
            # Should hold the group legend
            g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h", legend.grob = ggpubr::get_legend(gs))
          }
        } else {
          # No legend in these plots
          g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "none", align = "h")
        }
        g_list[[ii]] <- g
        ii <- ii+1
      }
    }
    if (is.null(self$loading_group_column) && self$n_col == 2) {
      # Set common group legend
      do.call(ggpubr::ggarrange, list(plotlist = rev(g_list),
                                      ncol = self$n_col,
                                      nrow = length(g_list)/self$n_col,
                                      labels = self$labels,
                                      common.legend = TRUE,
                                      legend = "bottom",
                                      legend.grob = ggpubr::get_legend(gs_legend))
      )
    } else {
      # The legends are already in the subfigures
      do.call(ggpubr::ggarrange, list(plotlist = rev(g_list),
                                      ncol = self$n_col,
                                      nrow = length(g_list)/self$n_col,
                                      labels = self$labels)
      )
    }
  } else {
    # Return only a single effect/component
    gs <- self$plot_effect_score(effect_i = effect_i, component = component)
    gl <- self$plot_effect_loading(effect_i = effect_i, component = component)
    if (is.null(self$loading_group_column)) {
      ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h", labels = self$labels)
    } else {
      ggpubr::ggarrange(gs, gl, align = "h", labels = self$labels)
    }
  }
}

plot_effect_score <- function(effect_i = 1, component = 1) {
  
  effect_terms <- self$model$effect_list$terms[[effect_i]]
  data_to_plot <- self$model$get_scores(effect_i = effect_i, component = component)
  colnames(data_to_plot)[colnames(data_to_plot) == effect_terms[[1]]] <- "x_data"
  
  if (length(effect_terms) == 1 && self$model$method == "LMM") {
    # Use some reference
    data_to_plot[, group_data := self$model$get_ref(self$model$get_plot_group)]
  } else {
    if (self$make_group_column) {
      data_to_plot[, group_data := do.call(paste, c(.SD, sep = "-")), .SDcols = self$model$effect_terms[-1]]
    } else {
      if (self$model$get_plot_group %in% colnames(data_to_plot)) {
        data_to_plot[, group_data := get(self$model$get_plot_group)]
      } else {
        data_to_plot[, group_data := x_data]
      }
    }
  }
  
  if (self$validate) {
    # Validated model
    
    if(self$model$method == "LMM") {
      g <- ggplot2::ggplot(data_to_plot,
                           ggplot2::aes_string(x = "x_data",
                                               y = "score",
                                               group = "group_data",
                                               color = "group_data",
                                               linetype = "group_data",
                                               ymin = "low",
                                               ymax = "high")) +
        ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = self$dodgewidth)) +
        ggplot2::geom_line(position = ggplot2::position_dodge(width = self$dodgewidth))
    } else {
      g <- ggplot2::ggplot(data_to_plot,
                           ggplot2::aes_string(x = "x_data",
                                               y = "score",
                                               color = "group_data",
                                               ymin = "low",
                                               ymax = "high")) +
        ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = self$dodgewidth))
    }
      
    if (self$ribbon && self$model$method == "LMM") {
      g <- g + ggplot2::geom_ribbon(ggplot2::aes_string(fill = "group_data"),
                                    alpha = .1,
                                    position = ggplot2::position_dodge(width = self$dodgewidth), color = NA
      ) + ggplot2::scale_fill_manual(values = self$get_plot_palette()) + ggplot2::labs(fill = self$group_label)
    }
  } else {
    # No validation
    
    if (self$model$method == "LMM") {
      g <- ggplot2::ggplot(data_to_plot,
                           ggplot2::aes_string(x = "x_data",
                                               y = "score",
                                               linetype = "group_data",
                                               group = "group_data",
                                               color = "group_data")) +
        ggplot2::geom_point() +
        ggplot2::geom_line()
    } else {
      g <- ggplot2::ggplot(data_to_plot,
                           ggplot2::aes_string(x = "x_data",
                                               y = "score",
                                               color = "group_data")) +
        ggplot2::geom_point()
    }
  }
  g <- g +
    ggplot2::scale_color_manual(values = self$get_plot_palette()) +
    ggplot2::scale_linetype_manual(values = self$get_plot_linetypes()) +
    ggplot2::labs(color = self$group_label,
                  linetype = self$group_label,
                  x = self$x_label,
                  y = self$get_explained_label(effect_i = effect_i, component = component, type= "Score")) +
    self$my_theme + self$xflip(flip = FALSE)
  return(g)
}

plot_effect_loading <- function(effect_i = 1, component = 1) {
  if (self$n_limit > 0) {
    self$model$log(paste("Showing", self$n_limit*2, "of",length(self$model$get_levels("variable")),"variables. Adjust the number with `n_limit`"), level = "WARN")
  }
  data_to_plot <- self$model$get_loadings(effect_i = effect_i, component = component, n_limit = self$n_limit)[[1]]
  data_to_plot[, covars := factor(covars, levels = data_to_plot[order(loading), covars])]
  
  if(!is.null(self$loading_group_column)) {
    data_to_plot <- merge(data_to_plot, self$model$variable_labels, by = "covars")
    if (self$sort_by_loading_group) {
      data_to_plot[, covars := factor(covars, levels = covars[order(covargroup, loading, decreasing = TRUE)])]
    }
  }
  
  if (self$validate) {
    # Validated model
    
    if (is.null(self$loading_group_column)) {
      g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covars", y = "loading", ymin = "low", ymax = "high"))
    } else {
      g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covars", y = "loading", ymin = "low", ymax = "high", shape = "covargroup", color = "covargroup"))
    }
    
    g <- g + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_pointrange()
    
  } else {
    # Unvalidated model
    
    if (is.null(self$loading_group_column)) {
      g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covars", y = "loading"))
    } else {
      g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covars", y = "loading", shape = "covargroup", color = "covargroup"))
    }
    
    g <- g + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_point()
  }
  
  g <- g +
    ggplot2::labs(x = self$variable_label,
                  y = self$get_explained_label(effect_i = effect_i, component = component, type= "Loading")) +
    self$my_theme + self$xflip()
  
  if (!is.null(self$loading_group_column)) {
    g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
      ggplot2::labs(color = self$loading_group_label, shape = self$loading_group_label)
  }
  
  return(g)
}

plot_effect_validation <- function() {
  
  effect_i <- self$effect_i
  component <- self$component
  
  if (length(effect_i) > 1 || length(component) > 1) {
    g_list <- list()
    ii <- 1
    # 
    # There are several cases to consider:
    # - No variable group
    #   - n_col = 1: Share group legend from last plot
    #   - n_col = 2: Share group legend from last plot
    # - Defined variable group
    #   - n_col = 1: Share combined legend from last plot
    #   - n_col = 2: Split group legend and variable legend
    
    for (effect_k in rev(effect_i)) {
      for (component_k in rev(component)) {
        
        gs <- self$plot_effect_validation_score(effect_i = effect_k, component = component_k)
        gl <- self$plot_effect_validation_loading(effect_i = effect_k, component = component_k)
        
        if (ii == 1 && (self$n_col == 1 || is.null(self$loading_group_column))) {
          # The first iteration is, in fact, the last, so add legend here
          if (is.null(self$loading_group_column)) {
            if (self$n_col == 2) {
              gs_legend <- gs # use this later
              g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "none", align = "h")
            } else {
              g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h")
            }
            
          } else {
            g <- ggpubr::ggarrange(gs, gl, align = "h")
          }
        } else if (ii <= 2 && self$n_col == 2 && !is.null(self$loading_group_column)) {
          if (ii == 1) {
            # Should hold the variable legend
            g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h", legend.grob = ggpubr::get_legend(gl))
          } else {
            # Should hold the group legend
            g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h", legend.grob = ggpubr::get_legend(gs))
          }
        } else {
          # No legend in these plots
          g <- ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "none", align = "h")
        }
        g_list[[ii]] <- g
        ii <- ii+1
      }
    }
    if (is.null(self$loading_group_column) && self$n_col == 2) {
      # Set common group legend
      do.call(ggpubr::ggarrange, list(plotlist = rev(g_list),
                                      ncol = self$n_col,
                                      nrow = length(g_list)/self$n_col,
                                      labels = self$labels,
                                      common.legend = TRUE,
                                      legend = "bottom",
                                      legend.grob = ggpubr::get_legend(gs_legend))
      )
    } else {
      # The legends are already in the subfigures
      do.call(ggpubr::ggarrange, list(plotlist = rev(g_list),
                                      ncol = self$n_col,
                                      nrow = length(g_list)/self$n_col,
                                      labels = self$labels)
      )
    }
  } else {
    # Return only a single effect/component
    gs <- self$plot_effect_validation_score(effect_i = effect_i, component = component)
    gl <- self$plot_effect_validation_loading(effect_i = effect_i, component = component)
    if (is.null(self$loading_group_column)) {
      ggpubr::ggarrange(gs, gl, common.legend = TRUE, legend = "bottom", align = "h", labels = self$labels)
    } else {
      ggpubr::ggarrange(gs, gl, align = "h", labels = self$labels)
    }
  }
}

plot_effect_validation_score <- function(effect_i = 1, component = 1) {
  
  if (!self$model$validate) {
    self$model$log("The model has not been validated", level = "ERROR")
    stop()
  }
  
  effect_terms <- self$model$effect_list$terms[[effect_i]]
  data_to_plot <- self$model$get_scores(effect_i = effect_i, component = component)
  data_to_plot$model <- 0
  data_to_plot$high <- NULL
  data_to_plot$low <- NULL
  data_to_add <- self$model$get_validation_scores(effect_i = effect_i)[PC == component]
  data_to_plot <- rbind(data_to_plot, data_to_add)
  colnames(data_to_plot)[colnames(data_to_plot) == effect_terms[[1]]] <- "x_data"
  
  if (length(effect_terms) == 1 && self$model$method == "LMM") {
    # Use some reference
    data_to_plot[, group_data := self$model$get_ref(self$model$get_plot_group)]
  } else {
    if (self$make_group_column) {
      data_to_plot[, group_data := do.call(paste, c(.SD, sep = "-")), .SDcols = self$model$effect_terms[-1]]
    } else {
      if (self$model$get_plot_group %in% colnames(data_to_plot)) {
        data_to_plot[, group_data := get(self$model$get_plot_group)]
      } else {
        data_to_plot[, group_data := x_data]
      }
    }
  }
  
  data_to_plot$grouping <- paste(data_to_plot$model, "-", data_to_plot$group_data)
  
  # Validated model
  g <- ggplot2::ggplot(data_to_plot[model != 0],
                       ggplot2::aes_string(x = "x_data",
                                           y = "score",
                                           group = "grouping",
                                           color = "group_data",
                                           linetype = "group_data")) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_point(data = data_to_plot[model == 0], color = "black") +
    ggplot2::geom_line(data = data_to_plot[model == 0], color = "black") +
    ggplot2::scale_alpha(range = c(0.1, 1), guide = "none") +
    ggplot2::scale_color_manual(values = self$get_plot_palette()) +
    ggplot2::scale_linetype_manual(values = self$get_plot_linetypes()) +
    ggplot2::labs(color = self$group_label,
                  linetype = self$group_label,
                  x = self$x_label,
                  y = self$get_explained_label(effect_i = effect_i, component = component, type= "Score")) +
    self$my_theme + self$xflip(flip = FALSE)
  return(g)
}

plot_effect_validation_loading <- function(effect_i = 1, component = 1) {
  
  if (!self$model$validate) {
    self$model$log("The model has not been validated", level = "ERROR")
    stop()
  }
  
  if (self$n_limit > 0) {
    self$model$log(paste("Showing", self$n_limit*2, "of",length(self$model$get_levels("variable")),"variables. Adjust the number with `n_limit`"), level = "WARN")
  }
  data_to_plot <- self$model$get_loadings(effect_i = effect_i, component = component, n_limit = self$n_limit)[[1]]
  data_to_plot$model <- 0
  data_to_add <- self$model$get_validation_loadings(effect_i = effect_i)[PC == component]
  data_to_add <- data_to_add[covars %in% unique(data_to_plot$covars)]
  data_to_plot <- rbind(data_to_plot, data_to_add, fill = TRUE)
  data_to_plot[, covars := factor(covars, levels = data_to_plot[model == 0][order(loading), covars])]
  
  if(!is.null(self$loading_group_column)) {
    data_to_plot <- merge(data_to_plot, self$model$variable_labels, by = "covars")
    if (self$sort_by_loading_group) {
      data_to_plot[, covars := factor(covars, levels = covars[order(covargroup, loading, decreasing = TRUE)])]
    }
  }
  data_to_plot[, nval := as.numeric(covars)]
  
    if (is.null(self$loading_group_column)) {
      g <- ggplot2::ggplot(data_to_plot[model != 0], ggplot2::aes_string(x = "covars", y = "loading", group = "model"))
    } else {
      g <- ggplot2::ggplot(data_to_plot[model != 0], ggplot2::aes_string(x = "covars", y = "loading", group = "model", shape = "covargroup", color = "covargroup"))
    }
    
    g <- g + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_point(alpha = 0.3) +
      ggplot2::geom_pointrange(data = data_to_plot[model == 0], ggplot2::aes(ymin = low, ymax = high), alpha = 1) +
      ggplot2::geom_errorbarh(data = data_to_plot[model == 0], ggplot2::aes(xmin = nval-0.6, xmax = nval+0.6), alpha = 1) +
      ggplot2::labs(x = self$variable_label,
                    y = self$get_explained_label(effect_i = effect_i, component = component, type= "Loading")) +
      ggplot2::scale_alpha(range = c(0, 1)) +
      self$my_theme + self$xflip()
  
  if (!is.null(self$loading_group_column)) {
    g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
      ggplot2::labs(color = self$loading_group_label, shape = self$loading_group_label)
  }
  
  return(g)
}

#' Prints a marco object
#' 
#' @param object A marco object.
#' @param variable A marco object.
#' 
#' @return A data table with residuals
#' 
#' @export
residuals.AlascaModel <- function(object, variable = NULL) {
  object$get_residuals(variable = variable)
}

#' Returns model loadings
#' 
#' @param object An ALASCA object.
#' @param effect The number of the effect(s) of interest
#' @param component The principal component(s) of interest
#' @param n_limit Limit the number of variables returned
#' 
#' @return A data table with loadings
#' 
#' @export
get_loadings.AlascaModel <- function(object, effect = 1, component = 1, n_limit = 0) {
  object$log(paste("Loadings for effect number", effect, "and principal component(s)", component))
  object$get_loadings(effect_i = effect, component = component, n_limit = n_limit)
}

#' Returns model scores
#' 
#' @inheritParams get_loadings.AlascaModel
#' 
#' @return A data table with loadings
#' 
#' @export
get_loadings <- function(object, ...) {
  UseMethod("get_loadings")
}

#' Returns model scores
#' 
#' @param object An ALASCA object.
#' @param effect The number of the effect of interest
#' @param component The principal component(s) of interest
#' 
#' @return A data table with scores
#' 
#' @export
get_scores.AlascaModel <- function(object, effect = 1, component = 1) {
  object$log(paste("Scores for effect number", effect, "and principal component(s)", component))
  object$get_scores(effect_i = effect, component = component)
}

#' Returns model scores
#' 
#' @inheritParams get_scores.AlascaModel
#' 
#' @return A data table with scores
#' 
#' @export
get_scores <- function(object, ...) {
  UseMethod("get_scores")
}

#' Get coefficients for covariates
#'
#'
#' @param object An ALASCA object
#' @param n_limit Return only the `n_limit` largest and smallest coefficients
#' @return A data table with regression coefficients
#'
#' @export
get_covars.AlascaModel <- function(object, n_limit = 0) {
  object$get_covars(n_limit = n_limit)
}

#' Get coefficients for covariates
#' 
#' @inheritParams get_covars.AlascaModel
#' 
#' @return A data table with regression coefficients
#' 
#' @export
get_covars <- function(object, ...) {
  UseMethod("get_covars")
}

#' Returns model predictions
#' 
#' @param object An ALASCA object
#' 
#' @return A data table with model predictions
#' 
#' @export
get_predictions.AlascaModel <- function(object) {
  object$get_predictions()
}

#' Get coefficients for covariates
#' 
#' @inheritParams get_predictions.AlascaModel
#' 
#' @return A data table with regression coefficients
#' 
#' @export
get_predictions <- function(object) {
  UseMethod("get_predictions")
}

plot_2D <- function() {
  
  if (!self$model$validate) {
    self$model$log("Please validate the model first", level = "ERROR")
    stop()
  }
  if (length(self$component) != 2) {
    self$model$log("Please provide two validated components, e.g. `component = c(1,2)`", level = "ERROR")
    stop()
  }
  
  gs <- self$plot_2D_score()
  gl1 <- self$plot_2D_loading_1()
  gl2 <- self$plot_2D_loading_2()
  tmp <- self$component
  self$component <- 10
  gscree <- self$plot_scree()
  g <- ggpubr::ggarrange(gs, gl2, gl1, gscree,
                    heights = c(3,2),
                    widths = c(3,2),
                    common.legend = TRUE,
                    legend = "top",
                    labels = self$labels)
  self$component <- tmp
  return(g)
}

plot_2D_score <- function() {
  
  effect_terms <- self$model$effect_list$terms[[self$effect_i]]
  data_to_plot <- self$model$get_scores(effect_i = self$effect_i, component = c(1,2))[, -c("high", "low")]
  data_to_plot$model <- 0
  if (self$model$validate) {
    data_to_add <- self$model$get_validation_scores(effect_i = self$effect_i)[PC %in% c(1,2)]
    data_to_plot <- rbind(data_to_plot, data_to_add)
  }
  data_to_plot$PC <- paste0("PC_", data_to_plot$PC)
  colnames(data_to_plot)[colnames(data_to_plot) == effect_terms[[1]]] <- "x_data"
  
  if (length(effect_terms) == 1 && self$model$method == "LMM") {
    # Use some reference
    data_to_plot[, group_data := self$model$get_ref(self$model$get_plot_group)]
  } else {
    if (self$make_group_column) {
      data_to_plot[, group_data := do.call(paste, c(.SD, sep = "-")), .SDcols = self$model$effect_terms[-1]]
    } else {
      if (self$model$get_plot_group %in% colnames(data_to_plot)) {
        data_to_plot[, group_data := get(self$model$get_plot_group)]
      } else {
        data_to_plot[, group_data := x_data]
      }
    }
  }
  
  data_to_plot <- dcast(data = data_to_plot, ...~PC, value.var = "score")
  
  g <- ggplot2::ggplot(data_to_plot,
                       ggplot2::aes_string(x = "PC_1",
                                           y = "PC_2",
                                           shape = "x_data",
                                           color = "group_data")) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(group = paste(group_data, model)), alpha = 0.3) +
    ggplot2::stat_ellipse(ggplot2::aes(linetype = x_data)) +
    ggplot2::scale_color_manual(values = self$get_plot_palette()) +
    ggplot2::scale_alpha(range = c(0, 1)) +
    ggplot2::labs(color = self$group_label,
                  linetype = self$x_label,
                  shape = self$x_label,
                  x = self$get_explained_label(effect_i = self$effect_i, component = self$component[[1]], type= "Score"),
                  y = self$get_explained_label(effect_i = self$effect_i, component = self$component[[2]], type= "Score")) +
    self$my_theme
}

plot_2D_loading_1 <- function() {
  
  data_to_plot <- self$model$get_loadings(effect_i = self$effect_i, component = self$component[[1]], n_limit = self$n_limit)[[1]]
  
  data_to_plot[, covars := factor(covars, levels = data_to_plot[order(loading), covars])]
  data_to_plot$facet <- ifelse(data_to_plot$loading > median(data_to_plot$loading), "Upper loadings \u2192", "\u2190 Lower loadings")
  
  ggpubr::ggarrange(
    ggplot2::ggplot(data_to_plot[facet == "\u2190 Lower loadings"][, covars := factor(covars, levels = covars[order(loading, decreasing = TRUE)])], 
                    ggplot2::aes(loading, covars, xmin = low, xmax = high)) + 
      ggplot2::geom_pointrange() +
      ggplot2::facet_grid(~facet) +
      ggplot2::scale_y_discrete(position = "left") +
      ggplot2::scale_x_continuous(limits = c(NA, max(0, data_to_plot[facet == "\u2190 Lower loadings"]$loading))) +
      self$my_theme +
      ggplot2::labs(x = "\u2190 Loading PC1") +
      ggplot2::theme(axis.title.y = element_blank()),
    ggplot2::ggplot(data_to_plot[facet == "Upper loadings \u2192"], ggplot2::aes(loading, covars, xmin = low, xmax = high)) +
      ggplot2::geom_pointrange() +
      ggplot2::facet_grid(~facet) +
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::scale_x_continuous(limits = c(min(0, data_to_plot[facet == "Upper loadings \u2192"]$loading), NA)) +
      self$my_theme +
      ggplot2::labs(x = "Loading PC1 \u2192") +
      ggplot2::theme(axis.title.y = element_blank())
  )
}

plot_2D_loading_2 <- function() {
  
  data_to_plot <- self$model$get_loadings(effect_i = self$effect_i, component = self$component[[2]], n_limit = self$n_limit)[[1]]
  
  data_to_plot[, covars := factor(covars, levels = data_to_plot[order(loading), covars])]
  data_to_plot$facet <- ifelse(data_to_plot$loading < median(data_to_plot$loading), "\u2190 Lower loadings", "Upper loadings \u2192")
  
  ggpubr::ggarrange(
    ggplot2::ggplot(data_to_plot[facet == "Upper loadings \u2192"], ggplot2::aes(loading, covars, xmin = low, xmax = high)) + 
      ggplot2::geom_pointrange() +
      ggplot2::facet_grid(rows = vars(facet), switch = "y") +
      ggplot2::scale_x_continuous(limits = c(min(0, data_to_plot[facet == "Upper loadings \u2192"]$loading), NA), position = "top") +
      self$my_theme +
      ggplot2::labs(x = "Loading PC2 \u2192") +
      ggplot2::theme(axis.title.y = element_blank()),
    ggplot(data_to_plot[facet == "\u2190 Lower loadings"], ggplot2::aes(loading, covars, xmin = low, xmax = high)) +
      ggplot2::geom_pointrange() +
      ggplot2::facet_grid(rows = vars(facet), switch = "y") +
      ggplot2::scale_x_continuous(limits = c(NA, max(0, data_to_plot[facet == "\u2190 Lower loadings"]$loading))) +
      self$my_theme +
      ggplot2::labs(x = "\u2190 Loading PC2") +
      ggplot2::theme(axis.title.y = element_blank()),
    ncol = 1, align = "v")
}

plot_histogram_score <- function() {
  
  if (!self$model$validate) {
    self$model$log("The model has not been validated", level = "ERROR")
    stop()
  }
  
  effect_terms <- self$model$effect_list$terms[[self$effect_i]]
  data_to_plot <- self$model$get_scores(effect_i = self$effect_i, component = self$component)
  data_to_plot$model <- 0
  data_to_add <- self$model$get_validation_scores(effect_i = self$effect_i)[PC == self$component]
  data_to_plot <- rbind(data_to_plot, data_to_add, fill = TRUE)
  colnames(data_to_plot)[colnames(data_to_plot) == effect_terms[[1]]] <- "x_data"
  data_to_plot[, x_data := paste0(self$x_label, ": ", x_data)]
  
  if (length(effect_terms) == 1 && self$model$method == "LMM") {
    # Use some reference
    data_to_plot[, group_data := self$model$get_ref(self$model$get_plot_group)]
  } else {
    if (self$make_group_column) {
      data_to_plot[, group_data := do.call(paste, c(.SD, sep = "-")), .SDcols = self$model$effect_terms[-1]]
    } else {
      if (self$model$get_plot_group %in% colnames(data_to_plot)) {
        data_to_plot[, group_data := get(self$model$get_plot_group)]
      } else {
        data_to_plot[, group_data := x_data]
      }
    }
  }
  
  # Lower alpha to the validation runs
  data_to_plot$alpha <- ifelse(data_to_plot$model > 0, 0.7, 1)
  
  data_to_plot$grouping <- paste(data_to_plot$model, "-", data_to_plot$group_data)
  
  g <- ggplot2::ggplot(data_to_plot, ggplot2::aes(score, fill = group_data)) +
    ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = self$n_bins) +
    ggplot2::geom_vline(data = data_to_plot[model == 0], ggplot2::aes(xintercept = score, color = group_data, linetype = group_data)) +
    ggplot2::scale_fill_manual(values = self$get_plot_palette()) +
    ggplot2::scale_color_manual(values = self$get_plot_palette()) +
    ggplot2::scale_linetype_manual(values = self$get_plot_linetypes()) +
    ggplot2::facet_wrap(~x_data, ncol = self$n_col) +
    ggplot2::labs(color = self$group_label,
                  fill = self$group_label,
                  linetype = self$group_label,
                  x = self$get_explained_label(effect_i = self$effect_i, component = self$component, type= "Score"),
                  y = "Count") +
    self$my_theme
  
} 

plot_histogram_loading <- function() {
  
  if (!self$model$validate) {
    self$model$log("The model has not been validated", level = "ERROR")
    stop()
  }
  
  if (self$n_limit > 0) {
    self$model$log(paste("Showing", self$n_limit*2, "of",length(self$model$get_levels("variable")),"variables. Adjust the number with `n_limit`"), level = "WARN")
  }
  data_to_plot <- self$model$get_loadings(effect_i = self$effect_i, component = self$component, n_limit = self$n_limit)[[1]]
  data_to_plot$model <- 0
  data_to_add <- self$model$get_validation_loadings(effect_i = self$effect_i)[PC == self$component]
  data_to_add <- data_to_add[covars %in% unique(data_to_plot$covars)]
  data_to_plot <- rbind(data_to_plot, data_to_add, fill = TRUE)
  data_to_plot[, covars := factor(covars, levels = data_to_plot[model == 0][order(loading, decreasing = TRUE), covars])]
  
  if(!is.null(self$loading_group_column)) {
    data_to_plot <- merge(data_to_plot, self$model$variable_labels, by = "covars")
    if (self$sort_by_loading_group) {
      data_to_plot[, covars := factor(covars, levels = covars[order(covargroup, loading, decreasing = TRUE)])]
    }
  }
  data_to_plot[, nval := as.numeric(covars)]
  
  ggplot2::ggplot(data_to_plot, ggplot2::aes(loading)) +
    ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = self$n_bins) +
    ggplot2::geom_vline(data = data_to_plot[model == 0], ggplot2::aes(xintercept = loading)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~covars) +
    ggplot2::labs(x = self$get_explained_label(effect_i = self$effect_i, component = self$component, type= "Loading"), y = "Count") +
    self$my_theme
  
}

plot_histogram <- function() {
  gl <- self$plot_histogram_loading()
  gs <- self$plot_histogram_score()
  ggpubr::ggarrange(gs, gl, widths = c(1,2), labels = self$labels, common.legend = TRUE, legend = "bottom")
}

