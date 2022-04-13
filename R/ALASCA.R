#' Get an ALASCA object
#'
#' `ALASCA` initializes an ALASCA model and returns an ALASCA object
#'
#' This function builds your ALASCA model. It needs a data frame containing at least a column identifying participants, a column called `time` contining time information, a column `group` containing group information, a column `variable` containing variable names, and a value column. In addition you need to specify the model you want, and whether you want to separate group and time effects (defaults to `TRUE`).
#'
#' @param df Data frame to be analyzed
#' @param formula Regression model
#' @param separate_time_and_group Logical: should time and group effect be separated?
#' @param p_adjust_method Method for correcting p values for multiple testing, see p.adjust.methods
#' @param validate Logical. If `TRUE`, give estimates for robustness
#' @param participant_column String. Name of the column containing participant identification
#' @param minimize_object Logical. If `TRUE`, remove unnecessary clutter, optimize for validation
#' @param scale_function Either a custom function or string to define scaling function: `sdall`, `sdref`, `sdt1`, `sdreft1`
#' @param scale_function.center Boolean. Mean centering as part of scaling
#' @param equal_baseline Set to `TRUE` (default) to remove interaction between group and first time point
#' @param sum_coding Set to `TRUE` to use sum coding instead of contrast coding for group (defaults to `FALSE`)
#' @param plot.x_label Defaults to "Time"
#' @param plot.group_label Defaults to "Group"
#' @param plot.figsize A vector containing `c(width,height,dpi)` (default: `c(180, 120, 300)`)
#' @param plot.figunit Defaults to "mm"
#' @param plot.filetype Which filetype you want to save the figure to (default: `png`)
#' @param plot.palette List of colors, named by group
#' @param keep_terms Additional terms to keep in the model matrix
#' @param stratification_vector Vector of same length as `df` that specifies stratification groups during validate_ Defaults to `NA`, where the group column is used.
#' @param validate_regression Whether to validate regression predictions or not (only if `validate` is `TRUE`)
#' @param do_debug Print what happens (default: `FALSE`)
#' @param save Save models and plots automatically (default: `FALSE`)
#' @param filename File name to save model and plots (when `save = TRUE`)
#' @param method Defaults to `NA` where method is either LM or LMM, depending on whether your formula contains a random effect or not
#' @param use_Rfast Boolean. Defaults to `TRUE`
#' @param n_validation_folds Partitions when validating
#' @param n_validation_runs number of validation runs
#' @param validation_method among  `loo` (leave-one-out, default)
#' @param validation_object Don't worry about me:)
#' @param validation_participants Don't worry about me:)
#' @return An ALASCA object
#'
#' @examples
#' load("PE.Rdata")
#' model <- ALASCA(df = df, formula = value ~ time * group + (1 | ID))
#' @export
ALASCA <- function(df,
                   formula,
                   ...) {

  object <- AlascaModel$new(df, formula, ...)
  
  # Validate the model ----
  if (object$validate) {
    log4r::debug(object$log, "Starting validation")
    object$do_validate()
    log4r::debug(object$log, "Completing validation")
  } else {
    if (object$do_debug) currentTs <- Sys.time()
    object$clean_alasca()
    if (object$do_debug) cat("* clean_alasca:", Sys.time() - currentTs, "s\n")
  }
  
  object$run_time <- Sys.time() - object$init_time

  # Save the model
  if (object$save & !object$minimize_object) {
    object$saveALASCA()
  }

  if (object$save_to_disk & !object$minimize_object) {
    DBI::dbDisconnect(object$db.con)
    DBI::dbUnloadDriver(object$db.driver)
  }

  return(object)
}

#' Run RMASCA
#'
#' Same as calling ALASCA with `method = "LMM"`
#'
#' @inheritParams ALASCA
#' @return An ALASCA object
#' @export
RMASCA <- function(...) {
  object <- ALASCA(..., method = "LMM")
  return(object)
}

#' Print ALASCA version
#'
#' @return String
#' @export
print_version <- function(object = FALSE, get = NA, print = TRUE) {
  ALASCA.version <- "0.0.0.92"
  ALASCA.version.date <- "2022-03-06"
  if (is.list(object)) {
    ALASCA.version <- object$ALASCA.version
    ALASCA.version.date <- object$ALASCA.version.date
  }
  if (is.na(get)) {
    if (print) {
      cat("\n\n====== ALASCA ======\n\n")
      cat(paste0(ALASCA.version, " (", ALASCA.version.date, ")\n\n"))
    }
  } else {
    if (get == "date") {
      return(ALASCA.version.date)
    } else {
      return(ALASCA.version)
    }
  }
}

#' Run SMASCA
#'
#' Same as calling ALASCA with `method = "LM"`
#'
#' @inheritParams ALASCA
#' @return An ALASCA object
#' @export
SMASCA <- function(...) {
  object <- ALASCA(..., method = "LM")
  return(object)
}

#' Sanitize an ALASCA object
#'
#' This function checks that the input to an ALASCA object is as expected
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
sanitize_object <- function() {
  if (!self$minimize_object) {
    
    log4r::debug(self$log, "Starting sanitation")
    self$get_info_from_formula()
    self$rename_columns_to_standard()
    self$wide_to_long()
    
    # Keep variable labels
    if (!is.na(self$plot.loadinggroupcolumn)) {
      self$variable_labels <- unique(self$df[, .SD, .SDcols = c("variable", self$plot.loadinggroupcolumn)])
      colnames(self$variable_labels) <- c("covars", "covargroup")
    }
    
    # Remove surplus data for efficiency
    self$df <- self$df[, .SD, .SDcols = c(self$all_formula_terms, "variable", "value")]
    
    self$check_that_columns_are_valid()
    
    log4r::info(self$log, "Making factors for time, group and variable")
    self$df[, time := factor(time, levels = self$timelist), ]
    self$df[, group := factor(group, levels = self$grouplist), ]
    self$df[, variable := factor(variable, levels = self$variablelist), ]
    self$stratification_vector <- self$df[, get(self$stratification_column)]
    
    # List of levels
    self$variablelist <- unique(self$df$variable)
    self$timelist <- levels(self$df$time)
    self$grouplist <- levels(self$df$group)
    
    self$adjust_design_matrix()
    
    self$get_scaling_function()
    self$get_pca_function()

    if (self$save_to_disk) {
      self$db.driver <- RSQLite::dbDriver("SQLite")
      self$db.filename <- get_filename(self, prefix = "validation/", filetype = "db")
      self$db.con <- DBI::dbConnect(self$db.driver, dbname = self$db.filename)
    }
    
    log4r::info(self$log, "Checking for missing information")
    self$find_missing_predictor_variables()
    self$find_missing_response_variables()
    
    # Use sum coding?
    if (self$sum_coding) {
      log4r::info(self$log, "Use sum coding")
      contrasts(self$df$group) <- contr.sum(length(unique(self$df$group)))
    }
  }

  # Keep a copy of unscaled data
  self$df_raw <- AlascaDataset$new(df = self$df)
  self$my_df_rows <- self$df_raw$rows_to_serve
  
  # Scale data
  self$df <- self$scale_function(self$df)
  
  log4r::debug(self$log, "Completed sanitation")
  #invisible(self)
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
rename_columns_to_standard <- function() {
  
  if (self$valCol != "value") {
    if ("value" %in% colnames(self$df)) {
      log4r::error(self$log, "Sorry, the value column is reserved by ALASCA; please give it another name or change `valCol`")
      stop()
    }
    log4r::warn(self$log, paste0("Changing ", self$valCol, " to `value`."))
    self$df[, value := get(self$valCol)]
    self$formula <- formula(paste(
      "value ~",
      as.character(self$formula)[3]
    ))
    self$get_info_from_formula()
  }
  
  if (self$x_column == "time" && ! "time" %in% self$all_formula_terms) {
    self$x_column <- self$formula_terms[1]
    log4r::info(self$log, paste("Will use", self$x_column, "for abscissa"))
  }
  
  if (self$x_column != "time") {
    if ("time" %in% colnames(self$df)) {
      self$time <- NULL
      log4r::warn(self$log, "Overwriting the `time` column")
    }
    self$df[, time := factor(get(self$x_column))]
    self$timelist <- levels(self$df$time)
    self$replace_term_in_formula(old_term = self$x_column, new_term = "time")
  }
  
  
  if (self$method == "LM" && !"group" %in% colnames(self$df)) {
    self$grouplist <- self$timelist
    self$df[, group := time]
  } else if (!"group" %in% colnames(self$df)) {
    self$df[, group := factor("NA")]
  } 
  
  self$variablelist <- unique(self$df$variable)
  self$timelist <- levels(self$df$time)
  self$grouplist <- levels(self$df$group)
  
  log4r::info(self$log, paste0("Using `",self$stratification_column,"` for stratification"))

  #invisible(self)
}

#' Check if response variables are missing
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
find_missing_response_variables <- function() {
  if(self$df[, uniqueN(variable), by = .(ID, time)][, uniqueN(V1)] > 1) {
    if (self$ignore_missing) {
      log4r::warn(self$log, "Response variables missing for some samples! Continue with caution!")
    } else {
      log4r::error(self$log, "Response variables missing for some samples! To ignore this, use `ignore_missing = TRUE`")
      stop()
    }
  }
}

#' Check if predictor variables are missing
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
find_missing_predictor_variables <- function() {
  if (any(is.na(self$df))) {
    if (self$ignore_missing_covars) {
      log4r::warn(self$log, "Predictor variables missing for some samples! Continue with caution!")
    } else {
      log4r::error(self$log, "Predictor variables missing for some samples! To ignore this, use `ignore_missing_covars = TRUE`")
      stop()
    }
  }
}

#' Check if response variables are missing
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
replace_term_in_formula <- function(old_term, new_term) {
  old_formula <- as.character(self$formula)[3]
  new_formula <- gsub(old_term, new_term, old_formula)
  log4r::info(self$log, paste0("New formula: ", new_formula))
  self$formula <- formula(paste(
    "value ~",
    new_formula
  ))
  self$get_info_from_formula()
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
check_that_columns_are_valid <- function() {
  if (any(!self$all_formula_terms %in% colnames(self$df))) {
    log4r::error(self$log, paste0("Column(s) missing:\n", paste0(self$all_formula_terms[!self$all_formula_terms %in% colnames(self$df)], collapse = "\n* "), "\nYou may want to use `keep_columns`"))
    stop()
  }
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
adjust_design_matrix <- function() {
  if (self$method %in% c("LMM")) {
    if (self$participant_column != "ID") {
      
      # The user has specified a column
      self$df[, ID := get(self$participant_column)]
      self <- replace_term_in_formula(self, old_term = self$participant_column, new_term = "ID")
      
    } else if (any(grepl("\\|ID", self$formula_terms))) {
      
      # Use ID for participants!
      
    } else if (sum(grepl("\\|", self$formula_terms)) > 1) {
      log4r::error(self$log, "Multiple random effects, couldn't determine participant-id. Please specify `participant_column`")
      stop()
    } else {
        
        # Try to find ID column from formula
        tmp <- self$formula_terms[grepl("\\|", self$formula_terms)]
        tmp <- gsub(" ", "", tmp)
        tmp <- strsplit(tmp, "\\|")
        self$participant_column <- tmp[[1]][2]
        self$df[, ID := get(self$participant_column)]
        self$replace_term_in_formula(old_term = self$participant_column, new_term = "ID")
    }
    
    if (self$use_Rfast) {
      # Using Rfast
      fixed_terms <- self$formula_terms[!grepl("\\|", self$formula_terms)]
      self$new_formula <- formula(paste("value ~ ", paste(fixed_terms, collapse = "+")))
    } else {
      # Using lme4
      rterms <- self$formula_terms[grepl("\\|", self$formula_terms)]
      rterms <- paste0("(", rterms, ")")
      self$new_formula <- formula(paste("value ~ modmat+", paste(rterms, collapse = "+")))
    }
  }else if (self$method %in% c("LM")) {
    self$new_formula <- value ~ modmat
  } else {
    log4r::error(self$log, "Sorry, an error occurred! Please check your model")
    stop()
  }
  #invisible(self)
}

#' Get information from formula
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_info_from_formula <- function() {
  
  # Response variable
  self$valCol <- as.character(self$formula)[2]
  
  # Terms in the regression model
  self$formula_terms <- gsub(" ", "", colnames(attr(terms.formula(self$formula), "factors")))
  
  # Get a list of all predictors
  self$all_formula_terms <- unlist(strsplit(self$formula_terms, split = "\\:|\\+|\\||\\*"))
  self$all_formula_terms <- gsub(" ", "", self$all_formula_terms)
  self$all_formula_terms <- unique(self$all_formula_terms[!self$all_formula_terms %in% c("1", "")])
  if (self$stratification_column == "group" && !"group" %in% colnames(self$df)) {
    self$stratification_column <- self$all_formula_terms[1]
    log4r::warn(self$log, paste("The `",self$all_formula_terms[1],"` column is used for stratification"))
  }
  if (!"group" %in% self$all_formula_terms) {
    self$all_formula_terms <- c(self$all_formula_terms, "group")
    if (self$stratification_column == "group") {
      log4r::warn(self$log, "The `group` column is used for stratification")
    }
  }
  self$all_formula_terms <- c(self$all_formula_terms, self$participant_column, self$stratification_column, self$keep_columns)
  self$all_formula_terms <- gsub(" ", "", self$all_formula_terms)
  self$all_formula_terms <- unique(self$all_formula_terms[!self$all_formula_terms %in% c("1", "")])
  
  ## We need to keep original IDs to have a unique identifier later on
  if (self$validation_method == "bootstrap") {
    self$all_formula_terms <- unique(c(self$all_formula_terms, "originalIDbeforeBootstrap"))
    self$df[, originalIDbeforeBootstrap := -1]
    self$all_formula_terms <- unique(c(self$all_formula_terms, "uniqueIDforBootstrap"))
    self$df[, uniqueIDforBootstrap := -1]
  }
  
  # Check what terms that is present in formula
  self$hasGroupTerm <- ifelse(any(self$formula_terms == "group"), TRUE, FALSE)
  self$hasInteractionTerm <- ifelse(any(self$formula_terms == "group:time" | self$formula_terms == "time:group"), TRUE, FALSE)
  self$covars <- self$formula_terms[!(self$formula_terms %in% c("time", "group", "group:time", "time:group", self$keep_terms))]
  if (self$do_debug) {
    cat(".... Group term in formula? ", self$hasGroupTerm, "\n")
    cat(".... Interaction term in formula? ", self$hasInteractionTerm, "\n")
    cat(".... Identified the following covariates in addition to time and group: ", self$covars, "\n")
  }
  
  self$LM_or_LMM()
  
  #invisible(self)
}

#' Determine whether to use LM or LMM
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
LM_or_LMM <- function() {
  if (is.na(self$method)) {
    # Find which default method to use
    if (any(grepl("\\|", self$formula_terms))) {
      self$method <- "LMM"
      log4r::info(self$log, "Will use linear mixed models!")
      if (sum(grepl("\\|", self$formula_terms)) > 1 && self$use_Rfast) {
        log4r::error(self$log, "Cannot use Rfast with multiple random effects. Use lme4 with `use_Rfast = FALSE` instead!")
        stop()
      }
      if (!any(grepl("1\\|ID", self$formula_terms)) && self$use_Rfast) {
        log4r::error(self$log, "Rfast only supports a single random intercept. Use lme4 with `use_Rfast = FALSE` instead!")
        stop()
      }
    } else {
      self$method <- "LM"
      log4r::info(self$log, "Will use linear models!")
    }
  } else {
    # The user has specified a method to use
    if (self$method == "LMM") {
      if (!any(grepl("\\|", self$formula_terms))) {
        log4r::error(self$log, "The model must contain at least one random effect. Are you sure you wanted linear mixed models?")
        stop()
      }
    } else if (self$method == "LM") {
      if (any(grepl("\\|", self$formula_terms))) {
        log4r::error(self$log, "The model contains at least one random effect. Are you sure you wanted linear models?")
        stop()
      }
    } else {
      log4r::error(self$log, "You entered an undefined method. Use `LMM` or `LM`!")
      stop()
    }
  }
  if (self$use_Rfast) {
    log4r::info(self$log, "Will use Rfast!")
  }
  #invisible(self)
}

#' Convert wide df to long df
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
wide_to_long <- function() {
  if (self$wide) {
    log4r::info(self$log, "Converting from wide to long!")
    self$df <- melt(self$df, id.vars = c(self$all_formula_terms[!self$all_formula_terms %in% c("variable", "value")]), value.name = self$valCol)
    log4r::info(self$log, paste0("Found ",length(unique(self$df$variable))," variables"))
    self$wide <- FALSE
    self$variablelist <- unique(self$df$variable)
    self$stratification_vector = self$df[, get(self$stratification_column)]
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
  self$df_raw <- NULL
  self$parts <- NULL
  self$validation_participants <- NULL
  self$stratification_vector <- NULL
  self$parts_with_variable <- NULL
  self$validation_self <- NULL
  self$regression_model <- NULL
  # self$regression_coefficients <- NULL
  self$effect_matrix <- NULL

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
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value[group == levels(group)[1]]), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1]]), by = variable]
      }
    }
  } else if (scale_function_string == "sdt1") {
    if (scale_function.center) {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[time == levels(time)[1]]), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := value / sd(value[time == levels(time)[1]]), by = variable]
      }
    }
  } else if (scale_function_string == "sdreft1") {
    if (scale_function.center) {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
      }
    } else {
      self$scale_function <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
      }
    }
  } else {
    log4r::error(self$log, "Unknown scaling method. Please use of one the following: `none`, `sdall`, `sdref`, `sdreft1`, `sdt1`")
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
      log4r::info(self$log, "Scaling data with custom function...")
    }
  } else if (self$scale_function == "none") {
    # The user do not want to scale
    if (!self$minimize_object) {
      log4r::warn(self$log, "Not scaling data...")
    }
    self$scale_function <- identity()
  } else if (is.character(self$scale_function)) {
    # Use a deafult scaling
    if (!self$minimize_object) {
      log4r::info(self$log, paste("Scaling data with",self$scale_function,"..."))
    }
    self$get_default_scaling_function()
  } else {
    log4r::error(self$log, "Unknown scaling function")
    stop()
  }
  #invisible(self)
}

#' Flip an ALASCA object
#'
#' Changes the sign of loadings and scores
#'
#' @param object An ALASCA object
#' @param component Components to be flipped, `NA` flips all (default)
#' @param effect Specify `time` or `group` to only flip subplot
#' @return An ALASCA object
#' @export
flip <- function(object, component = NA, effect = "both") {
  if (any(is.na(component))) {
    component <- unique(object$ALASCA$score$time$PC)
  }

  if (effect %in% c("both", "time")) {
    object$ALASCA$score$time[PC %in% component, score := -score]
    object$ALASCA$loading$time[PC %in% component, loading := -loading]
  }
  if (object$separate_time_and_group & effect %in% c("both", "group")) {
    object$ALASCA$score$group[PC %in% component, score := -score]
    object$ALASCA$loading$group[PC %in% component, loading := -loading]
  }
  if (object$validate) {
    if (effect %in% c("both", "time")) {
      object$ALASCA$score$time[PC %in% component, low := -low]
      object$ALASCA$score$time[PC %in% component, high := -high]
      object$ALASCA$loading$time[PC %in% component, low := -low]
      object$ALASCA$loading$time[PC %in% component, high := -high]
    }
    if (object$separate_time_and_group & effect %in% c("both", "group")) {
      object$ALASCA$score$group[PC %in% component, low := -low]
      object$ALASCA$score$group[PC %in% component, high := -high]
      object$ALASCA$loading$group[PC %in% component, low := -low]
      object$ALASCA$loading$group[PC %in% component, high := -high]
    }
  }

  if (object$validate && !object$save_to_disk) {
    for (i in seq_along(object$validation$temp_objects)) {
      object$validation$temp_objects[[i]] <- flip(object$validation$temp_objects[[i]], component = component, effect = effect)
    }
  }
  return(object)
}

#' Summary
#'
#' Gives some general information
#'
#' @param object
#' @param file
#' @param sessioninfo
#' @export
summary.ALASCA <- function(object, file = "", sessioninfo = FALSE) {
  cat("================ ALASCA ================\n", file = file, append = TRUE)
  cat("Model initialized ", as.character(object$init_time), " using ", object$method, " on ", length(unique(object$regression_coefficients$covar)), " variables. ", sep = "", file = file, append = TRUE)
  if (object$validate) {
    cat("The model been validated with ", object$validation_method, ".\n", file = file, append = TRUE)
  } else {
    cat("The model has *not* been validated yet.\n", file = file, append = TRUE)
  }
  cat("\nRegression model: ", as.character(object$rawFormula)[2], as.character(object$rawFormula)[1], as.character(object$rawFormula)[3], "\n", file = file, append = TRUE)
  cat("Separating time and group effects: ", object$separate_time_and_group, "\n", file = file, append = TRUE)
  cat("Force equal baseline: ", object$equal_baseline, "\n", file = file, append = TRUE)
  cat("Using Rfast: ", object$use_Rfast, "\n", file = file, append = TRUE)
  if (!object$use_Rfast) {
    cat("Adjustment of p-values: ", object$p_adjust_method, "\n", file = file, append = TRUE)
  }
  if (!is.null(object$reduce_dimensions.nComps)) {
    cat("Kept",object$reduce_dimensions.nComps,"components from initial PCA, explaining",100*cumsum(object$limm.explanatory_power)[object$reduce_dimensions.nComps],"% of variation\n", file = file, append = TRUE)
  }
  cat("\n\nPCs explaining at least 5% of variation:\n   Time: ",
    paste(get_relevant_pcs(object = object, effect = "time"), collapse = ", "), " (",
    paste(round(100 * object$pca$score$explained$time[get_relevant_pcs(object = object, effect = "time")], 2), collapse = "%, "), "%)",
    ifelse(object$separate_time_and_group, paste0(
      "\n   Group: ",
      paste(get_relevant_pcs(object = object, effect = "group"), collapse = ", "), " (",
      paste(round(100 * object$pca$score$explained$group[get_relevant_pcs(object = object, effect = "group")], 2), collapse = "%, "), "%)"
    ), "\n"),
    sep = "", file = file, append = TRUE
  )
  if (is.function(object$scale_function)) {
    cat("\nScaling function:\n", file = file, append = TRUE)
    cat(deparse(object$scale_function), file = file, append = TRUE)
  } else {
    cat("\nNo scaling performed.\n", file = file, append = TRUE)
  }
  if (sessioninfo & file != "") {
    write(capture.output(devtools::session_info()), file = file, append = TRUE)
  }
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
      log4r::error(self$log, "Unknown PCA function")
      stop()
    }
  } else {
    log4r::error(self$log, "Unknown PCA function")
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
  
  if (!self$minimize_object) {
    # This is not a validation run
    log4r::info(self$log, "Finished calculating regression coefficients!")
  }
  
  if (self$do_debug) currentTs <- Sys.time()
  self$remove_covars()
  if (self$do_debug) cat("* removeCovars:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$get_effect_matrix()
  if (self$do_debug) cat("* getEffectMatrix:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$do_pca()
  if (self$do_debug) cat("* doPCA:", Sys.time() - currentTs, "s\n")
  if (self$do_debug) currentTs <- Sys.time()
  self$clean_pca()
  if (self$do_debug) cat("* clean_pca:", Sys.time() - currentTs, "s\n")
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
  log4r::debug(self$log, "Find rows by variable")
  rows_by_variable <- lapply(self[["variablelist"]], function(x) self[["df"]][variable == x, which = TRUE] )
  names(rows_by_variable) <- self[["variablelist"]]
  
  if (self$use_Rfast && self$method %in% c("LMM")) {
    # start.time <- Sys.time()
    if (any(is.na(self[["df"]][, value]))) {
      log4r::error(self$log, "Rfast does NOT like NA's! Check your scaling function or value column.")
      stop()
    }
    log4r::debug(self$log, "Make model matrix")
    if (!self[["minimize_object"]]) {
      self[["modmat"]] <- model.matrix(self[["new_formula"]], data = self$df)
      if (self[["equal_baseline"]]) {
        self[["modmat"]] <- self[["modmat"]][, !grepl(paste0("time", self[["timelist"]][1]), colnames(self[["modmat"]]))]
      }
      self[["cnames_modmat"]] <- colnames(self[["modmat"]])
    }
    log4r::debug(self$log, "Finished model matrix")
    
    log4r::debug(self$log, "Starting regression")
    # https://stackoverflow.com/questions/61013078/fastest-way-to-convert-a-matrix-to-a-data-table
    self$regression_coefficients <- setDT(as.data.frame(
      vapply(self[["variablelist"]], function(x) {
        Rfast::rint.reg(
          y = self[["df"]][ rows_by_variable[[x]], value],
          x = self[["modmat"]][rows_by_variable[[x]], -1],
          id = as.numeric(factor(self[["df"]][ rows_by_variable[[x]], ID])),
          ranef = FALSE
        )$be
      }, FUN.VALUE = numeric(ncol(self[["modmat"]])))
    ))[]
    log4r::debug(self$log, "Finished regression")
    colnames(self$regression_coefficients) <- as.character(self[["variablelist"]])
    self$regression_coefficients[, variable := self[["cnames_modmat"]]]
    self$regression_coefficients <- melt(self$regression_coefficients, id.vars = "variable", variable.name = "covar", variable.factor = FALSE, value.name = "estimate")
    self$regression_coefficients[, pvalue := NA]
    
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
    names(self$regression_model) <- self$variablelist
  } else if (self$use_Rfast & self$method %in% c("LM")) {
    # start.time <- Sys.time()
    if (any(is.na(self$df[, value]))) {
      log4r::error(self$log, "Rfast does NOT like NA's! Check your scaling function or value column.")
      stop()
    }
    if (!self$minimize_object) {
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
    names(self$regression_model) <- self$variablelist
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
  if (!self$minimize_object) {
    log4r::info(self$log, "Calculating effect matrix")
  }
  
  reg_coefs <- dcast(self$regression_coefficients, variable~covar, value.var = "estimate")
  rownames(reg_coefs) <- reg_coefs$variable
  self$effect_list$effect_matrix <- lapply(self$set_design_matrices(), function(mm) {
      as.matrix(mm) %*% as.matrix(reg_coefs[colnames(mm), -1])
    }
  )
  
  if (!self$minimize_object) {
    log4r::info(self$log, "Finished calculating effect matrix!")
  }
  #invisible(self)
}

#' Save ALASCA object
#'
#' @inheritParams saveALASCAModel
#' @inheritParams save_to_csv
#' @return An ALASCA object
#' @export
saveALASCA <- function(object, filename = NA, filepath = NA, saveCSV = TRUE, saveScores = TRUE, saveLoadings = TRUE, saveCovars = TRUE, csv = "csv", ...) {
  saveALASCAModel(object = object, filename = filename, filepath = NA)
  save_to_csv(object = object, filename = filename, filepath = filepath, saveCSV = saveCSV, saveScores = saveScores, saveLoadings = saveLoadings, saveCovars = saveCovars, csv = "csv", ...)
  summary.ALASCA(object = object, file = get_filename(object = object, filetype = "txt"), sessioninfo = TRUE)
}

#' Save ALASCA object
#'
#' @inheritParams saveALASCAModel
#' @inheritParams save_to_csv
#' @return An ALASCA object
#' @export
save_bootstrap_ids <- function(object) {
  if (!(object$validate && object$validation_method == "bootstrap")) {
    stop("Please validate with bootstrapping")
  }
  for (i in seq_along(object$validation$temp_object)) {
    write(paste0(object$validation$temp_object[[i]]$validation$original_ids, collapse = ";"),
          file = get_filename(object = object, prefix = "bootstrapID_", filetype = ".csv", overwrite = TRUE), append = TRUE
    )
  }
}

#' Save ALASCA object
#'
#' @inheritParams saveALASCAModel
#' @inheritParams save_to_csv
#' @return An ALASCA object
#' @export
save_jackknife_ids <- function(object) {
  if (!(object$validate && object$validation_method %in% c("loo", "jack-knife", "jackknife"))) stop("Please validate with jack-knife")
  for (i in seq_along(object$validation$temp_object)) {
    write(paste0(unique(object$validation$temp_object[[i]]$partID), collapse = ";"),
          file = get_filename(object = object, prefix = "jackknifeID_", filetype = ".csv", overwrite = TRUE), append = TRUE
    )
  }
}

#' Save ALASCA object to csv
#'
#' @param object An ALASCA object
#' @param filepath
#' @param filename
#' @export
save_to_csv <- function(object, filename = NA, filepath = NA, saveCSV = TRUE, saveScores = TRUE, saveLoadings = TRUE, saveCovars = TRUE, csv = "csv", ...) {
  if (saveCSV) {
    if (csv == "csv") {
      if (object$separate_time_and_group) {
        if (saveScores) {
          write.csv(get_scores(object)$time,
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_time_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv(get_loadings(object)$time,
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_time_loadings", filetype = ".csv"), ...
          )
        }
        if (saveScores) {
          write.csv(get_scores(object)$group,
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_group_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv(get_loadings(object)$group,
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_group_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv(get_covars(object),
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      } else {
        if (saveScores) {
          write.csv(get_scores(object)$time,
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv(get_loadings(object)$time,
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv(get_covars(object),
                    file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      }
    } else if (csv == "csv2") {
      if (object$separate_time_and_group) {
        if (saveScores) {
          write.csv2(get_scores(object)$time,
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_time_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv2(get_loadings(object)$time,
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_time_loadings", filetype = ".csv"), ...
          )
        }
        if (saveScores) {
          write.csv2(get_scores(object)$group,
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_scores_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv2(get_loadings(object)$group,
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_scores_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv2(get_covars(object),
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      } else {
        if (saveScores) {
          write.csv2(get_scores(object)$time,
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv2(get_loadings(object)$time,
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv2(get_covars(object),
                     file = get_filename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      }
    } else {
      stop("Unkown input to `csv`")
    }
    cat(paste0("- Saved csv's to ", get_filename(object = object, filename = filename, filepath = filepath), "\n"))
  }
}

#' Save ALASCA object
#'
#' @param object An ALASCA object
#' @param filepath
#' @param filename
#' @return Full file name of the saved object (String)
#' @export
saveALASCAModel <- function(object, filename = NA, filepath = NA) {
  fname <- get_filename(object = object, filename = filename, filepath = filepath, filetype = ".rds")
  saveRDS(object, file = fname)
  cat(paste0("- Saved model to ", fname, "\n"))
  return(fname)
}

#' Save figure
#'
#' @param g ggplot-object
#' @return A ggplot2 objects.
#'
#' @export
saveALASCAPlot <- function(object, g, filetype = NA, figsize = NA, prefix = "plot/", suffix = "", figunit = NA) {
  if (any(is.na(filetype))) {
    filetype <- object$plot.filetype
  }
  if (any(is.na(figsize))) {
    figsize <- object$plot.figsize
  }
  if (is.na(figunit)) {
    figunit <- object$plot.figunit
  }
  
  for (i in filetype) {
    fname <- get_filename(object = object, prefix = prefix, suffix = suffix, filetype = i)
    ggplot2::ggsave(plot = g, filename = fname, width = figsize[1], height = figsize[2], dpi = figsize[3], unit = figunit)
    cat(paste0("- Saved ", fname, "\n"))
  }
}

#' Get file path
#'
#' @param object An ALASCA object
#' @param filename File name
#' @param filepath File path
#' @param prefix Prefix
#' @param suffix Suffix
#'
#' @return A ggplot2 objects.
#'
#' @export
get_filename <- function(object, filename = NA, filepath = NA, prefix = "", suffix = "", filetype = "", overwrite = FALSE) {
  # Use arguments if defined
  if (any(!is.na(filename))) {
    object$filename <- filename
  }
  if (any(!is.na(filepath))) {
    object$filepath <- filepath
  }
  
  # Set defaults if undefined
  if (any(is.na(object$filepath))) {
    object$filepath <- paste0("ALASCA/", strftime(object$initTime, format = "%Y%m%d_%H%M%S"), "/")
  }
  if (any(is.na(object$filename))) {
    object$filename <- "ALASCA"
  }
  
  # Create folder
  if (!dir.exists(object$filepath)) {
    dir.create(object$filepath, recursive = TRUE)
  }
  
  # Check if prefix implies subfolder
  if (grepl("/", prefix)) {
    if (!dir.exists(paste0(object$filepath, prefix))) {
      dir.create(paste0(object$filepath, prefix), recursive = TRUE)
    }
  }
  
  fname <- paste0(object$filepath, prefix, object$filename)
  
  # Check if file already exists
  if (!overwrite) {
    cnt <- 1
    while (file.exists(paste0(fname, suffix, ifelse(substr(filetype, 1, 1) == ".", filetype, paste0(".", filetype))))) {
      fname <- paste0(object$filepath, prefix, object$filename, "_", cnt)
      cnt <- cnt + 1
    }
  }
  fname <- paste0(fname, suffix, ifelse(substr(filetype, 1, 1) == "." | filetype == "", filetype, paste0(".", filetype)))
  return(fname)
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
    function(x) prcomp(x, scale = FALSE, center = !self$scale_function.center)
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
  
  log4r::debug(self$log, "Starting reduce_dimensions")
  
  wide_data <- dcast(data = self$df, ...~variable
  )
  if (self$do_debug) currentTs <- Sys.time()
  temp_pca_values <- self$function.pca(
    wide_data[, .SD, .SDcols = -self$all_formula_terms],
    center = !self$scale_function.center
  )
  
  self$reduce_dimensions.explanatory_power <- temp_pca_values$sdev^2 / sum(temp_pca_values$sdev^2)
  
  # Remove surplus columns
  if (is.null(self$reduce_dimensions.nComps)) {
    self$reduce_dimensions.nComps <- which(cumsum(self$reduce_dimensions.explanatory_power) >= self$reduce_dimensions.limit)[1]
    log4r::info(self$log, paste("Keeping",
                                self$reduce_dimensions.nComps,
                                "components from initial PCA, explaining",
                                100*cumsum(self$reduce_dimensions.explanatory_power)[self$reduce_dimensions.nComps],
                                "% of variation"))
  }
  if(ncol(temp_pca_values$rotation) > self$reduce_dimensions.nComps){
    temp_pca_values$rotation <- temp_pca_values$rotation[,-c((self$reduce_dimensions.nComps+1):ncol(temp_pca_values$rotation))]
    temp_pca_values$x <- temp_pca_values$x[, -c((self$reduce_dimensions.nComps+1):ncol(temp_pca_values$x))]
  }
  
  # Check if the pca model needs reflection to better fit the main model
  if (self$do_debug) currentTs <- Sys.time()
  for (i in seq_len(ncol(temp_pca_values$rotation))) {
    V1 <- sum((temp_pca_values$rotation[,i] - self$Limm$main$pca$rotation[,i])^2)
    V2 <- sum((-temp_pca_values$rotation[,i] - self$Limm$main$pca$rotation[,i])^2)
    if(V2 < V1){
      temp_pca_values$rotation[,i] = -temp_pca_values$rotation[,i]
      temp_pca_values$x[,i] = -temp_pca_values$x[,i]
    }
  }
  log4r::info(self$log, "Removing surplus PCs")
  
  self$Limm$loadings <- temp_pca_values$rotation
  self$Limm$pca <- temp_pca_values
  self$Limm$df <- self$df
  self$df <- melt(data = cbind(wide_data[, .SD, .SDcols = self$all_formula_terms], temp_pca_values$x), id.vars = self$all_formula_terms, variable.factor = FALSE)
  self$variablelist <- unique(self$df$variable)
  self$stratification_vector <- self$df[, get(self$stratification_column)]
  log4r::info(self$log, "Completed reduce_dimensions")
  #invisible(self)
}

#' Clean the PCA data
#'
#' This function makes the pca output more useful
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
clean_pca <- function() {
  log4r::debug(self$log, "Starting clean_pca")
  
  self$ALASCA$score <- lapply(seq_along(self$effect_list$pca), function(i){
    unique(
      cbind(
        self$df[variable == self$variablelist[[1]], .SD, .SDcols = self$effect_list$terms[[i]]],
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
    
  if(self$reduce_dimensions){
    # Loadings must be back-transformed
    self$pca$loading$time <- setDT(
      as.data.frame(
        as.matrix(self$Limm$loadings) %*% as.matrix(self$pca$loading$time[order(as.numeric(substr(covars, 3, nchar(covars)))), !"covars"])
      ),
      keep.rownames="covars")
    self$variablelist <- unique(self$pca$loading$time$covars)
    self$regression_coefficients <- rbindlist(
      lapply(unique(self$regression_coefficients$variable), function(x){
        setDT(
          data.frame(
            variable = x,
            pvalue = NA,
            estimate = as.matrix(self$Limm$loadings) %*% as.matrix(self$regression_coefficients[variable == x & order(as.numeric(substr(covar, 3, nchar(covar)))), "estimate"])
          ),
          keep.rownames="covar")
      })
    )
  }
  
  # Ensure that the highest loading has positive sign
  for (i in seq_along(self$ALASCA$loading)) {
    setkeyv(self$ALASCA$loading[[i]], cols = "covars")
    setkeyv(self$ALASCA$score[[i]], cols = self$effect_list$terms[[i]])
    for(k in colnames(self$ALASCA$loading[[i]][, -1])) {
      nVar <- self$ALASCA$loading[[i]][, .I[which.max(abs(get(k)))]]
      sVar <- self$ALASCA$loading[[i]][nVar, sign(get(k))]
      
      set(self$ALASCA$loading[[i]], j = (k), value = self$ALASCA$loading[[i]][, get(k)] * sVar)
      set(self$ALASCA$score[[i]], j = (k), value = self$ALASCA$score[[i]][, get(k)] * sVar)
    }
  }
  
  log4r::debug(self$log, "Completed clean_pca")
  
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
do_validate <- function(participant_column = FALSE, validate_regression = FALSE) {
  if (self$validate) {
    # stop("The object has already been validated")
  }
  self$validate <- TRUE
  if (validate_regression) {
    self$validate_regression <- TRUE
  }
  
  if (self$method == "LMM") {
    if (participant_column == FALSE) {
      if (self$participant_column == FALSE) {
        log4r::error(self$log, "You need to specify the column containing participant id in `participant_column`")
        stop()
      }
    } else {
      self$participant_column <- participant_column
    }
  }
  
  log4r::info(self$log, "Starting validation")
  
  start.time.all <- Sys.time()
  
  if (any(is.na(self$validation_ids))) {
    # Generate random samples
    if (self$save_to_disk) {
      limPC_time <- self$get_relevant_pcs(effect = "time")
      if (object$separate_time_and_group) {
        limPC_group <- self$get_relevant_pcs(effect = "group")
      }
      temp_object <- lapply(seq_len(self$n_validation_runs), FUN = function(ii) {
        log4r::info(self$log, paste0("- Run ", ii, " of ", self$n_validation_runs))
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- self$prepare_validation_run()
        
        # Rotate new loadings/scores to the original model
        if (self$optimize_score) {
          temp_object$rotate_matrix_optimize_score(target = self)
        } else {
          temp_object$rotate_matrix(target = self)
        }
        temp_object$clean_alasca()
        
        # Save to disk
        fname <- paste0("val_", ii)
        loading <- data.frame(temp_object$ALASCA$loading$time[temp_object$ALASCA$loading$time$PC %in% limPC_time, ], model = ii)
        # loading$covars <- object$numvariablelist[loading$covars]
        DBI::dbWriteTable(object$db.con, "time.loading", loading, append = T, overwrite = F)
        scores <- data.frame(temp_object$ALASCA$score$time[temp_object$ALASCA$score$time$PC %in% limPC_time, ], model = ii)
        # scores$time <- object$numtimelist[scores$time]
        # if(!object$separate_time_and_group){
        #  scores$group <- object$numgrouplist[scores$group]
        # }
        DBI::dbWriteTable(object$db.con, "time.score", scores, append = T, overwrite = F)
        if (nrow(get_covars(object)) > 0) DBI::dbWriteTable(object$db.con, "covars", data.frame(get_covars(temp_object), model = ii), append = T, overwrite = F)
        
        if (object$separate_time_and_group) {
          loading <- data.frame(temp_object$ALASCA$loading$group[temp_object$ALASCA$loading$group$PC %in% limPC_group, ], model = ii)
          # loading$covars <- object$numvariablelist[loading$covars]
          DBI::dbWriteTable(object$db.con, "group.loading", loading, append = T, overwrite = F)
          
          scores <- data.frame(temp_object$ALASCA$score$group[temp_object$ALASCA$score$group$PC %in% limPC_group, ], model = ii)
          # scores$time <- object$numtimelist[scores$time]
          # scores$group <- object$numgrouplist[scores$group]
          DBI::dbWriteTable(object$db.con, "group.score", scores, append = T, overwrite = F)
        }
        if (object$validate_regression) {
          DBI::dbWriteTable(object$db.con, "model_prediction", data.frame(temp_object$model_prediction, model = ii), append = T, overwrite = F)
        }
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(object$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((object$n_validation_runs - ii) * time_all, 2), " seconds"))
        fname
      })
    } else {
      temp_object <- lapply(seq_len(self$n_validation_runs), FUN = function(ii) {
        
        log4r::info(self$log, paste0("- Run ", ii, " of ", self$n_validation_runs))
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- self$prepare_validation_run()
        if (nrow(self$df) > 100000) gc()
        
        # Rotate new loadings/scores to the original model
        if (self$optimize_score) {
          temp_object$rotate_matrix_optimize_score(target = self)
        } else {
          temp_object$rotate_matrix(target = self)
        }
        temp_object$clean_alasca()
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(self$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((self$n_validation_runs - ii) * time_all, 2), " seconds"))
        temp_object
      })
    }
  } else {
    # Reuse previous samples
    log4r::info(self$log, "Using predefined samples")
    
    if (self$save_to_disk) {
      limPC_time <- get_relevant_pcs(object = object, effect = "time")
      if (object$separate_time_and_group) {
        limPC_group <- get_relevant_pcs(object = object, effect = "group")
      }
      temp_object <- lapply(seq_len(object$n_validation_runs), FUN = function(ii) {
        log4r::info(object$log, paste0("- Run ", ii, " of ", object$n_validation_runs))
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- prepare_validation_run(object, runN = ii)
        gc()
        
        # Rotate new loadings/scores to the original model
        if (object$optimize_score) {
          temp_object <- rotate_matrix_optimize_score(object = temp_object, target = object)
        } else {
          temp_object <- rotate_matrix(object = temp_object, target = object)
        }
        
        # Save to disk
        fname <- paste0("val_", ii)
        loading <- data.frame(temp_object$ALASCA$loading$time[temp_object$ALASCA$loading$time$PC %in% limPC_time, ], model = ii)
        # loading$covars <- object$numvariablelist[loading$covars]
        DBI::dbWriteTable(object$db.con, "time.loading", loading, append = T, overwrite = F)
        scores <- data.frame(temp_object$ALASCA$score$time[temp_object$ALASCA$score$time$PC %in% limPC_time, ], model = ii)
        # scores$time <- object$numtimelist[scores$time]
        # if(!object$separate_time_and_group){
        #  scores$group <- object$numgrouplist[scores$group]
        # }
        DBI::dbWriteTable(object$db.con, "time.score", scores, append = T, overwrite = F)
        if (nrow(get_covars(object)) > 0) DBI::dbWriteTable(object$db.con, "covars", data.frame(get_covars(temp_object), model = ii), append = T, overwrite = F)
        
        if (object$separate_time_and_group) {
          loading <- data.frame(temp_object$ALASCA$loading$group[temp_object$ALASCA$loading$group$PC %in% limPC_group, ], model = ii)
          # loading$covars <- object$numvariablelist[loading$covars]
          DBI::dbWriteTable(object$db.con, "group.loading", loading, append = T, overwrite = F)
          
          scores <- data.frame(temp_object$ALASCA$score$group[temp_object$ALASCA$score$group$PC %in% limPC_group, ], model = ii)
          # scores$time <- object$numtimelist[scores$time]
          # scores$group <- object$numgrouplist[scores$group]
          DBI::dbWriteTable(object$db.con, "group.score", scores, append = T, overwrite = F)
        }
        if (object$validate_regression) {
          DBI::dbWriteTable(object$db.con, "model_prediction", data.frame(temp_object$model_prediction, model = ii), append = T, overwrite = F)
        }
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(object$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((object$n_validation_runs - ii) * time_all, 2), " seconds"))
        fname
      })
    } else {
      temp_object <- lapply(seq_len(self$n_validation_runs), FUN = function(ii) {
        log4r::info(self$log, paste0("- Run ", ii, " of ", self$n_validation_runs))
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- self$prepare_validation_run(runN = ii)
        gc()
        
        # Rotate new loadings/scores to the original model
        if (self$optimize_score) {
          temp_object$rotate_matrix_optimize_score(target = self)
        } else {
          temp_object$rotate_matrix(target = self)
        }
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(self$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((self$n_validation_runs - ii) * time_all, 2), " seconds"))
        temp_object
      })
    }
  }
  
  self$clean_alasca()
  
  log4r::info(self$log, "Calculating percentiles for score and loading")
  self$get_validation_percentiles(objectlist = temp_object)
  
  if (self$keep_validation_objects) {
    self$validation$temp_objects <- temp_object
  }
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

#' Rotate PCA
#'
#' This function rotates loadings and scores during validation
#'
#' Optimizes the rotation for lowest possible difference in score
#'
#' @param object ALASCA object to be rotated (and returned)
#' @param target ALASCA object acting as target
#' @return An ALASCA object
rotate_matrix_optimize_score <- function(target) {
  log4r::debug(target$log, "Starting rotation")
  
  # PCA can give loadings with either sign, so we have to check whether switching signs improves the rotation
  for (effect_i in seq_along(self$ALASCA$loading)) {
    # Number of components to look at
    PCs_to_look_at <- 1:2
    cols_to_look_at <- paste0("PC", PCs_to_look_at)
    N <- length(PCs_to_look_at)
    
    # Create matrix with all possible combinations of signs
    vec <- c(-1, 1)
    lst <- lapply(numeric(N), function(x) vec)
    signMatrix <- as.matrix(expand.grid(lst))
    
    # Test all combinations and calculate residuals
    signVar <- vapply(seq_len(nrow(signMatrix) / 2), function(i) {
      c <- .procrustes(
        loadings = as.matrix(
          t(t(self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], ..cols_to_look_at]) * signMatrix[i, ])
          ),
        target = as.matrix(target$ALASCA$loading[[effect_i]][, ..cols_to_look_at])
      )
      sum((target$ALASCA$score[[effect_i]][, ..cols_to_look_at] - 
             as.matrix(t(t(self$ALASCA$score[[effect_i]][target$ALASCA$score[[effect_i]], ..cols_to_look_at]) * signMatrix[i, ])) %*% solve(c$t1))^2)
    }, FUN.VALUE = numeric(1))
    
    # Find the combination that minimizes the sum of squares
    minSignVar <- which(signVar == min(signVar))[1]
    
    # Switch signs
    for (i in PCs_to_look_at){
      set(self$ALASCA$loading[[effect_i]], j = cols_to_look_at[i], value = self$ALASCA$loading[[effect_i]][, get(cols_to_look_at[i])] * signMatrix[minSignVar, i])
      set(self$ALASCA$score[[effect_i]],   j = cols_to_look_at[i], value = self$ALASCA$score[[effect_i]][, get(cols_to_look_at[i])] * signMatrix[minSignVar, i])
    }
    
    # Final rotation
    c <- .procrustes(
      loadings = as.matrix(self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], ..cols_to_look_at]),
      target = as.matrix(target$ALASCA$loading[[effect_i]][, ..cols_to_look_at])
    )
    
    self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], (cols_to_look_at) := as.data.table(c$procrust)]
    self$ALASCA$score[[effect_i]][target$ALASCA$score[[effect_i]], (cols_to_look_at) := as.data.table(as.matrix(.SD) %*% solve(c$t1)), .SDcols = cols_to_look_at]
    
  }
  
  log4r::debug(target$log, "Completed rotation")

  #invisible(self)
}

#' Rotate PCA
#'
#' This function rotates loadings and scores during validation
#'
#' @param object ALASCA object to be rotated (and returned)
#' @param target ALASCA object acting as target
#' @return An ALASCA object
rotate_matrix <- function(target) {
  PCloading <- target$get_relevant_pcs(effect = "time")
  PCloading_t <- paste0("PC", PCloading)
  c <- .procrustes(
    loadings = as.matrix(self$pca$loading$time[target$pca$loading$time, ..PCloading_t]),
    target = as.matrix(target$pca$loading$time[, ..PCloading_t])
  )
  
  self$pca$loading$time[target$pca$loading$time, (PCloading_t) := as.data.frame(c$procrust)]
  self$pca$score$time[target$pca$score$time, (PCloading_t) := as.data.frame(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]
  
  if (self$separate_time_and_group) {
    PCloading <- taret$get_relevant_pcs(effect = "group")
    PCloading_t <- paste0("PC", PCloading)
    c <- .procrustes(
      loadings = as.matrix(self$pca$loading$group[target$pca$loading$group, ..PCloading_t]),
      target = as.matrix(target$pca$loading$group[, ..PCloading_t])
    )
    
    self$pca$loading$group[target$pca$loading$group, (PCloading_t) := as.data.frame(c$procrust)]
    self$pca$score$group[target$pca$score$group, (PCloading_t) := as.data.frame(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]
  }
  
  #invisible(self)
}

#' Extract percentiles
#'
#' This function extract percentiles during validation
#'
#' @param object ALASCA object
#' @param objectlist List of ALASCA objects
#' @return An ALASCA object
get_validation_percentiles <- function(objectlist) {
  
  log4r::debug(self$log, message = "Starting get_validation_percentiles_loading")
  self$get_validation_percentiles_loading(objectlist)
  log4r::debug(self$log, message = "Completed get_validation_percentiles_loading")
  log4r::debug(self$log, message = "Starting get_validation_percentiles_score")
  self$get_validation_percentiles_score(objectlist)
  log4r::debug(self$log, message = "Completed get_validation_percentiles_score")
  if (self$validate_regression) {
    log4r::debug(self$log, message = "Starting get_validation_percentiles_regression")
    self$get_validation_percentiles_regression(objectlist)
    log4r::debug(self$log, message = "Completed get_validation_percentiles_regression")
  }
  if (nrow(self$get_covars()) > 0) {
    log4r::debug(self$log, message = "Starting get_validation_percentiles_covars")
    self$get_validation_percentiles_covars(objectlist)
    log4r::debug(self$log, message = "Completed get_validation_percentiles_covars")
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
    res <- DBI::dbSendQuery(object$db.con, "SELECT * FROM 'model_prediction'")
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df <- rbindlist(lapply(objectlist, function(x) x$model_prediction))
  }
  df <- df[, as.list(quantile(pred, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(group, time, variable)]
  colnames(df) <- c("group", "time", "variable", "low", "high")
  
  self$model_prediction <- merge(self$model_prediction, df)
  #invisible(self)
}

#' Extract percentiles for covariates
#'
#' This function extract percentiles for validation of covariates
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_covars <- function(objectlist) {
  if (self$save_to_disk) {
    res <- DBI::dbSendQuery(self$db.con, paste0("SELECT * FROM 'covars'"))
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df <- rbindlist(lapply(objectlist, function(x) get_covars(x)))
  }
  df <- df[, as.list(quantile(estimate, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(covar, variable)]
  colnames(df) <- c("covar", "variable", "low", "high")
  
  self$covar_coefficients <- merge(self$covar_coefficients, df)
  #invisible(self)
}

#' Extract percentiles for loading
#'
#' This function extract percentiles during validation of loadings
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_loading <- function(objectlist) {
  
  PC_time <- 1:2 #self$get_relevant_pcs(effect = "time")
  
  for (effect_i in seq_along(self$ALASCA$loading)){
    if (self$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.loading' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_time <- rbindlist(
        lapply(objectlist, function(x) x$ALASCA$loading[[effect_i]][x$ALASCA$loading[[effect_i]]$PC %in% PC_time, ]),
        fill = TRUE
        )
    }
    
    tmp <- df_time[, as.list(
      quantile(loading, probs = self$limitsCI, type = self$validation_quantile_method)
      ), by = .(PC, covars)]
    colnames(tmp) <- c("PC", "covars", "low", "high")
    self$ALASCA$loading[[effect_i]] <- merge(tmp, self$ALASCA$loading[[effect_i]], by = c("PC", "covars"))
  }
  #invisible(self)
}

#' Extract percentiles for score
#'
#' This function extract percentiles during validation of scores
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_score <- function(objectlist) {
  
  PC_time <- 1:2 #self$get_relevant_pcs(effect = "time")
  
  for (effect_i in seq_along(self$ALASCA$score)) {
    
    if (self$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score[[effect_i]][x$ALASCA$score[[effect_i]]$PC %in% PC_time, ]), fill = TRUE)
    }
    
    tmp <- df_time[, as.list(
      quantile(score, probs = self$limitsCI, type = self$validation_quantile_method)
      ), by = c("PC", self$effect_list$terms[[effect_i]])]
    colnames(tmp) <- c("PC", self$effect_list$terms[[effect_i]], "low", "high")
    self$ALASCA$score[[effect_i]] <- merge(tmp, self$ALASCA$score[[effect_i]], by = c("PC",  self$effect_list$terms[[effect_i]]))
  }

  #invisible(self)
}

#' Get relevant PCs
#'
#' This function extract percentiles during validation
#'
#' @param x Explanatory power of PC
#' @return A vector with relevant PCs
get_relevant_pcs <- function(effect = "time") {
  if (effect == "time") {
    PC <- self$ALASCA$loading$explained$time >= self$explanatorylimit
  } else {
    PC <- self$ALASCA$loading$explained$group >= self$explanatorylimit
  }
  PC[1:2] <- TRUE
  return(which(PC))
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
    log4r::info(self$log, message = "Calculating predictions from regression models")
  }
  regCoeffAll <- dcast(data = self[["regression_coefficients"]], variable~covar, value.var = "estimate")
  rownames(regCoeffAll) <- regCoeffAll$variable
  regModel <- unique(model.matrix(self$formula, data = self$df))
  if (self$equal_baseline) {
    regModel <- regModel[, !grepl(paste0("time", levels(self$df$time)[1]), colnames(regModel), fixed = TRUE), drop = FALSE]
  }
  regModel <- regModel[, !grepl("|", colnames(regModel), fixed = TRUE)]
  if (self$keep_terms != "") {
    regModel <- regModel[, grepl(paste0(c("time", "group", self$keep_terms), collapse = "|"), colnames(regModel)), drop = FALSE]
  } else {
    regModel <- regModel[, grepl(paste0(c("time", "group"), collapse = "|"), colnames(regModel)), drop = FALSE]
  }
  regModel <- unique(regModel)
  
  if (self$keep_terms != "") {
    self$model_prediction <- melt(
      cbind(
        as.matrix(regModel) %*% as.matrix(regCoeffAll[colnames(regModel), -1]),
        self$df[as.numeric(rownames(regModel)), .SD, .SDcols = c("time", "group", self$keep_terms)]
      ),
      id.vars = c("time", "group", self$keep_terms), value.name = "pred"
    )
  } else {
    self$model_prediction <- melt(
      cbind(
        as.matrix(regModel) %*% as.matrix(regCoeffAll[colnames(regModel), -1]),
        self$df[as.numeric(rownames(regModel)), .SD, .SDcols = c("time", "group")]
      ),
      id.vars = c("time", "group"), value.name = "pred"
    )
  }
  
  if (self$keep_terms != "") {
    self$model_prediction[, group := apply(.SD, 1, paste, collapse = " - "), .SDcols = c("group", self$keep_terms)]
    self$model_prediction[, group := factor(group, levels = self$grouplist)]
  } 
  
  if (!self$minimize_object) {
    # This is not a validation run
    log4r::info(self$log, message = "Finished calculating predictions from regression models!")
  }
  #invisible(self)
}

#' Make a single validation run
#'
#' This function ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
prepare_validation_run <- function(runN = NA) {
  if (self$validation_method %in% c("loo", "jack-knife", "jackknife")) {
    # Use leave-one-out validation
    selectedParts <- data.frame()
    
    if (self$method %in% c("LMM")) {
      if (any(is.na(self$validation_ids))) {
        # For each group, divide the participants into n_validation_folds groups, and select n_validation_folds-1 of them
        selectedParts <- lapply(unique(self$stratification_vector), function(gr) {
          selectedParts_temp_all <- unique(self$df[self$stratification_vector == gr, ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% self$n_validation_folds
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
        bootobject <- self$clone()
        bootobject$my_df_rows <- unlist(lapply(unlist(selectedParts), function(x) self$df_raw$rows_by_ID[[as.character(x)]]))
        bootobject[["df"]] <- NULL
        bootobject[["df"]] <- self$scale_function(self$df_raw$df[bootobject$my_df_rows])
        bootobject[["modmat"]] <- self$modmat[bootobject$my_df_rows,]
        bootobject$update()
        return(bootobject)
      } else {
        bootobject <- self$clone()
        bootobject$my_df_rows <- unlist(lapply(object$validation_ids[runN, ], function(x) self$df_raw$rows_by_ID[[as.character(x)]]))
        bootobject[["df"]] <- NULL
        bootobject[["df"]] <- self$scale_function(self$df_raw$df[bootobject$my_df_rows])
        bootobject[["modmat"]] <- self$modmat[bootobject$my_df_rows,]
        bootobject$update()
        return(bootobject)
      }
    } else if (self$method %in% c("LM")) {
      self$df$ID <- c(seq_len(nrow(self$df)))
      if (any(is.na(self$validation_ids))) {
        # For each group, divide the participants into n_validation_folds groups, and select n_validation_folds-1 of them
        selectedParts <- lapply(unique(self$stratification_vector), function(gr) {
          selectedParts_temp_all <- unique(self$df[object$stratification_vector == gr, ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% self$n_validation_folds
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
        
        bootobject <- self$clone()
        bootobject$my_df_rows <- unlist(lapply(unlist(selectedParts), function(x) self$df_raw$rows_by_ID[[as.character(x)]]))
        bootobject[["df"]] <- NULL
        bootobject[["df"]] <- self$scale_function(self$df_raw$df[bootobject$my_df_rows])
        bootobject[["modmat"]] <- self$modmat[bootobject$my_df_rows,]
        bootobject$update()
        return(bootobject)
      } else {
        bootobject <- self$clone()
        bootobject$my_df_rows <- unlist(lapply(object$validation_ids[runN, ], function(x) self$df_raw$rows_by_ID[[as.character(x)]]))
        bootobject[["df"]] <- NULL
        bootobject[["df"]] <- self$scale_function(self$df_raw$df[bootobject$my_df_rows])
        bootobject[["modmat"]] <- self$modmat[bootobject$my_df_rows,]
        bootobject$update()
        return(bootobject)
      }
    }
  } else if (self$validation_method == "bootstrap") {
    # Use bootstrap validation
    # When using bootstrapping, we resample participants with replacement
    
    log4r::debug(self$log, message = "Preparing bootstrap sample")
    participants_in_bootstrap <- self$get_bootstrap_ids(runN = runN)
    
    # Create bootstrap object without data
    bootobject <- self$clone()
    bootobject[["df"]] <- NULL
    bootobject$get_bootstrap_data(df_raw = self$df_raw, participants_in_bootstrap)
    
    if (self$validation_assign_new_ids) {
      bootobject$df[, ID := uniqueIDforBootstrap] # Replace ID
    }
    
    if (self$save_validation_ids) {
      write(paste0(participants_in_bootstrap$old_id, collapse = ";"),
            file = get_filename(object = object, prefix = "bootstrapID_", filetype = ".csv", overwrite = TRUE), append = TRUE
      )
    }
    
    log4r::debug(self$log, message = "Completed preparation of bootstrap sample")
    
    bootobject$update()
    bootobject$validation$original_ids <- participants_in_bootstrap$old_id
    return(bootobject)
  }
}

#' Make bootstrap sample
#'
#' Get data frame with new and old IDs
#'
#' @param object An ALASCA object
#' @return A data frame
get_bootstrap_ids <- function(runN = NA) {
  if (is.na(runN)) {
    participants_in_bootstrap <- data.frame(
      new_id = seq(self$df[, uniqueN(ID)]),
      old_id = rbindlist(
        lapply(unique(self$stratification_vector), function(stratification_group){
          list(
            old_id = sample(
              unique(self$df[self$stratification_vector == stratification_group, ID]),
              replace = TRUE)
          )
        })
      )
    )
  } else {
    participants_in_bootstrap <- data.frame(
      new_id = seq(self$validation_ids[runN, ]),
      old_id = matrix(self$validation_ids[runN, ])
    )
  }
  return(participants_in_bootstrap)
}

#' Make bootstrap data set
#'
#' Get data frame
#'
#' @param object An ALASCA object
#' @return A data frame
get_bootstrap_data <- function(df_raw, participants_in_bootstrap) {
  selected_rows <- rbindlist(
    lapply(participants_in_bootstrap$new_id, function(participant){
      list(
        new_id = participant,
        row_nr = df_raw$rows_by_ID[[participant]]
      )
    })
  )
  self[["my_df_rows"]] <- selected_rows$row_nr
  self$df <- self$scale_function(df_raw$df[self$my_df_rows])
  self$df[, uniqueIDforBootstrap := selected_rows$new_id]
  self$df[, originalIDbeforeBootstrap := ID]
  self$modmat <- self$modmat[selected_rows$row_nr,]
  #invisible(self)
}

#' Get covariables
#'
#' This function returns the other covariables in an ALASCA model
#'
#' @inheritParams get_loadings
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
get_covars <- function(n_limit = 0) {
  if (n_limit > 0) {
    return(
      rbind(
        self$covar_coefficients[order(estimate, decreasing = TRUE), head(.SD, n_limit), by = variable],
        self$covar_coefficients[order(estimate, decreasing = FALSE), head(.SD, n_limit), by = variable]
      )
    )
  } else {
    return(self$covar_coefficients)
  }
}

#' Get linetypes
#'
#' This function returns a list with linetypes for plotting
#'
#' @param object An ALASCA object
#' @return A list with linetypes
#'
#' @export
get_plot_linetypes <- function() {
  if (is.null(self$plot.linetypes)) {
    self$plot.linetypes <- scales::linetype_pal()(length(self$grouplist))
    names(self$plot.linetypes) <- self$grouplist
  }
  return(self$plot.linetypes)
}

#' Get color palette
#'
#' This function returns a list with colors for plotting
#'
#' @param object An ALASCA object
#' @return A list with colors
#'
#' @export
get_plot_palette <- function() {
  if (is.null(self$plot.palette)) {
    self$plot.palette <- scales::viridis_pal(end = self$plot.palette.end)(length(self$grouplist))
    names(self$plot.palette) <- self$grouplist
  } 
  return(self$plot.palette)
}

#' Get exploratory power for plot label
#'
#' This function returns ...
#'
#' @param object An ALASCA object
#' @param comp Which two components to plot (default: `c(1, 2`)
#' @return A ggplot2 objects.
.get_explained_label <- function(component = 1, effect = "time", type = "Score") {
  if (effect == "time") {
    paste0(type, " PC", component, " (", round(100 * self$ALASCA$loading$explained$time[component], 2), "%)")
  } else {
    paste0(type, " PC", component, " (", round(100 * self$ALASCA$loading$explained$group[component], 2), "%)")
  }
}



#' Get screeplot
#'
#' This function returns a screeplot for an ALASCA model showing what proportion of the variance each component of the model explains
#'
#' @param object An ALASCA object
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(12, 8, 300)`)
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return An ggplot2 object (or a list og ggplot objects)
#'
#' @examples
#' load("PE.Rdata")
#' screeplot(model)
#' @export
screeplot.AlascaModel <- function(object,
                             effect = "both",
                             nComps = NA,
                             filename = "scree_plot",
                             filetype = object$plot.filetype,
                             figsize = object$plot.figsize,
                             figunit = object$plot.figunit,
                             my_theme = object$plot.my_theme) {
  explained <- as.data.frame(object$get_scores()$explained)
  explained$component <- seq_len(nrow(explained))
  if (!is.na(filename)) object$filename <- filename
  
  if (any(!is.na(nComps))) {
    if (length(nComps) == 1) {
      explained <- subset(explained, component <= nComps)
    } else {
      explained <- subset(explained, component %in% nComps)
      explained$component <- factor(explained$component)
    }
  }
  g <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = time, group = NA)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    my_theme +
    ggplot2::labs(x = "Principal Component", y = paste0("Relative Expl. of ", object$plot.x_label, " Var."))
  if (object$separate_time_and_group) {
    gg <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = group, group = NA)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      my_theme +
      ggplot2::labs(x = "Principal Component", y = "Relative Expl. of Group Var.")
    if (effect == "both") {
      g <- ggpubr::ggarrange(g, gg)
    }
  }
  if (effect == "group") {
    if (object$save) {
      object$saveALASCAPlot(g = gg, filetype = filetype, figsize = figsize, figunit = figunit)
    }
    return(gg)
  } else {
    if (object$save) {
      object$saveALASCAPlot(g = g, filetype = filetype, figsize = figsize, figunit = figunit)
    }
    return(g)
  }
}


#' Get an ALASCA object
#'
#' This function plots your ALASCA model
#'
#' @param object An [ALASCA()] object
#' @param component Integer stating which component to return (1 is default)
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param decreasing_loadings Sort the loadings in decreasing (`TRUE`, default) or increasing order (`FALSE`)
#' @param only String stating which plot to return; `both` (default), `score` or `loading`
#' @param enlist Logical. If `TRUE`, the plots are returned as a list and not as a composed figure (default)
#' @param too_dense Integer, If > 0, only name this number of covariables
#' @param x_label Defaults to "Time" if not specified here or during model setup
#' @param group_label Defaults to "Group" if not specified here or  during model setup
#' @param flip_axes When `TRUE` (default), list the variable loadings vertical instead of horizontal
#' @param plot_zeroline When `TRUE` (default), plot a zero line in the loading plot
#' @param limit_loading Only list robust loadings
#' @param filetype Which file type you want to save the figure to (default: `png`)
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit Unit for figure size (default: `mm`)
#' @param highlight Vector of strings with variables to highlight in the loadings plot
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return An ggplot2 object (or a list og ggplot objects)
#'
#' @examples
#' load("PE.Rdata")
#' plot(model.val)
#' plot(model, component = "PC2")
#' plot(model, only = "score", effect = "time")
#' plot(model, too_dense = 5)
#' plot(model, highlight = c("PlGF", "IL-1b", "IL-6"))
#' @export
plot_development <- function(object,
                             component = 1,
                             effect = "both",
                             decreasing_loadings = TRUE,
                             only = "both",
                             enlist = FALSE,
                             too_dense = NA,
                             n_limit = ifelse(length(object$variablelist) < 40, 0, 20),
                             highlight = NA,
                             x_label = NA,
                             group_label = NA,
                             flip_axes = TRUE,
                             plot_zeroline = TRUE,
                             filename = NA,
                             filetype = NA,
                             figsize = NA,
                             variables = NA,
                             dodgewidth = 0.35,
                             figunit = NA,
                             plot_ribbon = TRUE,
                             loadinggroup = NA,
                             limit_loading = FALSE,
                             sort_by_loadinggroup = FALSE,
                             my_theme = object$plot.my_theme) {
  
  if (!(effect %in% c("both", "time", "group"))) {
    log4r::error(object$log, "`effect` has to be `both`, `time` or `group`")
    stop()
  }
  if (!object$separate_time_and_group) effect <- "time"
  if (!is.na(x_label)) object$plot.x_label <- x_label
  if (!is.na(filename)) object$filename <- filename
  if (!is.na(group_label)) object$plot.group_label <- group_label
  if (!is.na(filetype)) object$plot.filetype <- filetype
  if (!is.na(figsize)) object$plot.figsize <- figsize
  if (!is.na(figunit)) object$plot.figunit <- figunit
  if (n_limit > 0 && length(object$variablelist) > 2*n_limit) log4r::warn(object$log, paste0("Will only plot the ",n_limit," upper and lower loadings. Use `n_limit = 0` to show all"))
  
  if (flip_axes) {
    plotwidths <- c(2, 3, 2, 3)
    plotalign <- "h"
    decreasing_loadings <- !decreasing_loadings
  } else {
    plotwidths <- c(1, 2, 1, 2)
    plotalign <- "hv"
  }
  if (only == "score") {
    if (effect == "both") {
      g_score_time <- get_score_plot(object, component = component, effect = "time", plot_ribbon = plot_ribbon, dodgewidth = dodgewidth, my_theme = my_theme)
      g_score_group <- get_score_plot(object, component = component, effect = "group", plot_ribbon = plot_ribbon, dodgewidth = dodgewidth, my_theme = my_theme)
      g <- list(g_score_time, g_score_group)
    } else {
      g <- get_score_plot(object, component = component, effect = effect, plot_ribbon = plot_ribbon, dodgewidth = dodgewidth, my_theme = my_theme)
    }
  } else if (only == "loading") {
    if (effect == "both") {
      g_loading_time <- get_loading_plot(object,
                                         component = component,
                                         effect = "time",
                                         decreasing_loadings = decreasing_loadings,
                                         flip_axes = flip_axes,
                                         plot_zeroline = plot_zeroline,
                                         too_dense = too_dense,
                                         n_limit = n_limit,
                                         highlight = highlight,
                                         variables = variables,
                                         loadinggroup = loadinggroup,
                                         limit_loading = limit_loading,
                                         sort_by_loadinggroup = sort_by_loadinggroup,
                                         my_theme = my_theme
      )
      g_loading_group <- get_loading_plot(object,
                                          component = component,
                                          effect = "group",
                                          decreasing_loadings = decreasing_loadings,
                                          flip_axes = flip_axes,
                                          plot_zeroline = plot_zeroline,
                                          too_dense = too_dense,
                                          n_limit = n_limit,
                                          variables = variables,
                                          highlight = highlight,
                                          limit_loading = limit_loading,
                                          loadinggroup = loadinggroup,
                                          sort_by_loadinggroup = sort_by_loadinggroup,
                                          my_theme = my_theme
      )
      g <- list(g_loading_time, g_loading_group)
    } else {
      g <- get_loading_plot(object,
                            component = component,
                            effect = effect,
                            decreasing_loadings = decreasing_loadings,
                            flip_axes = flip_axes,
                            plot_zeroline = plot_zeroline,
                            too_dense = too_dense,
                            n_limit = n_limit,
                            variables = variables,
                            limit_loading = limit_loading,
                            loadinggroup = loadinggroup,
                            sort_by_loadinggroup = sort_by_loadinggroup,
                            my_theme = my_theme
      )
    }
  } else {
    if (effect == "both") {
      g_loading_time <- get_loading_plot(object,
                                         component = component,
                                         effect = "time",
                                         decreasing_loadings = decreasing_loadings,
                                         flip_axes = flip_axes,
                                         plot_zeroline = plot_zeroline,
                                         too_dense = too_dense,
                                         n_limit = n_limit,
                                         limit_loading = limit_loading,
                                         highlight = highlight,
                                         loadinggroup = loadinggroup,
                                         sort_by_loadinggroup = sort_by_loadinggroup,
                                         my_theme = my_theme
      )
      g_loading_group <- get_loading_plot(object,
                                          component = component,
                                          effect = "group",
                                          decreasing_loadings = decreasing_loadings,
                                          flip_axes = flip_axes,
                                          plot_zeroline = plot_zeroline,
                                          too_dense = too_dense,
                                          n_limit = n_limit,
                                          variables = variables,
                                          highlight = highlight,
                                          limit_loading = limit_loading,
                                          loadinggroup = loadinggroup,
                                          sort_by_loadinggroup = sort_by_loadinggroup,
                                          my_theme = my_theme
      )
      g_score_time <- get_score_plot(object, component = component, effect = "time", plot_ribbon = plot_ribbon, dodgewidth = dodgewidth, my_theme = my_theme)
      g_score_group <- get_score_plot(object, component = component, effect = "group", plot_ribbon = plot_ribbon, dodgewidth = dodgewidth, my_theme = my_theme)
      if (enlist) {
        g <- list(g_score_time, g_loading_time, g_score_group, g_loading_group)
      } else {
        if (is.na(loadinggroup) & is.na(object$plot.loadinggroupcolumn)) {
          g <- ggpubr::ggarrange(
            g_score_time,
            g_loading_time,
            g_score_group,
            g_loading_group,
            nrow = 2, ncol = 2,
            widths = plotwidths,
            align = plotalign,
            common.legend = TRUE,
            legend.grob = ggpubr::get_legend(g_score_group),
            legend = "bottom"
          )
        } else {
          g <- ggpubr::ggarrange(
            ggpubr::ggarrange(
              g_score_time + ggplot2::theme(legend.position = "none"),
              g_loading_time + ggplot2::theme(legend.position = "none"),
              align = plotalign,
              widths = plotwidths[1:2]
            ),
            ggpubr::ggarrange(
              g_score_group,
              g_loading_group,
              align = plotalign,
              widths = plotwidths[1:2]
            ),
            nrow = 2, ncol = 1,
            align = plotalign
          )
        }
      }
    } else {
      g_loading <- get_loading_plot(object,
                                    component = component,
                                    effect = effect,
                                    decreasing_loadings = decreasing_loadings,
                                    flip_axes = flip_axes,
                                    plot_zeroline = plot_zeroline,
                                    too_dense = too_dense,
                                    n_limit = n_limit,
                                    variables = variables,
                                    limit_loading = limit_loading,
                                    highlight = highlight,
                                    loadinggroup = loadinggroup,
                                    sort_by_loadinggroup = sort_by_loadinggroup,
                                    my_theme = my_theme
      )
      g_score <- get_score_plot(object, component = component, effect = effect, plot_ribbon = plot_ribbon, dodgewidth = dodgewidth, my_theme = my_theme)
      if (enlist) {
        g <- list(g_score, g_loading)
      } else {
        if (is.na(loadinggroup) & is.na(object$plot.loadinggroupcolumn)) {
          g <- ggpubr::ggarrange(
            g_score,
            g_loading,
            nrow = 1,
            ncol = 2,
            widths = plotwidths,
            align = plotalign,
            common.legend = TRUE,
            legend = "bottom"
          )
        } else {
          g <- ggpubr::ggarrange(g_score,
                                 g_loading,
                                 nrow = 1,
                                 ncol = 2,
                                 widths = plotwidths,
                                 align = plotalign
          )
        }
      }
    }
  }
  if (object$save) {
    saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  }
  return(g)
}



#' Get loadings
#'
#' This function  returns the loadings for an ALASCA model
#'
#' @param object An ALASCA object
#' @param n_limit Returns the n highest and lowest loadings by PC (i.e., 2*n_limit loadings per PC)
#' @return A list with loadings for time (and group), and the exploratory power for each component
#' @export
get_loadings <- function(object, limit_loading = FALSE, n_limit = 0L, component = c(0)) {
  dfl <- list()
  if (component[[1]] > 0 || length(component) > 1) {
    dfl$time <- object$ALASCA$loading$time[PC %in% component]
    if (object$separate_time_and_group) {
      dfl$group <- object$ALASCA$loading$group[PC %in% component]
    }
  } else {
    dfl$time <- object$ALASCA$loading$time
    if (object$separate_time_and_group) {
      dfl$group <- object$ALASCA$loading$group
    }
  }
  if (limit_loading && object$validate) {
    dfl$time <- dfl$time[!is.na(low) & sign(low) == sign(high)]
    if (object$separate_time_and_group) {
      dfl$group <- dfl$group[!is.na(low) & sign(low) == sign(high)]
    }
  } else {
    if (n_limit > 0L) {
      index_head_and_tail <- c(seq(n_limit), length(object$variablelist)+1-seq(n_limit))
      dfl$time <- dfl$time[dfl$time[order(loading, decreasing = TRUE), .I[index_head_and_tail], by = PC]$V1]
      if (object$separate_time_and_group) {
        dfl$group <- dfl$group[dfl$group[order(loading, decreasing = TRUE), .I[index_head_and_tail], by = PC]$V1]
      }
    }
  }
  
  return(dfl)
}

#' Get scores
#'
#' This function returns the scores for an ALASCA model
#'
#' @inheritParams get_loadings
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
get_scores <- function(component = 0) {
  if(any(component > 0)){
    self$ALASCA$score$time <- self$ALASCA$score$time[ PC %in% component]
    if (self$separate_time_and_group) {
      self$ALASCA$score$group <- self$ALASCA$score$group[ PC %in% component]
    }
  }
  return(self$ALASCA$score)
}



#' Get loading plot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An ALASCA object
#' @param component Which component to plot?
#' @param effect Plot time or group
#' @param decreasing_loadings Logical. Should loadings be sorted in decreasing order?
#' @param flip_axes When `TRUE` (default), list the variable loadings vertical instead of horizontal
#' @param plot_zeroline When `TRUE` (default), plot a zero line in the loading plot
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(12, 8, 300)`)
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot object
get_loading_plot <- function(object,
                             component = 1,
                             effect = "time",
                             decreasing_loadings = TRUE,
                             too_dense = NA,
                             n_limit = 0,
                             highlight = NA,
                             filetype = NA,
                             flip_axes = TRUE,
                             plot_zeroline = TRUE,
                             figsize = NA,
                             figunit = NA,
                             variables = NA,
                             loadinggroup = NA,
                             limit_loading = FALSE,
                             sort_by_loadinggroup = TRUE,
                             point_size = 0.4,
                             my_theme = NA) {
  
  if (any(is.na(my_theme))) my_theme <- object$plot.my_theme
  if (!is.na(filetype)) object$plot.filetype <- filetype
  if (any(is.na(variables))) variables <- object$variablelist
  if (!is.na(figsize)) object$plot.figsize <- figsize
  if (!is.na(figunit)) object$plot.figunit <- figunit
  
  if (effect == "time") {
    loadings <- subset(get_loadings(object, limit_loading = limit_loading, n_limit = n_limit, component = component)$time, covars %in% variables)
  } else {
    loadings <- subset(get_loadings(object, limit_loading = limit_loading, n_limit = n_limit, component = component)$group, covars %in% variables)
  }
  if (!is.na(object$plot.loadinggroupcolumn)) {
    loadings <- merge(loadings, object$variable_labels, by = "covars")
  } else {
    loadings$covargroup <- NA
  }
  if (sort_by_loadinggroup & !is.na(object$plot.loadinggroupcolumn)) {
    loadings$covars <- factor(loadings$covars, levels = unique(loadings$covars[order(loadings$covargroup, loadings$loading, decreasing = decreasing_loadings)]))
  } else {
    loadings$covars <- factor(loadings$covars, levels = unique(loadings$covars[order(loadings$loading, decreasing = decreasing_loadings)]))
  }
  if (object$validate) {
    if (any(colnames(loadings) == "model")) {
      g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high, shape = model)) +
        ggplot2::geom_pointrange(size = point_size)
    } else {
      if (is.na(object$plot.loadinggroupcolumn)) {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high)) +
          ggplot2::geom_pointrange(size = point_size)
      } else {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high, color = covargroup, shape = covargroup)) +
          ggplot2::geom_pointrange(size = point_size)
      }
    }
  } else {
    if (any(colnames(loadings) == "model")) {
      g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, shape = model)) +
        ggplot2::geom_point()
    } else {
      if (is.na(object$plot.loadinggroupcolumn)) {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading)) +
          ggplot2::geom_point()
      } else {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, color = covargroup, shape = covargroup)) +
          ggplot2::geom_point()
      }
    }
  }
  g <- g + my_theme +
    ggplot2::labs(
      x = "Variable",
      y = .get_explained_label(object, component = component, effect = effect, type = "Loading")
    )
  if (plot_zeroline) {
    g <- g + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  }
  if (n_limit > 0 & !sort_by_loadinggroup) {
    g <- g + ggplot2::geom_vline(xintercept = n_limit + 0.5, linetype = "dotted")
  }
  if (flip_axes) {
    g <- g + ggplot2::coord_flip()
  } else {
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = c(0.8, 0.8))
  }
  if (!is.na(object$plot.loadinggroupcolumn)) {
    g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
      ggplot2::labs(color = object$plot.loadinggroup_label, shape = object$plot.loadinggroup_label) +
      ggplot2::theme(legend.position = "bottom") # ggplot2::scale_color_brewer(palette = "Dark2")
  }
  if (!any(is.na(highlight))) {
    g <- g + ggplot2::geom_point(color = ifelse(loadings$covars %in% highlight, "red", "grey")) +
      ggrepel::geom_text_repel(data = subset(loadings, covars %in% highlight), ggplot2::aes(label = covars), max.iter = 5000) +
      my_theme +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )
  } else if (!is.na(too_dense) & too_dense > 0) {
    limUpper <- unique(loadings$loading[order(loadings$loading, decreasing = TRUE)[too_dense]])
    limLower <- unique(loadings$loading[order(loadings$loading, decreasing = FALSE)[too_dense]])
    g <- g + ggplot2::geom_point(color = ifelse(loadings$loading <= limLower | loadings$loading >= limUpper, "red", "grey")) +
      ggrepel::geom_text_repel(
        data = subset(loadings, loading <= limLower | loading >= limUpper),
        ggplot2::aes(label = covars),
        max.iter = 5000
      ) +
      my_theme
    if (flip_axes) {
      g <- g + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      )
    } else {
      g <- g + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )
    }
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit, suffix = "_loading")
  
  return(g)
}

#' Get score plot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An ALASCA object
#' @param component Which component to plot?
#' @param effect Plot time or group
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::scale_color_viridis_d(end = 0.9) + ggplot2::theme_bw()`
#' @return A ggplot object
get_score_plot <- function(object,
                           component = 1,
                           effect = "time",
                           filetype = NA,
                           figsize = NA,
                           figunit = NA,
                           plot_ribbon = TRUE,
                           dodgewidth = 0.35,
                           point_size = 2,
                           my_theme = NA) {
  if (any(is.na(my_theme))) my_theme <- object$plot.my_theme
  if (!is.na(filetype)) object$plot.filetype <- filetype
  if (any(!is.na(figsize))) object$plot.figsize <- figsize
  if (!is.na(figunit)) object$plot.figunit <- figunit
  
  if (effect == "time") {
    if (object$separate_time_and_group) {
      score <- get_scores(object, component = component)$time
      if (object$validate) {
        # Show error bars
        if (grepl("permutation", object$validation_method)) {
          pvals <- object$pvals
          score <- merge(score, pvals, by.x = "time", by.y = "effect", all.x = TRUE, all.y = FALSE)
          score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, label = p.value.str)) +
            ggplot2::geom_point(position = ggplot2::position_dodge(width = dodgewidth)) +
            ggplot2::geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = dodgewidth), show.legend = FALSE) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
        } else {
          if (any(colnames(score) == "model")) {
            g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = model, linetype = model, shape = model, ymin = low, ymax = high)) +
              ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
              ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
          } else {
            g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, ymin = low, ymax = high)) +
              ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth))
            if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
            if (plot_ribbon && object$method %in% c("LMM")) {
              g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                                            alpha = .1,
                                            position = ggplot2::position_dodge(width = dodgewidth), color = NA
              ) +
                ggplot2::scale_fill_manual(values = get_plot_palette(object)) +
                ggplot2::labs(fill = object$plot.group_label)
            }
          }
        }
      } else {
        # No validation - no error bars
        if (any(colnames(score) == "model")) {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = model, linetype = model)) +
            ggplot2::geom_point() +
            ggplot2::geom_line()
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group)) +
            ggplot2::geom_point()
          if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line()
        }
      }
    } else {
      score <- get_scores(object, component = component)$time
      if (object$validate) {
        if (grepl("permutation", object$validation_method)) {
          pvals <- object$pvals
          score$effect <- paste(score$time, score$group)
          score <- merge(score, pvals, by.x = "effect", by.y = "effect", all.x = TRUE, all.y = FALSE)
          score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, label = p.value.str)) +
            ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
            geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth))
          if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
          if (plot_ribbon && object$method %in% c("LMM")) {
            g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                                          alpha = .1,
                                          position = ggplot2::position_dodge(width = dodgewidth), color = NA
            ) +
              ggplot2::scale_fill_manual(values = get_plot_palette(object)) + ggplot2::labs(fill = object$plot.group_label)
          }
        }
      } else {
        if (any(colnames(score) == "model")) {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = model)) +
            ggplot2::geom_point() +
            ggplot2::geom_line()
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group)) +
            ggplot2::geom_point()
          if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line()
        }
      }
    }
    g <- g + my_theme +
      ggplot2::scale_color_manual(values = get_plot_palette(object))
    if (object$method %in% c("LMM")) g <- g + ggplot2::scale_linetype_manual(values = get_plot_linetypes(object))
    g <- g + ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        x = object$plot.x_label,
        group = object$plot.group_label, color = object$plot.group_label, linetype = object$plot.group_label,
        y = .get_explained_label(object, component = component, effect = "time")
      )
  } else {
    # Group effect
    score <- get_scores(object, component = component)$group
    if (object$validate) {
      if (grepl("permutation", object$validation_method)) {
        pvals <- object$pvals
        score$effect <- paste(score$time, score$group)
        score <- merge(score, pvals, by.x = "effect", by.y = "effect", all.x = TRUE, all.y = FALSE)
        score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, label = p.value.str)) +
          ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
          geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
      } else {
        if (any(colnames(score) == "model")) {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, linetype = model, shape = model, color = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
          if (plot_ribbon) {
            g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                                          alpha = .1,
                                          position = ggplot2::position_dodge(width = dodgewidth), color = NA
            ) +
              ggplot2::scale_fill_manual(values = get_plot_palette(object)) +
              ggplot2::labs(fill = object$plot.group_label)
          }
        }
      }
    } else {
      if (any(colnames(score) == "model")) {
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = model)) +
          ggplot2::geom_point() +
          ggplot2::geom_line()
      } else {
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group)) +
          ggplot2::geom_point() +
          ggplot2::geom_line()
      }
    }
    g <- g + my_theme +
      ggplot2::scale_color_manual(values = get_plot_palette(object)) +
      ggplot2::scale_linetype_manual(values = get_plot_linetypes(object)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        x = object$plot.x_label,
        group = object$plot.group_label, color = object$plot.group_label, linetype = object$plot.group_label,
        y = .get_explained_label(object, component = component, effect = "group")
      )
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit, suffix = "_score")
  
  return(g)
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
        log4r::error(object$log, "You need to specify participant and value columns")
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
          log4r::error(object$log, "You need to specify participant column")
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
      log4r::error(object$log, "Wrong input object: must be a ALASCA model or a data frame")
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

#' Plot model predictions
#'
#' This function returns the scores for an ALASCA model
#'
#' @param object An ALASCA object or a data frame. If a data frame, you need to specify the column names for participant and value. This also applies if you have not specified the participant column in the ALASCA model before.
#' @param variable List of variable names to print. If `NA`, return all (default).
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm",
#' @return A list with ggplot2 objects.
#'
#' @examples
#' load("PE.Rdata")
#' plotPred(model, variable = "IL-6")[[1]] + ggplot2::theme_bw()
#' do.call(
#'   ggpubr::ggarrange,
#'   c(plotPred(model, variable = c("PlGF", "IL-6", "IL-1b", "IFN-g", "Eotaxin-2", "Eotaxin")),
#'     common.legend = TRUE, legend = "bottom"
#'   )
#' )
#' @export
plot_prediction <- function(object,
                            variables = object$variablelist,
                            filetype = NA,
                            figsize = NA,
                            figunit = NA,
                            filename = NA,
                            dodgewidth = 0.35,
                            plot_ribbon = TRUE,
                            as.list = FALSE,
                            plot.ylabel = "value",
                            my_theme = object$plot.my_theme) {
  
  if (!is.na(filename)) object$filename <- filename
  
  if ( as.list ){
    if (object$validate) {
      gg <- lapply(variables, function(x) {
        g <- ggplot2::ggplot(subset(object$model_prediction, variable == x), ggplot2::aes(
          x = time,
          y = pred,
          color = group,
          group = group,
          linetype = group,
          ymin = low,
          ymax = high
        )) +
          ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth)) +
          ggplot2::scale_color_manual(values = get_plot_palette(object)) +
          my_theme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = object$plot.x_label, y = plot.ylabel, color = object$plot.group_label, linetype = object$plot.group_label)
        if (plot_ribbon) {
          g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                                        alpha = .1,
                                        position = ggplot2::position_dodge(width = dodgewidth), color = NA
          ) +
            ggplot2::scale_fill_manual(values = get_plot_palette(object)) + ggplot2::labs(fill = object$plot.group_label)
        }
        g
      })
    } else {
      gg <- lapply(variables, function(x) {
        g <- ggplot2::ggplot(subset(object$model_prediction, variable == x), ggplot2::aes(x = time, y = pred, color = group, linetype = group, group = group)) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_color_manual(values = get_plot_palette(object)) +
          my_theme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = object$plot.x_label, y = plot.ylabel, color = object$plot.group_label, linetype = object$plot.group_label)
        g
      })
    }
    names(gg) <- variable
    if (object$save) {
      for (i in seq_along(gg)) {
        saveALASCAPlot(object = object, g = gg[[i]], filetype = filetype, figsize = figsize, figunit = figunit, suffix = names(gg)[i])
      }
    }
  } else {
    if (object$validate) {
      g <- ggplot2::ggplot(object$model_prediction[object$model_prediction$variable %in% variables,], ggplot2::aes(
        x = time,
        y = pred,
        color = group,
        group = group,
        linetype = group,
        ymin = low,
        ymax = high
      )) +
        ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
        ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth)) +
        ggplot2::scale_color_manual(values = get_plot_palette(object)) +
        my_theme +
        ggplot2::facet_wrap(~variable, scales = "free_y") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = object$plot.x_label, y = plot.ylabel, color = object$plot.group_label, linetype = object$plot.group_label)
      if (plot_ribbon) {
        g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                                      alpha = .1,
                                      position = ggplot2::position_dodge(width = dodgewidth), color = NA
        ) +
          ggplot2::scale_fill_manual(values = get_plot_palette(object)) + ggplot2::labs(fill = object$plot.group_label)
      }
    } else {
      g <- ggplot2::ggplot(object$model_prediction[object$model_prediction$variable %in% variables,], ggplot2::aes(x = time, y = pred, color = group, linetype = group, group = group)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(values = get_plot_palette(object)) +
        my_theme +
        ggplot2::facet_wrap(~variable, scales = "free_y") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = object$plot.x_label, y = plot.ylabel, color = object$plot.group_label, linetype = object$plot.group_label)
    }
    if (object$save) {
      saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
    }
    return(g)
  }
  
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

#' Plot validations models
#'
#' This function returns a plot of the validation models
#'
#' @param object A validated ALASCA object
#' @param component Which component to plot (default: 1)
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param my_theme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#'
#' @export
plot_validation <- function(object,
                            component = 1,
                            filetype = NA,
                            filename = NA,
                            figsize = NA,
                            figunit = NA,
                            n_limit = 0,
                            decreasing_loadings = FALSE,
                            plot_zeroline = TRUE,
                            my_theme = NA,
                            flip_axes = TRUE,
                            plot.alpha = 0.3) {
  if (!object$validate) {
    log4r::error(object$log, message = "You must validate the model first!")
    stop()
  }
  if (any(is.na(my_theme)))  my_theme <- object$plot.my_theme
  if (!is.na(filename)) object$filename <- filename
  
  if (object$separate_time_and_group) {
    # Score plots
    ## Time
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    }))
    dfm <- get_scores(object, component = component)$time
    gst <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = model, color = group, linetype = group)) +
      ggplot2::geom_point(alpha = plot.alpha)
    if (object$method %in% c("LMM")) gst <- gst + ggplot2::geom_line(alpha = plot.alpha)
    gst <- gst + ggplot2::geom_point(data = dfm, group = NA, alpha = 1, color = "black")
    if (object$method %in% c("LMM")) gst <- gst + ggplot2::geom_line(data = dfm, group = dfm$group, alpha = 1, color = "black")
    gst <- gst + ggplot2::labs(x = object$plot.x_label,
                               color = object$plot.group_label, linetype = object$plot.group_label,
                               y = .get_explained_label(object, component = component, effect = "time", type = "Score")) +
      ggplot2::scale_color_manual(values = get_plot_palette(object)) +
      ggplot2::scale_linetype_manual(values = get_plot_linetypes(object)) +
      my_theme + ggplot2::theme(legend.position = "none")
    
    ## Group
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$group, PC == component),
        model = x
      )
    }))
    dfm <- get_scores(object, component = component)$group
    dff$plotGroup <- paste0(dff$model, "-", dff$group)
    gsg <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = plotGroup, color = group, linetype = group)) +
      ggplot2::geom_point(alpha = plot.alpha)
    if (object$method %in% c("LMM")) gsg <- gsg + ggplot2::geom_line(alpha = plot.alpha)
    gsg <- gsg + ggplot2::geom_point(data = dfm, group = NA, alpha = 1, color = "black")
    if (object$method %in% c("LMM")) gsg <- gsg + ggplot2::geom_line(data = dfm, group = dfm$group, alpha = 1, color = "black")
    gsg <- gsg +   ggplot2::scale_color_manual(values = get_plot_palette(object)) +
      ggplot2::scale_linetype_manual(values = get_plot_linetypes(object)) +
      ggplot2::labs(x = object$plot.x_label,
                    y = .get_explained_label(object, component = component, effect = "group", type = "Score")) +
      my_theme + ggplot2::theme(legend.position = "bottom")
    
    # Loading plot
    ## Time
    dfm <- get_loadings(object, component = component, n_limit = n_limit)$time
    dfm$covars <- factor(dfm$covars, levels = unique(dfm$covars[order(dfm$loading, decreasing = decreasing_loadings)]))
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_loadings(object$validation$temp_objects[[x]], component = component)$time, covars %in% dfm$covars)
      )
    }))
    dff$covars <- factor(dff$covars, levels = unique(dff$covars[order(dff$loading, decreasing = decreasing_loadings)]))
    
    glt <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) +
      ggplot2::geom_point(data = dfm, alpha = 1, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(ymax = high, ymin = low), alpha = 0.3, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(xmax = as.numeric(covars) + 0.5, xmin = as.numeric(covars) - 0.5), alpha = 1, color = "black") +
      ggplot2::geom_point(alpha = 0.2, color = "black")
    if (plot_zeroline) glt <- glt + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    if (flip_axes) glt <- glt + ggplot2::coord_flip()
    glt <- glt + ggplot2::labs(x = "Variable",
                               y = .get_explained_label(object, component = component, effect = "time", type = "Loading")) +
      ggplot2::theme(legend.position = "none") +
      my_theme
    
    ## Group
    dfm <- get_loadings(object, component = component, n_limit = n_limit)$group
    dfm$covars <- factor(dfm$covars, levels = unique(dfm$covars[order(dfm$loading, decreasing = decreasing_loadings)]))
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_loadings(object$validation$temp_objects[[x]], component = component)$group, covars %in% dfm$covars)
      )
    }))
    dff$covars <- factor(dff$covars, levels = unique(dff$covars[order(dff$loading, decreasing = decreasing_loadings)]))
    glg <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) +
      ggplot2::geom_point(data = dfm, alpha = 1, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(ymax = high, ymin = low), alpha = 0.3, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(xmax = as.numeric(covars) + 0.5, xmin = as.numeric(covars) - 0.5), alpha = 1, color = "black") +
      ggplot2::geom_point(alpha = 0.2, color = "black")
    if (plot_zeroline) glg <- glg + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    if (flip_axes) glg <- glg + ggplot2::coord_flip()
    glg <- glg + ggplot2::labs(x = "Variable",
                               y = .get_explained_label(object, component = component, effect = "group", type = "Loading")) +
      ggplot2::theme(legend.position = "bottom") +
      my_theme
    
    g <- ggpubr::ggarrange(gst, glt, gsg, glg, nrow = 2, ncol = 2, widths = c(2, 3, 2, 3), common.legend = TRUE, legend = "bottom")
  } else {
    # Score plot
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    }))
    dff$plotGroup <- paste0(dff$model, "-", dff$group)
    dfm <- get_scores(object, component = component)$time
    gs <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = plotGroup, color = group, linetype = group)) +
      ggplot2::geom_point(alpha = plot.alpha)
    if (object$method %in% c("LMM")) gs <- gs + ggplot2::geom_line(alpha = plot.alpha)
    gs <- gs + ggplot2::geom_point(data = dfm, group = NA, alpha = 1, color = "black")
    if (object$method %in% c("LMM")) gs <- gs + ggplot2::geom_line(data = dfm, group = dfm$group, alpha = 1, color = "black")
    gs <- gs + ggplot2::scale_color_manual(values = get_plot_palette(object)) +
      ggplot2::scale_linetype_manual(values = get_plot_linetypes(object)) +
      ggplot2::labs(x = object$plot.x_label,
                    color = object$plot.group_label, linetype = object$plot.group_label,
                    y = .get_explained_label(object, component = component, effect = "time", type = "Score")) +
      my_theme +
      ggplot2::theme(legend.position = "bottom")
    
    # Loading plot
    dfm <- get_loadings(object, component = component, n_limit = n_limit)$time
    dfm$covars <- factor(dfm$covars, levels = unique(dfm$covars[order(dfm$loading, decreasing = decreasing_loadings)]))
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_loadings(object$validation$temp_objects[[x]], component = component)$time, covars %in% dfm$covars)
      )
    }))
    dff$covars <- factor(dff$covars, levels = unique(dff$covars[order(dff$loading, decreasing = decreasing_loadings)]))
    
    gl <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) +
      ggplot2::geom_point(data = dfm, alpha = 1, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(ymax = high, ymin = low), alpha = 0.3, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(xmax = as.numeric(covars) + 0.5, xmin = as.numeric(covars) - 0.5), alpha = 1, color = "black") +
      ggplot2::geom_point(alpha = 0.2, color = "black")
    if (plot_zeroline) gl <- gl + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    if (flip_axes) gl <- gl + ggplot2::coord_flip()
    gl <- gl + ggplot2::labs(x = "Variable",
                             color = object$plot.group_label,
                             y = .get_explained_label(object, component = component, effect = "time", type = "Loading")) +
      my_theme
    
    g <- ggpubr::ggarrange(gs, gl, nrow = 1, ncol = 2, widths = c(2, 3), common.legend = TRUE, legend = "bottom")
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  
  return(g)
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

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param x_label Alternative names for the covariables
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm",
#' @param return_data Set to `TRUE` to return data instead of plot
#'
#' @param my_theme A ggplot2 theme
#' @return A ggplot2 objects\.
#'
#' @export
plot_covars <- function(object,
                        covar = NA,
                        x_label = NA,
                        variables = NA,
                        n_limit = 0,
                        filename = "covars",
                        return_data = FALSE,
                        filetype = NA,
                        figsize = NA,
                        figunit = NA,
                        my_theme = NA,
                        pvalue = "star") {
  if (any(is.na(my_theme))) my_theme <- object$plot.my_theme
  if (!is.na(filename)) object$filename <- filename
  
  df <- get_covars(object, n_limit = n_limit)
  if ( nrow(df) == 0 ) {
    log4r::error(object$log, "No covariates to plot")
    stop()
  }
  df$covar <- factor(df$covar, levels = unique(df$covar[order(df$estimate)]))
  if (any(is.na(covar))) {
    covar <- unique(df$variable)
  }
  if (any(is.na(variables))) {
    variables <- object$variablelist
  }
  df <- subset(df, covar %in% variables)
  if (any(is.na(x_label))) {
    x_label <- covar
  }
  for (i in seq_len(length(covar))) {
    df$x_label[df$variable == covar[i]] <- x_label[i]
  }
  if (!is.na(object$plot.loadinggroupcolumn)) {
    df <- merge(df, object$variable_labels, by.x = "covar", by.y = "covars")
  } else {
    df$covargroup <- NA
  }
  if (!object$use_Rfast) {
    # lmer and lm provide p values
    df$pvalue_label <- ifelse(df$pvalue >= 0.05, "Not significant", ifelse(df$pvalue < 0.001, "< 0.001", ifelse(df$pvalue < 0.01, "< 0.01", "< 0.05")))
    df$pvalue_sign <- ifelse(df$pvalue >= 0.05, "", ifelse(df$pvalue < 0.001, "***", ifelse(df$pvalue < 0.01, "**", "*")))
    if (return_data) {
      return(df)
    } else {
      if (pvalue == "shape") {
        object$plot.loadinggroupcolumn <- NA # Cannot use shape for both loadinggroup and significance
        g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, shape = pvalue_label)) +
          ggplot2::scale_shape_manual(values = c("Not significant" = 3, "< 0.05" = 15, "< 0.01" = 16, "< 0.001" = 17, "Baseline" = 5))
      } else if (pvalue == "star" | pvalue == "asterisk" | pvalue == "stars") {
        if (is.na(object$plot.loadinggroupcolumn)) {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, label = pvalue_sign)) +
            ggplot2::geom_text(hjust = 0.5, vjust = 0, show.legend = FALSE)
        } else {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, label = pvalue_sign, color = covargroup, shape = covargroup)) +
            ggplot2::geom_text(hjust = 0.5, vjust = 0, show.legend = FALSE)
        }
      } else {
        g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar))
      }
      g <- g +
        ggplot2::geom_point() +
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::facet_wrap(~x_label, scales = "free_y") + ggplot2::labs(x = "Coefficient", y = "", shape = "P value") +
        my_theme + ggplot2::theme(legend.position = "bottom", legend.box = "vertical", legend.margin = ggplot2::margin())
      
      if (!is.na(object$plot.loadinggroupcolumn)) {
        g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
          ggplot2::labs(color = object$plot.loadinggroup_label, shape = object$plot.loadinggroup_label) +
          ggplot2::theme(legend.position = "bottom") # ggplot2::scale_color_brewer(palette = "Dark2")
      }
      
      if (object$save) {
        saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
      }
      return(g)
    }
  } else {
    # Rfast does not provide p values, use bootstrap intervals
    df$covar <- factor(df$covar, levels = unique(df$covar[order(df$estimate)]))
    if (return_data) {
      return(df)
    } else {
      if (object$validate) {
        if (is.na(object$plot.loadinggroupcolumn)) {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, xmin = low, xmax = high)) +
            ggplot2::geom_pointrange()
        } else {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, xmin = low, xmax = high, color = covargroup, shape = covargroup)) +
            ggplot2::geom_pointrange()
        }
      } else {
        if (is.na(object$plot.loadinggroupcolumn)) {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar)) +
            ggplot2::geom_point()
        } else {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, color = covargroup, shape = covargroup)) +
            ggplot2::geom_point()
        }
      }
      
      g <- g + 
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::facet_wrap(~x_label, scales = "free_y") +
        ggplot2::labs(x = "Coefficient", y = "") +
        my_theme
      
      if (!is.na(object$plot.loadinggroupcolumn)) {
        g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
          ggplot2::labs(color = object$plot.loadinggroup_label, shape = object$plot.loadinggroup_label) +
          ggplot2::theme(legend.position = "bottom") # ggplot2::scale_color_brewer(palette = "Dark2")
      }
      
      if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
      
      return(g)
    }
  }
}

#' Plot PCs
#'
#' This function returns a plot of ...
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param xlab Alternative names for the covariables
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm"
#' @param plottext If `TRUE`, plot time point as text
#' @param validationshape Either `NA`, "ellipse", or "cross"
#' @param validationlevel Defaults to 0.95
#' @param return_data Set to `TRUE` to return data instead of plot
#'
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot2 object.
#'
#' @export
plot_components <- function(object,
                            comps = c(1, 2),
                            filename = NA,
                            filetype = NA,
                            figsize = NA,
                            figunit = NA,
                            ...) {
  g_score <- plot_components_score(object = object, comps = comps, ...)
  g_loading <- plot_components_loadings(object = object, comps = comps, ...)
  g <- ggpubr::ggarrange(
    g_score, g_loading,
    ncol = 1, nrow = 2
  )
  if (!is.na(filename)) object$filename <- filename
  if (object$save) {
    saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  }
  return(g)
}

plot_components_loadings <- function(object,
                                     comps = c(1, 2),
                                     filetype = NA,
                                     figsize = NA,
                                     figunit = NA,
                                     validationshape = NA,
                                     myTheme = ggplot2::theme_classic(),
                                     ...) {
  if (any(!comps %in% get_relevant_pcs(object = object, effect = "time"))) {
    warning("Please note: Some components have low explanatory power and HAVE NOT BEEN rotated during rotation. Proceed with care.")
  }
  
  dff <- subset(get_loadings(object)$time, PC %in% comps)
  dff$PC <- paste0("PC", dff$PC)
  dff <- reshape2::melt(data = dff, id.vars = c("PC", "covars"))
  dff <- reshape2::dcast(data = dff, covars ~ PC + variable, value.var = "value")
  dff$color <- "black"
  dff$lab <- dff$covars
  if (object$validate == TRUE) {
    dff$color <- ifelse(sign(dff[paste0("PC", comps[1], "_low")]) == sign(dff[paste0("PC", comps[1], "_high")]) & sign(dff[paste0("PC", comps[2], "_low")]) == sign(dff[paste0("PC", comps[2], "_high")]), "black", "gray")
    dff$lab <- ifelse(dff$color == "black", dff$covars, NA)
  }
  
  
  g <- ggplot2::ggplot(dff, ggplot2::aes_string(
    x = paste0("PC", comps[1], "_loading"),
    y = paste0("PC", comps[2], "_loading"),
    label = "lab"
  ))
  if (object$validate == TRUE) {
    g <- g +
      ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high")), color = dff$color) +
      ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")), color = dff$color)
  } else {
    g <- g + ggplot2::geom_point()
  }
  g <- g +
    ggrepel::geom_text_repel() +
    ggplot2::labs(
      x = .get_explained_label(object, component = comps[1], effect = "time", type = "Loading"),
      y = .get_explained_label(object, component = comps[2], effect = "time", type = "Loading")
    ) +
    myTheme
  
  if (object$separate_time_and_group) {
    dff <- subset(get_loadings(object)$group, PC %in% comps)
    dff$PC <- paste0("PC", dff$PC)
    dff <- reshape2::melt(data = dff, id.vars = c("PC", "covars"))
    dff <- reshape2::dcast(data = dff, covars ~ PC + variable, value.var = "value")
    dff$color <- "black"
    dff$lab <- dff$covars
    if (object$validate == TRUE) {
      dff$color <- ifelse(sign(dff[paste0("PC", comps[1], "_low")]) == sign(dff[paste0("PC", comps[1], "_high")]) & sign(dff[paste0("PC", comps[2], "_low")]) == sign(dff[paste0("PC", comps[2], "_high")]), "black", "gray")
      dff$lab <- ifelse(dff$color == "black", dff$covars, NA)
    }
    
    glg <- ggplot2::ggplot(dff, ggplot2::aes_string(
      x = paste0("PC", comps[1], "_loading"),
      y = paste0("PC", comps[2], "_loading"),
      label = "lab"
    ))
    if (object$validate == TRUE) {
      glg <- glg +
        ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high")), color = dff$color) +
        ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")), color = dff$color)
    } else {
      glg <- glg + ggplot2::geom_point()
    }
    glg <- glg +
      ggrepel::geom_text_repel() +
      ggplot2::labs(
        x = .get_explained_label(object, component = comps[1], effect = "group", type = "Loading"),
        y = .get_explained_label(object, component = comps[2], effect = "group", type = "Loading")
      ) +
      myTheme
    
    g <- ggpubr::ggarrange(g, glg)
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  
  return(g)
}

#' Plot scores along PCs
#'
#' This function returns a plot of ...
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param xlab Alternative names for the covariables
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm"
#' @param plottext If `TRUE`, plot time point as text
#' @param validationshape Either `NA`, "ellipse", or "cross"
#' @param validationlevel Defaults to 0.95
#' @param return_data Set to `TRUE` to return data instead of plot
#'
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot2 object.
#'
#' @export
plot_components_score <- function(object,
                                  comps = c(1, 2),
                                  xlabel = NA,
                                  return_data = FALSE,
                                  filename = NA,
                                  filetype = NA,
                                  figsize = NA,
                                  figunit = NA,
                                  validationshape = "ellipse",
                                  validationlevel = 0.95,
                                  plottext = TRUE,
                                  texthjust = -0.2,
                                  alphavalidate = 0.4,
                                  myTheme = ggplot2::theme_classic(),
                                  ...) {
  
  if (!is.na(filename)) object$filename <- filename
  if (!object$validate) validationshape <- NA
  if (any(!comps %in% get_relevant_pcs(object = object, effect = "time"))) warning("Please note: Some components have low explanatory power and HAVE NOT BEEN rotated during rotation. Proceed with care.")
  
  if (validationshape == "cross" & !is.na(validationshape)) {
    dff <- subset(get_scores(object)$time, PC %in% comps)
    dff$PC <- paste0("PC", dff$PC)
    dff <- reshape2::melt(data = dff, id.vars = c("PC", "time", "group"))
    dff <- reshape2::dcast(data = dff, time + group ~ PC + variable, value.var = "value")
    dff$time <- factor(dff$time, levels = object$timelist)
    dff$group <- factor(dff$group, levels = object$grouplist)
    g <- ggplot2::ggplot(dff, ggplot2::aes_string(
      x = paste0("PC", comps[1], "_score"),
      y = paste0("PC", comps[2], "_score"),
      shape = "time",
      color = "group",
      group = "group",
      linetype = "group"
    )) +
      ggplot2::geom_line() +
      ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high"))) +
      ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")))
    g <- .get_plot_handle(
      g = g,
      object = object,
      myTheme = myTheme,
      comp = comps,
      effect = "time"
    )
    
    if (object$separate_time_and_group) {
      dff <- subset(get_scores(object)$group, PC %in% comps)
      dff$PC <- paste0("PC", dff$PC)
      dff <- reshape2::melt(data = dff, id.vars = c("PC", "time", "group"))
      dff <- reshape2::dcast(data = dff, time + group ~ PC + variable, value.var = "value")
      dff$time <- factor(dff$time, levels = object$timelist)
      dff$group <- factor(dff$group, levels = object$grouplist)
      gsg <- ggplot2::ggplot(dff, ggplot2::aes_string(
        x = paste0("PC", comps[1], "_score"),
        y = paste0("PC", comps[2], "_score"),
        shape = "time",
        color = "group",
        group = "group",
        linetype = "group"
      )) +
        ggplot2::geom_line() +
        ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high"))) +
        ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")))
      gsg <- .get_plot_handle(
        g = gsg,
        object = object,
        myTheme = myTheme,
        comp = comps,
        effect = "group",
        legend = "bottom"
      )
      
      g <- ggpubr::ggarrange(g, gsg, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(gsg), align = "hv")
    }
  } else if (validationshape == "ellipse" & !is.na(validationshape)) {
    # Score plots
    ## Time
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC %in% comps),
        model = x
      )
    }))
    df_temp <- subset(get_scores(object)$time, PC %in% comps)
    df_temp$model <- 0
    df_temp$low <- NULL
    df_temp$high <- NULL
    dff <- rbind(dff, df_temp)
    
    dff$PC <- paste0("PC", dff$PC)
    dfff <- reshape2::dcast(data = dff, time + group + model ~ PC, value.var = "score")
    dfff$time <- factor(dfff$time, levels = object$timelist)
    dfff$group <- factor(dfff$group, levels = object$grouplist)
    dfff$alpha <- ifelse(dfff$model == 0, 1, alphavalidate)
    dfff <- dfff[order(dff$time), ]
    g <- ggplot2::ggplot(dfff, ggplot2::aes_string(paste0("PC", comps[1]), paste0("PC", comps[2]), shape = "time", color = "group")) +
      ggplot2::geom_line(data = subset(dfff, model == 0), ggplot2::aes(group = paste(model, group), linetype = group)) +
      ggplot2::geom_point(alpha = dfff$alpha) +
      ggplot2::stat_ellipse(level = validationlevel)
    if (plottext) {
      g <- g + ggplot2::geom_text(
        data = subset(dfff, model == 0),
        color = "black",
        ggplot2::aes(label = time), hjust = texthjust
      )
    }
    g <- .get_plot_handle(g = g, object = object, myTheme = myTheme, comp = comps, effect = "time")
    if (object$separate_time_and_group) {
      ## Group
      dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
        data.frame(
          subset(get_scores(object$validation$temp_objects[[x]])$group, PC %in% comps),
          model = x
        )
      }))
      df_temp <- subset(get_scores(object)$group, PC %in% comps)
      df_temp$model <- 0
      df_temp$low <- NULL
      df_temp$high <- NULL
      dff <- rbind(dff, df_temp)
      
      dff$PC <- paste0("PC", dff$PC)
      dfff <- reshape2::dcast(data = dff, time + group + model ~ PC, value.var = "score")
      dfff$time <- factor(dfff$time, levels = object$timelist)
      dfff$group <- factor(dfff$group, levels = object$grouplist)
      dfff$alpha <- ifelse(dfff$model == 0, 1, alphavalidate)
      dfff <- dfff[order(dff$time), ]
      gsg <- ggplot2::ggplot(dfff, ggplot2::aes_string(paste0("PC", comps[1]), paste0("PC", comps[2]), shape = "time", color = "group")) +
        ggplot2::geom_line(data = subset(dfff, model == 0), ggplot2::aes(group = paste(model, group), linetype = group)) +
        ggplot2::geom_point(alpha = dfff$alpha) +
        ggplot2::stat_ellipse(level = validationlevel)
      if (plottext) {
        gsg <- gsg + ggplot2::geom_text(data = subset(dfff, model == 0), color = "black", ggplot2::aes(label = time), hjust = texthjust)
      }
      
      gsg <- .get_plot_handle(g = gsg, object = object, myTheme = myTheme, comp = comps, effect = "group", legend = "bottom")
      
      g <- ggpubr::ggarrange(g, gsg, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(gsg), align = "hv")
    }
  } else {
    dff <- subset(get_scores(object)$time, PC %in% comps)
    dff$PC <- paste0("PC", dff$PC)
    dff <- reshape2::dcast(data = dff, time + group ~ PC, value.var = "score")
    dff$time <- factor(dff$time, levels = object$timelist)
    dff$group <- factor(dff$group, levels = object$grouplist)
    g <- ggplot2::ggplot(dff, ggplot2::aes_string(x = paste0("PC", comps[1]), y = paste0("PC", comps[2]), shape = "time", color = "group", group = "group", linetype = "group")) +
      ggplot2::geom_line() +
      ggplot2::geom_point()
    g <- .get_plot_handle(g = g, object = object, myTheme = myTheme, comp = comps, effect = "time")
    if (object$separate_time_and_group) {
      dff <- subset(get_scores(object)$group, PC %in% comps)
      dff$PC <- paste0("PC", dff$PC)
      dff <- reshape2::dcast(data = dff, time + group ~ PC, value.var = "score")
      dff$time <- factor(dff$time, levels = object$timelist)
      dff$group <- factor(dff$group, levels = object$grouplist)
      gsg <- ggplot2::ggplot(dff, ggplot2::aes_string(x = paste0("PC", comps[1]), y = paste0("PC", comps[2]), shape = "time", color = "group", "group" = "group", linetype = "group")) +
        ggplot2::geom_line() +
        ggplot2::geom_point()
      gsg <- .get_plot_handle(
        g = gsg,
        object = object,
        myTheme = myTheme,
        comp = comps,
        effect = "group",
        legend = "bottom"
      )
      
      g <- ggpubr::ggarrange(g, gsg, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(gsg))
    }
  }
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  
  return(g)
}

#' Get plotting features
#'
#' This function returns a list with colors for plotting
#'
#' @param object An ALASCA object
#' @return A list with colors
#'
#' @export
.get_plot_handle <- function(g, object, myTheme, effect = "time", comps = c(1, 2), legend = NA) {
  g <- g + ggplot2::scale_color_manual(values = get_plot_palette(object)) +
    ggplot2::labs(
      x = .get_explained_label(object, component = comps[1], effect = effect),
      y = .get_explained_label(object, component = comps[2], effect = effect)
    ) +
    myTheme
  
  if (!is.na(legend)) g <- g + ggplot2::theme(legend.position = legend)
  
  return(g)
}



#' Plot projection of participants
#'
#' This function returns a plot of...
#'
#' @param object An ALASCA object
#' @param comp Which two components to plot (default: `c(1, 2`)
#' @param return_data Set to `TRUE` to return data instead of plot
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param myTheme A ggplot2 theme to use
#' @return A ggplot2 objects.
#'
#' @export
plot_projection <- function(object,
                            comp = c(1, 2),
                            return_data = FALSE,
                            filename = NA,
                            filetype = NA,
                            figsize = NA,
                            figunit = NA,
                            myTheme = NA) {
  if (any(is.na(myTheme))) myTheme <- object$plot.myTheme
  if (!is.na(filename)) object$filename <- filename
  
  df <- object$df
  df$ID <- df[, ID]
  loadings_Time <- subset(get_loadings(object)$time, PC %in% comp)
  loadings_Time <- reshape2::dcast(data = loadings_Time, covars ~ paste0("PC", PC), value.var = "loading")
  df_time <- merge(df, loadings_Time, by.x = "variable", by.y = "covars")
  if (object$separate_time_and_group) {
    df_time <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x) {
      data.frame(
        part = subset(df_time, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_time, paste0(ID, time) == x)$PC1 * subset(df_time, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_time, paste0(ID, time) == x)$PC2 * subset(df_time, paste0(ID, time) == x)$value),
        time = subset(df_time, paste0(ID, time) == x)$time
      )
    }))
    df_time <- df_time[!duplicated(df_time), ]
    colnames(df_time) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time")
    if (!return_data) {
      g_t <- ggplot2::ggplot(df_time, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part")) +
        ggplot2::geom_point() +
        ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.20, "cm"))) +
        myTheme
    }
    loadings_group <- subset(get_loadings(object)$group, PC %in% comp)
    loadings_group <- reshape2::dcast(data = loadings_group, covars ~ paste0("PC", PC), value.var = "loading")
    df_group <- merge(df, loadings_group, by.x = "variable", by.y = "covars")
    df_group <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x) {
      data.frame(
        part = subset(df_group, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_group, paste0(ID, time) == x)$PC1 * subset(df_group, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_group, paste0(ID, time) == x)$PC2 * subset(df_group, paste0(ID, time) == x)$value),
        time = subset(df_group, paste0(ID, time) == x)$time,
        group = subset(df_group, paste0(ID, time) == x)$group
      )
    }))
    df_group <- df_group[!duplicated(df_group), ]
    colnames(df_group) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time", "group")
    if (!return_data) {
      g_g <- ggplot2::ggplot(df_group, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part", color = "group")) +
        ggplot2::geom_point() +
        ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.20, "cm"))) +
        myTheme
      g <- ggpubr::ggarrange(g_t, g_g, common.legend = TRUE, legend = "bottom")
    }
    g <- list(df_time, df_group)
    names(g) <- c("time", "group")
    if (object$save) {
      for (i in seq_along(g)) {
        saveALASCAPlot(object = object, g = g[[i]], filetype = filetype, figsize = figsize, figunit = figunit, suffix = names(g)[i])
      }
    }
    return(g)
  } else {
    df_time <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x) {
      data.frame(
        part = subset(df_time, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_time, paste0(ID, time) == x)$PC1 * subset(df_time, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_time, paste0(ID, time) == x)$PC2 * subset(df_time, paste0(ID, time) == x)$value),
        time = subset(df_time, paste0(ID, time) == x)$time,
        group = subset(df_time, paste0(ID, time) == x)$group
      )
    }))
    df_time <- df_time[!duplicated(df_time), ]
    colnames(df_time) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time", "group")
    g <- ggplot2::ggplot(df_time, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part", color = "group")) +
      ggplot2::geom_point() +
      ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.20, "cm"))) +
      myTheme
    if (return_data) {
      return(df_time)
    } else {
      
      if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
      
      return(g)
    }
  }
}
