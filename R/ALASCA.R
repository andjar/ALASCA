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
  rows_by_variable <- lapply(self[["variablelist"]], function(x) self[["df"]][, .I[variable == x]] )
  names(rows_by_variable) <- self[["variablelist"]]
  
  if (self$use_Rfast && self$method %in% c("LMM")) {
    # start.time <- Sys.time()
    if (any(is.na(self[["df"]][, value]))) {
      log4r::error("Rfast does NOT like NA's! Check your scaling function or value column.")
      stop()
    }
    if (!self[["minimize_object"]]) {
      self[["modmat"]] <- model.matrix(self[["new_formula"]], data = self$df)
      if (self[["equal_baseline"]]) {
        self[["modmat"]] <- self[["modmat"]][, !grepl(paste0("time", self[["timelist"]][1]), colnames(self[["modmat"]]))]
      }
    }
    
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
    colnames(self$regression_coefficients) <- as.character(self[["variablelist"]])
    self$regression_coefficients[, variable := colnames(self[["modmat"]])]
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
  self[["pca"]][["time"]] <- prcomp(
    self$effect_matrix[self$effect_matrix$comp == "TIME",
                       seq_len(ncol(self$effect_matrix) - 1)],
    scale = FALSE,
    center = !self$scale_function.center
  )
  if (self$separate_time_and_group) {
    self[["pca"]][["group"]] <- prcomp(
      self$effect_matrix[self$effect_matrix$comp == "GROUP",
                         seq_len(ncol(self$effect_matrix) - 1)],
      scale = FALSE,
      center = !self$scale_function.center
    )
  }
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
  
  # Clean scores ----
  PC_time <- as.data.frame(self$pca$time$x)
  PC_time$time <- self$parts$time
  if (self$separate_time_and_group) {
    PC_time$group <- self$grouplist[1]
  }else{
    PC_time$group <- self$parts$group
  }
  self$pca$score$time <- setDT(PC_time[!duplicated(paste(PC_time$time, PC_time$group)),])
  setkey(self$pca$score$time, time, group)
  self$pca$score$explained$time <- self$pca$time$sdev^2 / sum(self$pca$time$sdev^2)
  
  # Clean loadings ----
  self$pca$loading$time <- setDT(as.data.frame(self$pca$time$rotation), keep.rownames="covars")
  self$pca$loading$explained$time <- self$pca$score$explained$time
  
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
  
  setkey(self$pca$loading$time, covars)
  
  PCloading <- paste0("PC", self$get_relevant_pcs( effect = "time"))
  for (i in PCloading) {
    # Ensure that the highest loading has positive sign
    nVar <- self$pca$loading$time[, .I[which.max(abs(get(i)))]]
    sVar <- self$pca$loading$time[nVar, sign(get(i))]
    set(self$pca$loading$time, j = (i), value = self$pca$loading$time[, get(i)] * sVar)
    set(self$pca$score$time, j = (i), value = self$pca$score$time[, get(i)] * sVar)
  }
  
  if (self$separate_time_and_group) {
    # Clean scores ----
    PC_group <- as.data.frame(self$pca$group$x)
    PC_group$time <- rownames(PC_group)
    PC_group$group <- self$parts$group
    PC_group$time <- self$parts$time
    self$pca$score$group <- setDT(PC_group[!duplicated(paste(PC_group$time, PC_group$group)),])
    setkey(self$pca$score$group, time, group)
    self$pca$score$explained$group <- self$pca$group$sdev^2 / sum(self$pca$group$sdev^2)
    
    # Clean loadings ----
    self$pca$loading$group <- setDT(as.data.frame(self$pca$group$rotation), keep.rownames="covars")
    self$pca$loading$explained$group <- self$pca$score$explained$group
    
    if(self$reduce_dimensions){
      # Loadings must be back-transformed
      self$pca$loading$group <- setDT(
        as.data.frame(
          as.matrix(self$Limm$loadings) %*% as.matrix(self$pca$loading$group[order(as.numeric(substr(covars, 3, nchar(covars)))), !"covars"])
        ),
        keep.rownames="covars")
    }
    
    setkey(self$pca$loading$group, covars)
    
    PCloading <- paste0("PC", self$get_relevant_pcs( effect = "group"))
    for (i in PCloading) {
      # Ensure that the highest loading has positive sign
      nVar <- self$pca$loading$group[, .I[which.max(abs(get(i)))]]
      sVar <- self$pca$loading$time[nVar, sign(get(i))]
      set(self$pca$loading$group, j = (i), value = self$pca$loading$group[, get(i)] * sVar)
      set(self$pca$score$group, j = (i), value = self$pca$score$group[, get(i)] * sVar)
    }
  }
  
  log4r::debug(self$log, "Completed clean_pca")
  
  #invisible(self)
}

clean_alasca <- function() {
  
  log4r::debug(self$log, "Starting clean_alasca")
  # We need to create new group names for the combined group and keep_terms
  if (self$keep_terms != "") {
    if (self$separate_time_and_group) {
      self$grouplist <- unique(self$pca$score$group$group)
    } else {
      self$grouplist <- unique(self$pca$score$time$group)
    }
  }
  
  # Clean up a copy
  # Time effect
  self$ALASCA$loading$time <- melt(self$pca$loading$time, id.vars = "covars", variable.factor = FALSE)
  colnames(self$ALASCA$loading$time) <- c("covars", "PC", "loading")
  self$ALASCA$loading$time[, PC := as.integer(substr(PC, 3, nchar(PC))), ]
  
  self$ALASCA$score$time <- self$pca$score$time
  if (self$separate_time_and_group) {
    self$ALASCA$score$time$group <- self$grouplist[1]
  }
  self$ALASCA$score$time <- melt(self$ALASCA$score$time, id.vars = c("time", "group"), variable.factor = FALSE)
  colnames(self$ALASCA$score$time) <- c("time", "group", "PC", "score")
  self$ALASCA$score$time[, time := factor(time, levels = self$timelist), ]
  self$ALASCA$score$time[, group := factor(group, levels = self$grouplist), ]
  self$ALASCA$score$time[, PC := as.integer(substr(PC, 3, nchar(PC))), ]
  
  self$ALASCA$score$explained$time <- self$pca$score$explained$time
  self$ALASCA$loading$explained$time <- self$pca$loading$explained$time
  
  if (self$separate_time_and_group) {
    # Group effect
    self$ALASCA$loading$group <- melt(self$pca$loading$group, id.vars = "covars", variable.factor = FALSE)
    colnames(self$ALASCA$loading$group) <- c("covars", "PC", "loading")
    self$ALASCA$loading$group[, PC := as.integer(substr(PC, 3, nchar(PC))), ]
    
    self$ALASCA$score$group <- self$pca$score$group
    self$ALASCA$score$group <- melt(self$ALASCA$score$group, id.vars = c("time", "group"), variable.factor = FALSE)
    colnames(self$ALASCA$score$group) <- c("time", "group", "PC", "score")
    self$ALASCA$score$group[, time := factor(time, levels = self$timelist), ]
    self$ALASCA$score$group[, group := factor(group, levels = self$grouplist), ]
    self$ALASCA$score$group[, PC := as.integer(substr(PC, 3, nchar(PC))), ]
    
    self$ALASCA$score$explained$group <- self$pca$score$explained$group
    self$ALASCA$loading$explained$group <- self$pca$loading$explained$group
  }
  
  log4r::debug(self$log, "Finished clean_alasca")
  
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
        
        temp_object <- clean_alasca(temp_object)
        
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
        temp_object$clean_alasca()
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(self$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((self$n_validation_runs - ii) * time_all, 2), " seconds"))
        temp_object
      })
    }
  }
  
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
  log4r::debug(self$log, "Starting rotation of time components")
  # We are only looking at components explaining more than a predefined value
  PCloading <- target$get_relevant_pcs(effect = "time")
  PCloading_t <- paste0("PC", PCloading)
  
  # PCA can give loadings with either sign, so we have to check whether switching signs improves the rotation
  N <- length(PCloading) # Number of components to look at
  
  # Create matrix with all possible combinations of signs
  vec <- c(-1, 1)
  lst <- lapply(numeric(N), function(x) vec)
  signMatrix <- as.matrix(expand.grid(lst))
  
  # Test all combinations and calculate residuals
  signVar <- vapply(seq_len(nrow(signMatrix) / 2), function(i) {
    c <- .procrustes(
      loadings = as.matrix(t(t(self$pca$loading$time[target$pca$loading$time, ..PCloading_t]) * signMatrix[i, ])),
      target = as.matrix(target$pca$loading$time[, ..PCloading_t])
    )
    sum((target$pca$score$time[, ..PCloading_t] - 
           as.matrix(t(t(self$pca$score$time[target$pca$score$time, ..PCloading_t]) * signMatrix[i, ])) %*% solve(c$t1))^2)
  }, FUN.VALUE = numeric(1))
  
  # Find the combination that minimizes the sum of squares
  minSignVar <- which(signVar == min(signVar))[1]
  
  # Switch signs
  for (i in PCloading){
    set(self$pca$loading$time, j = PCloading_t[i], value = self$pca$loading$time[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
    set(self$pca$score$time, j = PCloading_t[i], value = self$pca$score$time[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
  }
  
  # Final rotation
  c <- .procrustes(
    loadings = as.matrix(self$pca$loading$time[target$pca$loading$time, ..PCloading_t]),
    target = as.matrix(target$pca$loading$time[, ..PCloading_t])
  )
  
  self$pca$loading$time[target$pca$loading$time, (PCloading_t) := as.data.table(c$procrust)]
  self$pca$score$time[target$pca$score$time, (PCloading_t) := as.data.table(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]
  
  log4r::debug(self$log, "Completed rotation of time components")
  if (self$separate_time_and_group) {
    log4r::debug(self$log, "Starting rotation of group components")
    # We are only looking at components explaining more than a set limit
    PCloading <- target$get_relevant_pcs(effect = "group")
    PCloading_t <- paste0("PC", PCloading)
    
    # PCA can give loadings with either sign, so we have to check whether swithcing signs improves the rotation
    N <- length(PCloading)
    vec <- c(-1, 1)
    lst <- lapply(numeric(N), function(x) vec)
    signMatrix <- as.matrix(expand.grid(lst))
    signVar <- vapply(seq_len(nrow(signMatrix) / 2), function(i) {
      c <- .procrustes(
        loadings = as.matrix(t(t(self$pca$loading$group[target$pca$loading$group, ..PCloading_t]) * signMatrix[i, ])),
        target = as.matrix(target$pca$loading$group[, ..PCloading_t])
      )
      sum((target$pca$score$group[, ..PCloading_t] - 
             as.matrix(t(t(self$pca$score$group[target$pca$score$group, ..PCloading_t]) * signMatrix[i, ])) %*% solve(c$t1))^2)
    }, FUN.VALUE = numeric(1))
    
    minSignVar <- which(signVar == min(signVar))[1]
    for (i in PCloading){
      set(self$pca$loading$group, j = PCloading_t[i], value = self$pca$loading$group[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
      set(self$pca$score$group, j = PCloading_t[i], value = self$pca$score$group[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
    }
    
    c <- .procrustes(
      loadings = as.matrix(self$pca$loading$group[target$pca$loading$group, ..PCloading_t]),
      target = as.matrix(target$pca$loading$group[, ..PCloading_t])
    )
    self$pca$loading$group[target$pca$loading$group, (PCloading_t) := as.data.table(c$procrust)]
    self$pca$score$group[target$pca$score$group, (PCloading_t) := as.data.table(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]
    log4r::debug(self$log, "Completed rotation of group components")
  }
  
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
  if ("low" %in% colnames(self$ALASCA$loading$time)) {
    self$ALASCA$loading$time$low <- NULL
    self$ALASCA$loading$time$high <- NULL
    self$ALASCA$score$time$low <- NULL
    self$ALASCA$score$time$high <- NULL
    if (self$separate_time_and_group) {
      self$ALASCA$loading$group$low <- NULL
      self$ALASCA$loading$group$high <- NULL
      self$ALASCA$score$group$low <- NULL
      self$ALASCA$score$group$high <- NULL
    }
  }
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
  PC_time <- self$get_relevant_pcs(effect = "time")
  if (self$save_to_disk) {
    res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.loading' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
    df_time <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$loading$time[x$ALASCA$loading$time$PC %in% PC_time, ]), fill = TRUE)
  }
  
  self$validation$time$loading <- df_time[, as.list(quantile(loading, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(PC, covars)]
  colnames(self$validation$time$loading) <- c("PC", "covars", "low", "high")
  self$ALASCA$loading$time <- merge(self$ALASCA$loading$time, self$validation$time$loading, all.x = TRUE)
  
  if (self$separate_time_and_group) {
    PC_group <- self$get_relevant_pcs(effect = "group")
    if (self$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'group.loading' WHERE PC IN(", paste(PC_group, collapse = ", "), ")"))
      df_group <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_group <- rbindlist(lapply(objectlist, function(x) x$ALASCA$loading$group[x$ALASCA$loading$group$PC %in% PC_group, ]), fill = TRUE)
    }
    
    self$validation$group$loading <- df_group[, as.list(quantile(loading, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(PC, covars)]
    colnames(self$validation$group$loading) <- c("PC", "covars", "low", "high")
    self$ALASCA$loading$group <- merge(self$ALASCA$loading$group, self$validation$group$loading, all.x = TRUE)
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
  if (self$separate_time_and_group) {
    # Separate time and group effects
    
    PC_time <- self$get_relevant_pcs(effect = "time")
    if (self$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$time[x$ALASCA$score$time$PC %in% PC_time, ]), fill = TRUE)
    }
    
    self$validation$time$score <- df_time[, as.list(quantile(score, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(PC, time)]
    colnames(self$validation$time$score) <- c("PC", "time", "low", "high")
    self$ALASCA$score$time <- merge(self$ALASCA$score$time, self$validation$time$score, all.x = TRUE)
    self$ALASCA$score$time[, time := factor(time, levels = self$timelist), ]
    self$ALASCA$score$time[, group := factor(group, levels = self$grouplist), ]
    
    PC_group <- self$get_relevant_pcs(effect = "group")
    if (self$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'group.score' WHERE PC IN(", paste(PC_group, collapse = ", "), ")"))
      df_group <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_group <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$group[x$ALASCA$score$group$PC %in% PC_group, ]), fill = TRUE)
    }
    
    self$validation$group$score <- df_group[, as.list(quantile(score, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(PC, time, group)]
    colnames(self$validation$group$score) <- c("PC", "time", "group", "low", "high")
    self$ALASCA$score$group <- merge(self$ALASCA$score$group, self$validation$group$score, all.x = TRUE)
    self$ALASCA$score$group[, time := factor(time, levels = self$timelist), ]
    self$ALASCA$score$group[, group := factor(group, levels = self$grouplist), ]
  } else {
    # Pooled time and groups effects
    PC_time <- self$get_relevant_pcs(effect = "time")
    if (self$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$time[x$ALASCA$score$time$PC %in% PC_time, ]), fill = TRUE)
    }
    
    self$validation$time$score <- df_time[, as.list(quantile(score, probs = self$limitsCI, type = self$validation_quantile_method)), by = .(PC, time, group)]
    colnames(self$validation$time$score) <- c("PC", "time", "group", "low", "high")
    self$ALASCA$score$time <- merge(self$ALASCA$score$time, self$validation$time$score, all.x = TRUE)
    self$ALASCA$score$time[, time := factor(time, levels = self$timelist), ]
    self$ALASCA$score$time[, group := factor(group, levels = self$grouplist), ]
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
    PC <- self$ALASCA$loading$explained$time >= object$explanatorylimit
  } else {
    PC <- self$ALASCA$loading$explained$group >= object$explanatorylimit
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
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% unlist(selectedParts)
        )
      } else {
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% object$validation_ids[runN, ]
        )
      }
    } else if (self$method %in% c("LM")) {
      self$df$ID <- c(seq_len(nrow(self$df)))
      if (any(is.na(self$validation_ids))) {
        # For each group, divide the participants into n_validation_folds groups, and select n_validation_folds-1 of them
        selectedParts <- lapply(unique(self$stratification_vector), function(gr) {
          selectedParts_temp_all <- unique(self$df[object$stratification_vector == gr, ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% self$n_validation_folds
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
        
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% unlist(selectedParts)
        )
      } else {
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% object$validation_ids[runN, ]
        )
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
        row_nr = df_raw$rows_by_ID[[participants_in_bootstrap$old_id[participant]]]
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