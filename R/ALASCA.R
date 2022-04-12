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
                   wide = FALSE,
                   separate_time_and_group = FALSE,
                   scale_function = "sdall",
                   reduce_dimensions = FALSE,
                   equal_baseline = FALSE,
                   sum_coding = FALSE,
                   method = NA,
                   use_Rfast = TRUE,
                   p_adjust_method = NA,
                   participant_column = "ID",
                   x_column = "time",
                   ignore_missing = FALSE,
                   ignore_missing_covars = FALSE,
                   silent = FALSE,
                   scale_function.center = FALSE,
                   stratification_column = "group",
                   stratification_vector = NA,
                   minimize_object = FALSE,
                   limitsCI = c(0.025, 0.975),
                   pca_function = "prcomp",
                   plot.x_label = "Time",
                   plot.group_label = "Group",
                   plot.figsize = c(180, 120, 300),
                   plot.figunit = "mm",
                   plot.filetype = "png",
                   plot.palette = NA,
                   plot.loadinggroupcolumn = NA,
                   plot.loadinggroup_label = "Variable group",
                   plot.palette.end = 0.8,
                   explanatorylimit = 0.05,
                   do_debug = FALSE,
                   n_validation_folds = 7,
                   n_validation_runs = 1000,
                   validation_quantile_method = 2,
                   plot.my_theme = ggplot2::theme_classic(),
                   keep_columns = c(""),
                   keep_terms = c(""),
                   filename = NA,
                   filepath = NA,
                   reduce_dimensions.nComps= NULL,
                   reduce_dimensions.limit = 0.95,
                   save = FALSE,
                   save_to_disk = FALSE,
                   save_validation_ids = FALSE,
                   optimize_score = TRUE,
                   validate = FALSE,
                   validate_regression = TRUE,
                   validation = FALSE,
                   validation_method = "bootstrap",
                   validation_ids = NA,
                   validation_object = NA,
                   validation_assign_new_ids = FALSE,
                   validation_participants = NA) {
  if (!is.na(validation_object[1])) {
    # This is a validation run

    object <- validation_object
    # Overwrite some data

    ## Unscaled values
    object$df <- validation_object$df_raw
    if(object$reduce_dimensions){
      object$Limm$main$pca <- object$Limm$pca
      object$Limm$pca <- NULL
    }

    ## Avoid recursion
    object$validate <- FALSE

    ## Save space
    object$minimize_object <- TRUE

    ## Selected participants for this run
    object$validation_participants <- validation_participants

    ## Keep original object?
    object$validation_object <- NULL # validation_object
    
    object$variablelist <- unique(object$df$variable)
    object$timelist <- levels(object$df$time)
    object$grouplist <- levels(object$df$group)
    
    object$df <- object$df[object$validation_participants]
    object$modmat <- object$modmat[object$validation_participants,]
    #object$df[, ID := factor(ID)]

    #object$do_debug <- FALSE
  } else {
    object <- as.list(environment())
    object$df <- setDT(object$df)
    object$validate <- object$validate || object$validation
    object$validation <- NULL
    object$n_validation_runs <- ifelse(any(is.na(object$validation_ids)), object$n_validation_runs, nrow(object$validation_ids))
    object$init_time <- Sys.time()
    object$filepath <- ifelse(is.na(object$filepath), NA, ifelse(substr(object$filepath, nchar(object$filepath), nchar(object$filepath)) == "/", object$filepath,  paste0(object$filepath, "/")))
    object$save_to_disk <- object$validate && object$save_to_disk
    object$rawFormula <- object$formula
    object$keep_validation_objects <- TRUE
    object$ALASCA.version <- print_version(get = "version")
    object$ALASCA.version.date <- print_version(get = "date")
    object$log_file <- tempfile()
    object$log_level <- ifelse(object$do_debug, "DEBUG", "INFO")
    if (object$silent) {
      object$log <- log4r::logger(object$log_level, appenders = list(log4r::file_appender(object$log_file)))
    } else {
      object$log <- log4r::logger(object$log_level, appenders = list(log4r::console_appender(), log4r::file_appender(object$log_file)))
    }
    
    log4r::info(object$log, "Initializing ALASCA.")
    class(object) <- "ALASCA"
    print_version()
  }
  
  # Clean input ----

  if (!object$minimize_object) {
    log4r::info(object$log, "Has initialized the ALASCA model. Next step is to clean it and check input")
  }
  object <- sanitize_object(object)

  # Build the ALASCA model ----
  
  object <- buildModel(object)

  # To save space, we remove unnecessary embedded data ----
  if (object$minimize_object) {
    log4r::debug(object$log, "Starting to remove embedded data")
    object <- remove_embedded_data(object)
    log4r::debug(object$log, "Completed remove embedded data")
  }

  # Validate the model ----
  if (object$validate) {
    log4r::debug(object$log, "Starting validation")
    object <- validate(object)
    log4r::debug(object$log, "Completing validation")
  }
  object$run_time <- Sys.time() - object$init_time

  # Save the model
  if (object$save & !object$minimize_object) {
    saveALASCA(object)
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
sanitize_object <- function(object) {
  if (!object$minimize_object) {
    
    log4r::debug(object$log, "Starting sanitation")
    object <- get_info_from_formula(object)
    object <- rename_columns_to_standard(object)
    object <- wide_to_long(object)
    
    # Keep variable labels
    if (!is.na(object$plot.loadinggroupcolumn)) {
      object$variable_labels <- unique(object$df[, .SD, .SDcols = c("variable", object$plot.loadinggroupcolumn)])
      colnames(object$variable_labels) <- c("covars", "covargroup")
    }
    
    # Remove surplus data for efficiency
    object$df <- object$df[, .SD, .SDcols = c(object$all_formula_terms, "variable", "value")]
    
    check_that_columns_are_valid(object)
    
    log4r::info(object$log, "Making factors for time, group and variable")
    object$df[, time := factor(time, levels = object$timelist), ]
    object$df[, group := factor(group, levels = object$grouplist), ]
    object$df[, variable := factor(variable, levels = object$variablelist), ]
    object$stratification_vector <- object$df[, get(object$stratification_column)]
    
    # List of levels
    object$variablelist <- unique(object$df$variable)
    object$timelist <- levels(object$df$time)
    object$grouplist <- levels(object$df$group)
    
    object <- adjust_design_matrix(object)
    
    object <- get_scaling_function(object)
    object <- get_pca_function(object)

    if (object$save_to_disk) {
      object$db.driver <- RSQLite::dbDriver("SQLite")
      object$db.filename <- get_filename(object, prefix = "validation/", filetype = "db")
      object$db.con <- DBI::dbConnect(object$db.driver, dbname = object$db.filename)
    }
    
    log4r::info(object$log, "Checking for missing information")
    find_missing_predictor_variables(object)
    find_missing_response_variables(object)
    
    # Use sum coding?
    if (object$sum_coding) {
      log4r::info(object$log, "Use sum coding")
      contrasts(object$df$group) <- contr.sum(length(unique(object$df$group)))
    }
  }

  # Keep a copy of unscaled data
  object$df_raw <- object$df
  
  # Scale data
  object$df <- object$scale_function(object$df)
  
  log4r::debug(object$log, "Completed sanitation")
  return(object)
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
rename_columns_to_standard <- function(object) {
  
  if (object$valCol != "value") {
    if ("value" %in% colnames(object$df)) {
      log4r::error(object$log, "Sorry, the value column is reserved by ALASCA; please give it another name or change `valCol`")
      stop()
    }
    log4r::warn(object$log, paste0("Changing ", object$valCol, " to `value`."))
    object$df[, value := get(object$valCol)]
    object$formula <- formula(paste(
      "value ~",
      as.character(object$formula)[3]
    ))
    object <- get_info_from_formula(object)
  }
  
  if (object$x_column == "time" && ! "time" %in% object$all_formula_terms) {
    object$x_column <- object$formula_terms[1]
    log4r::info(object$log, paste("Will use", object$x_column, "for abscissa"))
  }
  
  if (object$x_column != "time") {
    if ("time" %in% colnames(object$df)) {
      object$time <- NULL
      log4r::warn(object$log, "Overwriting the `time` column")
    }
    object$df[, time := factor(get(object$x_column))]
    object$timelist <- levels(object$df$time)
    object <- replace_term_in_formula(object, old_term = object$x_column, new_term = "time")
  }
  
  
  if (object$method == "LM" && !"group" %in% colnames(object$df)) {
    object$grouplist <- object$timelist
    object$df[, group := time]
  } else if (!"group" %in% colnames(object$df)) {
    object$df[, group := factor("NA")]
  } 
  
  object$variablelist <- unique(object$df$variable)
  object$timelist <- levels(object$df$time)
  object$grouplist <- levels(object$df$group)
  
  log4r::info(object$log, paste0("Using `",object$stratification_column,"` for stratification"))

  return(object)
}

#' Check if response variables are missing
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
find_missing_response_variables <- function(object) {
  if(object$df[, uniqueN(variable), by = .(ID, time)][, uniqueN(V1)] > 1) {
    if (object$ignore_missing) {
      log4r::warn(object$log, "Response variables missing for some samples! Continue with caution!")
    } else {
      log4r::error(object$log, "Response variables missing for some samples! To ignore this, use `ignore_missing = TRUE`")
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
find_missing_predictor_variables <- function(object) {
  if (any(is.na(object$df))) {
    if (object$ignore_missing_covars) {
      log4r::warn(object$log, "Predictor variables missing for some samples! Continue with caution!")
    } else {
      log4r::error(object$log, "Predictor variables missing for some samples! To ignore this, use `ignore_missing_covars = TRUE`")
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
replace_term_in_formula <- function(object, old_term, new_term) {
  old_formula <- as.character(object$formula)[3]
  new_formula <- gsub(old_term, new_term, old_formula)
  log4r::info(object$log, paste0("New formula: ", new_formula))
  object$formula <- formula(paste(
    "value ~",
    new_formula
  ))
  object <- get_info_from_formula(object)
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
check_that_columns_are_valid <- function(object) {
  if (any(!object$all_formula_terms %in% colnames(object$df))) {
    log4r::error(object$log, paste0("Column(s) missing:\n", paste0(object$all_formula_terms[!object$all_formula_terms %in% colnames(object$df)], collapse = "\n* "), "\nYou may want to use `keep_columns`"))
    stop()
  }
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
adjust_design_matrix <- function(object) {
  if (object$method %in% c("LMM")) {
    if (object$participant_column != "ID") {
      
      # The user has specified a column
      object$df[, ID := get(object$participant_column)]
      object <- replace_term_in_formula(object, old_term = object$participant_column, new_term = "ID")
      
    } else if (any(grepl("\\|ID", object$formula_terms))) {
      
      # Use ID for participants!
      
    } else if (sum(grepl("\\|", object$formula_terms)) > 1) {
      log4r::error(object$log, "Multiple random effects, couldn't determine participant-id. Please specify `participant_column`")
      stop()
    } else {
        
        # Try to find ID column from formula
        tmp <- object$formula_terms[grepl("\\|", object$formula_terms)]
        tmp <- gsub(" ", "", tmp)
        tmp <- strsplit(tmp, "\\|")
        object$participant_column <- tmp[[1]][2]
        object$df[, ID := get(object$participant_column)]
        object <- replace_term_in_formula(object, old_term = object$participant_column, new_term = "ID")
    }
    
    if (object$use_Rfast) {
      # Using Rfast
      fixed_terms <- object$formula_terms[!grepl("\\|", object$formula_terms)]
      object$new_formula <- formula(paste("value ~ ", paste(fixed_terms, collapse = "+")))
    } else {
      # Using lme4
      rterms <- object$formula_terms[grepl("\\|", object$formula_terms)]
      rterms <- paste0("(", rterms, ")")
      object$new_formula <- formula(paste("value ~ modmat+", paste(rterms, collapse = "+")))
    }
  }else if (object$method %in% c("LM")) {
    object$new_formula <- value ~ modmat
  } else {
    log4r::error(object$log, "Sorry, an error occurred! Please check your model")
    stop()
  }
  return(object)
}

#' Get information from formula
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_info_from_formula <- function(object) {
  
  # Response variable
  object$valCol <- as.character(object$formula)[2]
  
  # Terms in the regression model
  object$formula_terms <- gsub(" ", "", colnames(attr(terms.formula(object$formula), "factors")))
  
  # Get a list of all predictors
  object$all_formula_terms <- unlist(strsplit(object$formula_terms, split = "\\:|\\+|\\||\\*"))
  object$all_formula_terms <- gsub(" ", "", object$all_formula_terms)
  object$all_formula_terms <- unique(object$all_formula_terms[!object$all_formula_terms %in% c("1", "")])
  if (object$stratification_column == "group" && !"group" %in% colnames(object$df)) {
    object$stratification_column <- object$all_formula_terms[1]
    log4r::warn(object$log, paste("The `",object$all_formula_terms[1],"` column is used for stratification"))
  }
  if (!"group" %in% object$all_formula_terms) {
    object$all_formula_terms <- c(object$all_formula_terms, "group")
    if (object$stratification_column == "group") {
      log4r::warn(object$log, "The `group` column is used for stratification")
    }
  }
  object$all_formula_terms <- c(object$all_formula_terms, object$participant_column, object$stratification_column, object$keep_columns)
  object$all_formula_terms <- gsub(" ", "", object$all_formula_terms)
  object$all_formula_terms <- unique(object$all_formula_terms[!object$all_formula_terms %in% c("1", "")])
  
  ## We need to keep original IDs to have a unique identifier later on
  if (object$validation_method == "bootstrap") {
    object$all_formula_terms <- unique(c(object$all_formula_terms, "originalIDbeforeBootstrap"))
    object$df[, originalIDbeforeBootstrap := -1]
    object$all_formula_terms <- unique(c(object$all_formula_terms, "uniqueIDforBootstrap"))
    object$df[, uniqueIDforBootstrap := -1]
  }
  
  # Check what terms that is present in formula
  object$hasGroupTerm <- ifelse(any(object$formula_terms == "group"), TRUE, FALSE)
  object$hasInteractionTerm <- ifelse(any(object$formula_terms == "group:time" | object$formula_terms == "time:group"), TRUE, FALSE)
  object$covars <- object$formula_terms[!(object$formula_terms %in% c("time", "group", "group:time", "time:group", object$keep_terms))]
  if (object$do_debug) {
    cat(".... Group term in formula? ", object$hasGroupTerm, "\n")
    cat(".... Interaction term in formula? ", object$hasInteractionTerm, "\n")
    cat(".... Identified the following covariates in addition to time and group: ", object$covars, "\n")
  }
  
  object <- LM_or_LMM(object)
  
  return(object)
}

#' Determine whether to use LM or LMM
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
LM_or_LMM <- function(object) {
  if (is.na(object$method)) {
    # Find which default method to use
    if (any(grepl("\\|", object$formula_terms))) {
      object$method <- "LMM"
      log4r::info(object$log, "Will use linear mixed models!")
      if (sum(grepl("\\|", object$formula_terms)) > 1 && object$use_Rfast) {
        log4r::error(object$log, "Cannot use Rfast with multiple random effects. Use lme4 with `use_Rfast = FALSE` instead!")
        stop()
      }
      if (!any(grepl("1\\|ID", object$formula_terms)) && object$use_Rfast) {
        log4r::error(object$log, "Rfast only supports a single random intercept. Use lme4 with `use_Rfast = FALSE` instead!")
        stop()
      }
    } else {
      object$method <- "LM"
      log4r::info(object$log, "Will use linear models!")
    }
  } else {
    # The user has specified a method to use
    if (object$method == "LMM") {
      if (!any(grepl("\\|", object$formula_terms))) {
        log4r::error(object$log, "The model must contain at least one random effect. Are you sure you wanted linear mixed models?")
        stop()
      }
    } else if (object$method == "LM") {
      if (any(grepl("\\|", object$formula_terms))) {
        log4r::error(object$log, "The model contains at least one random effect. Are you sure you wanted linear models?")
        stop()
      }
    } else {
      log4r::error(object$log, "You entered an undefined method. Use `LMM` or `LM`!")
      stop()
    }
  }
  if (object$use_Rfast) {
    log4r::info(object$log, "Will use Rfast!")
  }
  return(object)
}

#' Convert wide df to long df
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
wide_to_long <- function(object) {
  if (object$wide) {
    log4r::info(object$log, "Converting from wide to long!")
    object$df <- melt(object$df, id.vars = c(object$all_formula_terms[!object$all_formula_terms %in% c("variable", "value")]), value.name = object$valCol)
    log4r::info(object$log, paste0("Found ",length(unique(object$df$variable))," variables"))
    object$wide <- FALSE
    object$variablelist <- unique(object$df$variable)
    object$stratification_vector = object$df[, get(object$stratification_column)]
  }
  return(object)
}

#' Remove df from object
#'
#' This function removes unnecessary data
#'
#' @param object An ALASCA object
#' @return An ALASCA object
remove_embedded_data <- function(object) {
  object$partID <- object$df$ID
  object$bootPartID <- object$df$originalIDbeforeBootstrap
  object$df <- NULL
  object$df_raw <- NULL
  object$parts <- NULL
  object$validation_participants <- NULL
  object$stratification_vector <- NULL
  object$parts_with_variable <- NULL
  object$validation_object <- NULL
  object$regression_model <- NULL
  # object$regression_coefficients <- NULL
  object$effect_matrix <- NULL

  attr(object$new_formula, ".Environment") <- NULL
  attr(object$formula, ".Environment") <- NULL
  return(object)
}

#' Get a scaling function
#'
#' Return scaling function
#'
#' @param scale_function_string String to define scaing function: `sdall`, `sdref`, `sdt1`, `sdreft1`
#' @param scale_function.center Boolean. Mean centering
#' @return A scaling function
get_default_scaling_function <- function(object) {
  scale_function_string <- object$scale_function
  scale_function.center <- object$scale_function.center
  if (scale_function_string == "sdall") {
    if (scale_function.center) {
      scale_function <- function(df) {
        # Scale by the SD of all rows
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value), by = variable]
      }
    } else {
      scale_function <- function(df) {
        # Scale by the SD of all rows
        df[, value := as.double(value)][, value := value / sd(value), by = variable]
      }
    }
  } else if (scale_function_string == "sdref") {
    if (scale_function.center) {
      scale_function <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value[group == levels(group)[1]]), by = variable]
      }
    } else {
      scale_function <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1]]), by = variable]
      }
    }
  } else if (scale_function_string == "sdt1") {
    if (scale_function.center) {
      scale_function <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[time == levels(time)[1]]), by = variable]
      }
    } else {
      scale_function <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := value / sd(value[time == levels(time)[1]]), by = variable]
      }
    }
  } else if (scale_function_string == "sdreft1") {
    if (scale_function.center) {
      scale_function <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
      }
    } else {
      scale_function <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
      }
    }
  } else {
    log4r::error(object$log, "Unknown scaling method. Please use of one the following: `none`, `sdall`, `sdref`, `sdreft1`, `sdt1`")
    stop()
  }
  return(scale_function)
}

#' Get a scaling function
#'
#' Return scaling function
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_scaling_function <- function(object) {
  if (is.function(object$scale_function)) {
    # The user provided a custom function
    if (!object$minimize_object) {
      log4r::info(object$log, "Scaling data with custom function...")
    }
  } else if (object$scale_function == "none") {
    # The user do not want to scale
    if (!object$minimize_object) {
      log4r::warn(object$log, "Not scaling data...")
    }
    object$scale_function <- identity()
  } else if (is.character(object$scale_function)) {
    # Use a deafult scaling
    if (!object$minimize_object) {
      log4r::info(object$log, paste("Scaling data with",object$scale_function,"..."))
    }
    object$scale_function <- get_default_scaling_function(object)
  } else {
    log4r::error(object$log, "Unknown scaling function")
    stop()
  }
  return(object)
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
get_pca_function <- function (object) {
  if (is.function(object$pca_function)) {
    object$function.pca <- object$pca_function
  } else if (is.character(object$pca_function)) {
    if (object$pca_function == "prcomp") {
      object$function.pca <- function(df, center = TRUE) {
        prcomp(df, scale = FALSE, center = center)
      }
    } else if (object$pca_function == "irlba") {
      object$function.pca <- function(df, center = TRUE) {
        k <- irlba::prcomp_irlba(df, scale = FALSE, center = center, n = floor(0.9*min(dim(df))))
        rownames(k$rotation) <- colnames(df)
        k
      }
    } else if (object$pca_function == "princomp") {
      object$function.pca <- function(df, center = TRUE) {
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
      log4r::error(object$log, "Unknown PCA function")
      stop()
    }
  } else {
    log4r::error(object$log, "Unknown PCA function")
    stop()
  }
  return(object)
}
