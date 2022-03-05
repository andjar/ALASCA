#' Get an ALASCA object
#'
#' `ALASCA` initializes an ALASCA model and returns an ALASCA object
#'
#' This function builds your ALASCA model. It needs a data frame containing at least a column identifying participants, a column called `time` contining time information, a column `group` containing group information, a column `variable` containing variable names, and a value column. In addition you need to specify the model you want, and whether you want to separate group and time effects (defaults to `TRUE`).
#'
#' @param df Data frame to be analyzed
#' @param formula Regression model
#' @param separateTimeAndGroup Logical: should time and group effect be separated?
#' @param pAdjustMethod Method for correcting p values for multiple testing, see p.adjust.methods
#' @param validate Logical. If `TRUE`, give estimates for robustness
#' @param participantColumn String. Name of the column containing participant identification
#' @param minimizeObject Logical. If `TRUE`, remove unnecessary clutter, optimize for validation
#' @param scaleFun Either a custom function or string to define scaling function: `sdall`, `sdref`, `sdt1`, `sdreft1`
#' @param scaleFun.center Boolean. Mean centering as part of scaling
#' @param forceEqualBaseline Set to `TRUE` (default) to remove interaction between group and first time point
#' @param useSumCoding Set to `TRUE` to use sum coding instead of contrast coding for group (defaults to `FALSE`)
#' @param plot.xlabel Defaults to "Time"
#' @param plot.grouplabel Defaults to "Group"
#' @param plot.figsize A vector containing `c(width,height,dpi)` (default: `c(180, 120, 300)`)
#' @param plot.figunit Defaults to "mm"
#' @param plot.filetype Which filetype you want to save the figure to (default: `png`)
#' @param plot.palette List of colors, named by group
#' @param keepTerms Additional terms to keep in the model matrix
#' @param stratificationVector Vector of same length as `df` that specifies stratification groups during validation. Defaults to `NA`, where the group column is used.
#' @param validateRegression Whether to validate regression predictions or not (only if `validate` is `TRUE`)
#' @param doDebug Print what happens (default: `FALSE`)
#' @param save Save models and plots automatically (default: `FALSE`)
#' @param lowerLimit A data frame with lower limits for every variable
#' @param filename File name to save model and plots (when `save = TRUE`)
#' @param method Defaults to `NA` where method is either LM or LMM, depending on whether your formula contains a random effect or not
#' @param useRfast Boolean. Defaults to `TRUE`
#' @param nValFold Partitions when validating
#' @param nValRuns number of validation runs
#' @param validationMethod among  `loo` (leave-one-out, default)
#' @param validationObject Don't worry about me:)
#' @param validationParticipants Don't worry about me:)
#' @return An ALASCA object
#'
#' @examples
#' load("PE.Rdata")
#' model <- ALASCA(df = df, formula = value ~ time * group + (1 | ID))
#' @export
ALASCA <- function(df,
                   formula,
                   wide = FALSE,
                   separateTimeAndGroup = FALSE,
                   pAdjustMethod = NA,
                   participantColumn = "ID",
                   validate = FALSE,
                   scaleFun = "sdall",
                   scaleFun.center = TRUE,
                   reduceDimensions = FALSE,
                   forceEqualBaseline = FALSE,
                   useSumCoding = FALSE,
                   method = NA,
                   xColumn = "time",
                   ignoreMissing = FALSE,
                   ignoreMissingCovars = FALSE,
                   useRfast = TRUE,
                   stratificationColumn = "group",
                   stratificationVector = NA,
                   minimizeObject = FALSE,
                   limitsCI = c(0.025, 0.975),
                   plot.xlabel = "Time",
                   plot.grouplabel = "Group",
                   plot.figsize = c(180, 120, 300),
                   plot.figunit = "mm",
                   plot.filetype = "png",
                   plot.palette = NA,
                   plot.loadinggroupcolumn = NA,
                   plot.loadinggrouplabel = "Variable group",
                   plot.palette.end = 0.8,
                   explanatorylimit = 0.05,
                   doDebug = FALSE,
                   saveValidationIDs = FALSE,
                   nValFold = 7,
                   nValRuns = 1000,
                   validationQuantileMethod = 2,
                   plot.myTheme = ggplot2::theme_classic(),
                   keepTerms = c(""),
                   keepColumn = c(""),
                   save = FALSE,
                   validation = FALSE,
                   filename = NA,
                   filepath = NA,
                   reduceDimensions.nComps= NULL,
                   reduceDimensions.limit = 0.95,
                   lowerLimit = NA,
                   savetodisk = FALSE,
                   optimizeScore = TRUE,
                   validateRegression = TRUE,
                   validationMethod = "bootstrap",
                   validationIDs = NA,
                   validationObject = NA,
                   validationAssignNewID = FALSE,
                   validationParticipants = NA) {
  if (!is.na(validationObject[1])) {
    # This is a validation run

    object <- validationObject
    # Overwrite some data

    ## Unscaled values
    object$df <- validationObject$dfRaw
    if(object$reduceDimensions){
      object$Limm$main$pca <- object$Limm$pca
      object$Limm$pca <- NULL
    }

    ## Avoid recursion
    object$validate <- FALSE

    ## Save space
    object$minimizeObject <- TRUE

    ## Selected participants for this run
    object$validationParticipants <- validationParticipants

    ## Keep original object?
    object$validationObject <- NULL # validationObject
    
    object$variablelist <- unique(object$df$variable)
    object$timelist <- levels(object$df$time)
    object$grouplist <- levels(object$df$group)
    
    object$df <- object$df[object$validationParticipants]
    #object$df[, ID := factor(ID)]

    #object$doDebug <- FALSE
  } else {
    object <- list(
      df = setDT(df),
      formula = formula,
      wide = wide,
      separateTimeAndGroup = separateTimeAndGroup,
      pAdjustMethod = pAdjustMethod,
      participantColumn = participantColumn,
      validate = ifelse(validate | validation, TRUE, FALSE),
      scaleFun = scaleFun,
      scaleFun.center = scaleFun.center,
      forceEqualBaseline = forceEqualBaseline,
      useSumCoding = useSumCoding,
      method = method,
      xColumn = xColumn,
      plot.xlabel = plot.xlabel,
      plot.grouplabel = plot.grouplabel,
      plot.figsize = plot.figsize,
      plot.figunit = plot.figunit,
      plot.filetype = plot.filetype,
      plot.palette = plot.palette,
      plot.palette.end = plot.palette.end,
      reduceDimensions.nComps = reduceDimensions.nComps,
      reduceDimensions = reduceDimensions,
      reduceDimensions.limit = reduceDimensions.limit,
      plot.loadinggroupcolumn = plot.loadinggroupcolumn,
      plot.loadinggrouplabel = plot.loadinggrouplabel,
      plot.myTheme = plot.myTheme,
      explanatorylimit = explanatorylimit,
      limitsCI = limitsCI,
      ignoreMissing = ignoreMissing,
      ignoreMissingCovars = ignoreMissingCovars,
      keepColumn = keepColumn,
      minimizeObject = minimizeObject,
      doDebug = doDebug,
      nValFold = nValFold,
      nValRuns = ifelse(any(is.na(validationIDs)), nValRuns, nrow(validationIDs)),
      useRfast = useRfast,
      keepTerms = keepTerms,
      initTime = Sys.time(),
      save = save,
      saveValidationIDs = saveValidationIDs,
      lowerLimit = lowerLimit,
      filename = filename,
      filepath = ifelse(is.na(filepath), NA, ifelse(substr(filepath, nchar(filepath), nchar(filepath)) == "/", filepath,  paste0(filepath, "/"))),
      savetodisk = ifelse(validate | validation, savetodisk, FALSE),
      rawFormula = formula,
      optimizeScore = optimizeScore,
      stratificationColumn = stratificationColumn,
      keepValidationObjects = TRUE,
      validateRegression = validateRegression,
      validationMethod = validationMethod,
      validationObject = validationObject,
      validationParticipants = validationParticipants,
      validationAssignNewID = validationAssignNewID,
      validationIDs = validationIDs,
      validationQuantileMethod = validationQuantileMethod,
      ALASCA.version = printVer(get = "version"),
      ALASCA.version.date = printVer(get = "date")
    )
    object$stratificationVector = object$df[, get(object$stratificationColumn)]
    class(object) <- "ALASCA"
    printVer()
  }
  
  # Clean input ----
  if (object$doDebug) {
    cat(".. Has initialized the ALASCA model. Next step is to clean it and check input\n")
    currentTs <- Sys.time()
  }
  object <- sanitizeObject(object)
  if (object$doDebug) cat("..* sanitizeObject:", Sys.time() - currentTs, "s\n")

  # Build the ALASCA model ----
  if (object$doDebug) {
    cat(".. Has cleaned the ALASCA model. Next step is building it\n")
    currentTs <- Sys.time()
  }
  
  object <- buildModel(object)
  if (object$doDebug) cat("..* buildModel:", Sys.time() - currentTs, "s\n")

  # To save space, we remove unnecessary embedded data ----
  if (object$minimizeObject) {
    if (object$doDebug) currentTs <- Sys.time()
    object <- removeEmbedded(object)
    if (object$doDebug) cat("..* removeEmbedded:", Sys.time() - currentTs, "s\n")
  }

  # Validate the model ----
  if (object$validate) {
    if (object$doDebug) {
      cat(".. You chose to validate the model. Starting validation\n")
      currentTs <- Sys.time()
    }
    object <- validate(object)
    if (object$doDebug) cat("..* validate:", Sys.time() - currentTs, "s\n")
  }
  object$runtime <- Sys.time() - object$initTime

  # Save the model
  if (object$save & !object$minimizeObject) {
    saveALASCA(object)
  }

  if (object$savetodisk & !object$minimizeObject) {
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
printVer <- function(object = FALSE, get = NA, print = TRUE) {
  ALASCA.version <- "0.0.0.91"
  ALASCA.version.date <- "2022-03-05"
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
sanitizeObject <- function(object) {
  if (!object$minimizeObject) {
    
    object <- get_info_from_formula(object)
    object <- rename_columns_to_standard(object)
    object <- wide_to_long(object)
    
    # Remove surplus data for efficiency
    object$df <- object$df[, .SD, .SDcols = c(object$allFormulaTerms, "variable", "value")]
    
    check_that_columns_are_valid(object)
    
    if (object$doDebug) {
      cat(".... Making factors for time, group and variable\n")
    }
    object$df[, time := factor(time, levels = object$timelist), ]
    object$df[, group := factor(group, levels = object$grouplist), ]
    object$df[, variable := factor(variable, levels = object$variablelist), ]
    
    # List of levels
    object$variablelist <- unique(object$df$variable)
    object$timelist <- levels(object$df$time)
    object$grouplist <- levels(object$df$group)
    
    object <- adjust_design_matrix(object)

    if (object$savetodisk) {
      object$db.driver <- RSQLite::dbDriver("SQLite")
      object$db.filename <- getFilename(object, prefix = "validation/", filetype = "db")
      object$db.con <- DBI::dbConnect(object$db.driver, dbname = object$db.filename)
    }
  }

  if (object$doDebug) {
    cat(".... Checking for missing information\n")
  }
  
  find_missing_predictor_variables(object)

  # Use sum coding?
  if (object$useSumCoding) {
    if (object$doDebug) {
      cat(".... Use sum coding\n")
    }
    contrasts(object$df$group) <- contr.sum(length(unique(object$df$group)))
  }

  # Keep a copy of unscaled data
  object$dfRaw <- object$df
  
  
  object <- get_scaling_function(object)

  return(object)
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
rename_columns_to_standard <- function(object) {
  
  if (object$valCol != "value") {
    if ("value" %in% colnames(object$df)) stop("Sorry, the value column is reserved by ALASCA; please give it another name or change `valCol`")
    warning("Changing", object$valCol, "to `value`.")
    object$df[, value := get(object$valCol)]
    object$formula <- formula(paste(
      "value ~",
      as.character(object$formula)[3]
    ))
    object <- get_info_from_formula(object)
  }
  
  if (object$xColumn != "time") {
    if ("time" %in% colnames(object$df)) stop("Sorry, the time column is reserved by ALASCA; please give it another name or change `xColumn`")
    object$df[, time := get(object$xColumn)]
    object <- replace_term_in_formula(object, oldTerm = object$xColumn, newTerm = "time")
  }
  
  if (!"group" %in% colnames(object$df)) {
    object$df[, group := factor("NA")]
  }
  
  object$variablelist <- unique(object$df$variable)
  object$timelist <- levels(object$df$time)
  object$grouplist <- levels(object$df$group)
  
  cat("Using",object$stratificationColumn,"for stratification.\n")
  
  
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
    if (object$ignoreMissing) {
      warning("Response variables missing for some samples! Continue with caution!")
    } else {
      stop("Response variables missing for some samples! To ignore this, use `ignoreMissing = TRUE`")
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
    if (object$ignoreMissingCovars) {
      warning("Predictor variables missing for some samples! Continue with caution!")
    } else {
      stop("Predictor variables missing for some samples! To ignore this, use `ignoreMissingCovars = TRUE`")
    }
  }
}

#' Check if response variables are missing
#'
#' ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
replace_term_in_formula <- function(object, oldTerm, newTerm) {
  old_formula <- as.character(object$formula)[3]
  new_formula <- gsub(oldTerm, newTerm, old_formula)
  cat("New formula:",new_formula, "\n")
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
  if (any(!object$allFormulaTerms %in% colnames(object$df))) {
    stop("Column(s) missing:\n", paste0(object$allFormulaTerms[!object$allFormulaTerms %in% colnames(object$df)], collapse = "\n* "), "\nYou may want to use `keepColumn`")
  }
}

#' Check if columns are missing
#'
#' ...
#'
#' @param object An ALASCA object
adjust_design_matrix <- function(object) {
  if (object$method %in% c("LMM")) {
    if (object$participantColumn != "ID") {
      
      # The user has specified a column
      object$df[, ID := get(object$participantColumn)]
      object <- replace_term_in_formula(object, oldTerm = object$participantColumn, newTerm = "ID")
      
    } else if (any(grepl("\\|ID", object$formulaTerms))) {
      
      # Use ID for participants!
      
    } else if (sum(grepl("\\|", object$formulaTerms)) > 1) {
        stop("Multiple random effects, couldn't determine participant-id. Please specify `participantColumn`")
    } else {
        
        # Try to find ID column from formula
        tmp <- object$formulaTerms[grepl("\\|", object$formulaTerms)]
        tmp <- gsub(" ", "", tmp)
        tmp <- strsplit(tmp, "\\|")
        object$participantColumn <- tmp[[1]][2]
        object$df[, ID := get(object$participantColumn)]
        object <- replace_term_in_formula(object, oldTerm = object$participantColumn, newTerm = "ID")
    }
    
    if (object$useRfast) {
      # Using Rfast
      fixed_terms <- object$formulaTerms[!grepl("\\|", object$formulaTerms)]
      object$newformula <- formula(paste("value ~ ", paste(fixed_terms, collapse = "+")))
    } else {
      # Using lme4
      rterms <- object$formulaTerms[grepl("\\|", object$formulaTerms)]
      rterms <- paste0("(", rterms, ")")
      object$newformula <- formula(paste("value ~ modmat+", paste(rterms, collapse = "+")))
    }
  }else if (object$method %in% c("LM")) {
    object$newformula <- value ~ modmat
  } else {
    stop("Sorry, an error occurred! Please check your model")
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
  object$formulaTerms <- gsub(" ", "", colnames(attr(terms.formula(object$formula), "factors")))
  
  # Get a list of all predictors
  object$allFormulaTerms <- unlist(strsplit(object$formulaTerms, split = "\\:|\\+|\\||\\*"))
  object$allFormulaTerms <- c(object$allFormulaTerms, object$participantColumn, object$stratificationColumn, object$keepColumn, "group")
  object$allFormulaTerms <- gsub(" ", "", object$allFormulaTerms)
  object$allFormulaTerms <- unique(object$allFormulaTerms[!object$allFormulaTerms %in% c("1", "")])
  
  ## We need to keep original IDs to have a unique identifier later on
  if (object$validationMethod == "bootstrap") {
    object$allFormulaTerms <- unique(c(object$allFormulaTerms, "originalIDbeforeBootstrap"))
    object$df[, originalIDbeforeBootstrap := -1]
    object$allFormulaTerms <- unique(c(object$allFormulaTerms, "uniqueIDforBootstrap"))
    object$df[, uniqueIDforBootstrap := -1]
  }
  
  # If variable groups are defined, keep them for later
  if (!is.na(object$plot.loadinggroupcolumn)) {
    object$allFormulaTerms <- unique(c(object$allFormulaTerms, object$plot.loadinggroupcolumn))
  }
  
  # Check what terms that is present in formula
  object$hasGroupTerm <- ifelse(any(object$formulaTerms == "group"), TRUE, FALSE)
  object$hasInteractionTerm <- ifelse(any(object$formulaTerms == "group:time" | object$formulaTerms == "time:group"), TRUE, FALSE)
  object$covars <- object$formulaTerms[!(object$formulaTerms %in% c("time", "group", "group:time", "time:group", object$keepTerms))]
  if (object$doDebug) {
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
    if (any(grepl("\\|", object$formulaTerms))) {
      object$method <- "LMM"
      cat("Will use linear mixed models!\n")
      if (sum(grepl("\\|", object$formulaTerms)) > 1 && object$useRfast) {
        stop("Cannot use Rfast with multiple random effects. Use lme4 with `useRfast = FALSE` instead!\n")
      }
      if (!any(grepl("1\\|ID", object$formulaTerms)) && object$useRfast) {
        stop("Rfast only supports a single random intercept. Use lme4 with `useRfast = FALSE` instead!\n")
      }
    } else {
      object$method <- "LM"
      cat("Will use linear models!\n")
    }
  } else {
    # The user has specified a method to use
    if (object$method == "LMM") {
      if (!any(grepl("\\|", object$formulaTerms))) {
        stop("The model must contain at least one random effect. Are you sure you wanted linear mixed models?")
      }
    } else if (object$method == "LM") {
      if (any(grepl("\\|", object$formulaTerms))) {
        stop("The model contains at least one random effect. Are you sure you wanted linear models?")
      }
    } else {
      stop("You entered an undefined method. Use `LMM` or `LM`!")
    }
  }
  if (object$useRfast) {
    cat("Will use Rfast!\n")
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
    cat("Converting from wide to long!\n")
    object$df <- melt(object$df, id.vars = c(object$allFormulaTerms[!object$allFormulaTerms %in% c("variable", "value")]), value.name = object$valCol)
    cat("Found",length(unique(object$df$variable)),"variables\n")
    object$wide <- FALSE
    object$variablelist <- unique(object$df$variable)
    object$stratificationVector = object$df[, get(object$stratificationColumn)]
  }
  return(object)
}

#' Remove df from object
#'
#' This function removes unnecessary data
#'
#' @param object An ALASCA object
#' @return An ALASCA object
removeEmbedded <- function(object) {
  object$partID <- object$df$ID
  object$bootPartID <- object$df$originalIDbeforeBootstrap
  object$df <- NULL
  object$dfRaw <- NULL
  object$parts <- NULL
  object$validationParticipants <- NULL
  object$stratificationVector <- NULL
  object$partsWithVariable <- NULL
  object$validationObject <- NULL
  object$regr.model <- NULL
  # object$RegressionCoefficients <- NULL
  object$effect.matrix <- NULL

  attr(object$newformula, ".Environment") <- NULL
  attr(object$formula, ".Environment") <- NULL
  return(object)
}

#' Get a scaling function
#'
#' Return scaling function
#'
#' @param scaleFun_string String to define scaing function: `sdall`, `sdref`, `sdt1`, `sdreft1`
#' @param scaleFun.center Boolean. Mean centering
#' @return A scaling function
get_default_scaling_function <- function(scaleFun_string = "sdall", scaleFun.center = TRUE) {
  if (scaleFun_string == "sdall") {
    if (scaleFun.center) {
      scaleFun <- function(df) {
        # Scale by the SD of all rows
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value), by = variable]
      }
    } else {
      scaleFun <- function(df) {
        # Scale by the SD of all rows
        df[, value := as.double(value)][, value := value / sd(value), by = variable]
      }
    }
  } else if (scaleFun_string == "sdref") {
    if (scaleFun.center) {
      scaleFun <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := (value-mean(value)) / sd(value[group == levels(group)[1]]), by = variable]
      }
    } else {
      scaleFun <- function(df) {
        # Scale by the SD of all rows in the refence group
        df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1]]), by = variable]
      }
    }
  } else if (scaleFun_string == "sdt1") {
    if (scaleFun.center) {
      scaleFun <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[time == levels(time)[1]]), by = variable]
      }
    } else {
      scaleFun <- function(df) {
        # Scale by the SD of all baseline rows
        df[, value := as.double(value)][, value := value / sd(value[time == levels(time)[1]]), by = variable]
      }
    }
  } else if (scaleFun_string == "sdreft1") {
    if (scaleFun.center) {
      scaleFun <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := (value - mean(value)) / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
      }
    } else {
      scaleFun <- function(df) {
        # Scale by the SD of all baseline rows in the reference group
        df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
      }
    }
  } else {
    stop("Unknown scaling method. Please use of one the following: `none`, `sdall`, `sdref`, `sdreft1`, `sdt1`")
  }
}

#' Get a scaling function
#'
#' Return scaling function
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_scaling_function <- function(object) {
  # The user provided a custom function
  if (is.function(object$scaleFun)) {
    if (!object$minimizeObject) {
      cat("Scaling data with custom function...\n")
    }
    object$df <- object$scaleFun(object$df)
    
    # The user do not want to scale
  } else if (object$scaleFun == "none") {
    if (!object$minimizeObject) {
      warning("Not scaling data...\n")
    }
    
    # Use a deafult scaling
  } else if (is.character(object$scaleFun)) {
    object$scaleFun <- get_scaling_function(scaleFun_string = object$scaleFun, scaleFun.center = object$scaleFun.center)
    if (!object$minimizeObject) {
      cat("Scaling data...\n")
    }
    object$df <- object$scaleFun(object$df)
    
  } else {
    stop("Unknown scaling function")
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
flipIt <- function(object, component = NA, effect = "both") {
  if (any(is.na(component))) {
    component <- unique(object$ALASCA$score$time$PC)
  }

  if (effect %in% c("both", "time")) {
    object$ALASCA$score$time[PC %in% component, score := -score]
    object$ALASCA$loading$time[PC %in% component, loading := -loading]
  }
  if (object$separateTimeAndGroup & effect %in% c("both", "group")) {
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
    if (object$separateTimeAndGroup & effect %in% c("both", "group")) {
      object$ALASCA$score$group[PC %in% component, low := -low]
      object$ALASCA$score$group[PC %in% component, high := -high]
      object$ALASCA$loading$group[PC %in% component, low := -low]
      object$ALASCA$loading$group[PC %in% component, high := -high]
    }
  }

  if (object$validate && !object$savetodisk) {
    for (i in seq_along(object$validation$temp_objects)) {
      object$validation$temp_objects[[i]] <- flipIt(object$validation$temp_objects[[i]], component = component, effect = effect)
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
  cat("Model initialized ", as.character(object$initTime), " using ", object$method, " on ", length(unique(object$RegressionCoefficients$covar)), " variables. ", sep = "", file = file, append = TRUE)
  if (object$validate) {
    cat("The model been validated with ", object$validationMethod, ".\n", file = file, append = TRUE)
  } else {
    cat("The model has *not* been validated yet.\n", file = file, append = TRUE)
  }
  cat("\nRegression model: ", as.character(object$rawFormula)[2], as.character(object$rawFormula)[1], as.character(object$rawFormula)[3], "\n", file = file, append = TRUE)
  cat("Separating time and group effects: ", object$separateTimeAndGroup, "\n", file = file, append = TRUE)
  cat("Force equal baseline: ", object$forceEqualBaseline, "\n", file = file, append = TRUE)
  cat("Using Rfast: ", object$useRfast, "\n", file = file, append = TRUE)
  if (!object$useRfast) {
    cat("Adjustment of p-values: ", object$pAdjustMethod, "\n", file = file, append = TRUE)
  }
  if (!is.null(object$reduceDimensions.nComps)) {
    cat("Kept",object$reduceDimensions.nComps,"components from initial PCA, explaining",100*cumsum(object$limm.explanatory_power)[object$reduceDimensions.nComps],"% of variation\n", file = file, append = TRUE)
  }
  cat("\n\nPCs explaining at least 5% of variation:\n   Time: ",
    paste(getRelevantPCs(object = object, effect = "time"), collapse = ", "), " (",
    paste(round(100 * object$pca$score$explained$time[getRelevantPCs(object = object, effect = "time")], 2), collapse = "%, "), "%)",
    ifelse(object$separateTimeAndGroup, paste0(
      "\n   Group: ",
      paste(getRelevantPCs(object = object, effect = "group"), collapse = ", "), " (",
      paste(round(100 * object$pca$score$explained$group[getRelevantPCs(object = object, effect = "group")], 2), collapse = "%, "), "%)"
    ), "\n"),
    sep = "", file = file, append = TRUE
  )
  if (is.function(object$scaleFun)) {
    cat("\nScaling function:\n", file = file, append = TRUE)
    cat(deparse(object$scaleFun), file = file, append = TRUE)
  } else {
    cat("\nNo scaling performed.\n", file = file, append = TRUE)
  }
  if (sessioninfo & file != "") {
    write(capture.output(devtools::session_info()), file = file, append = TRUE)
  }
}
