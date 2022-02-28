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
#' @param scaleFun If `TRUE` (default), each variable is scaled to unit SD and zero mean. If `FALSE`, no scaling is performed. You can also provide a custom scaling function that has the data frame `df` as input and output
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
#' @param method Defaults to `NA` where method is either LM or LMM, depending on whether your formula contains a random effect or not. Set to KM or KMM for survival analysis
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
                   separateTimeAndGroup = FALSE,
                   pAdjustMethod = NA,
                   participantColumn = "ID",
                   validate = FALSE,
                   scaleFun = "sdall",
                   forceEqualBaseline = FALSE,
                   useSumCoding = FALSE,
                   method = NA,
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
                   keepColumn = "",
                   save = FALSE,
                   validation = FALSE,
                   filename = NA,
                   filepath = NA,
                   limm.nComps= NULL,
                   limm.limit = 0.95,
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
    if(object$method %in% c("Limm", "Lim")){
      object$Limm$main$pca <- object$Limm$pca
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

    #object$doDebug <- FALSE
  } else {
    object <- list(
      df = setDT(df),
      formula = formula,
      separateTimeAndGroup = separateTimeAndGroup,
      pAdjustMethod = pAdjustMethod,
      participantColumn = participantColumn,
      validate = ifelse(validate | validation, TRUE, FALSE),
      scaleFun = scaleFun,
      forceEqualBaseline = forceEqualBaseline,
      useSumCoding = useSumCoding,
      method = method,
      plot.xlabel = plot.xlabel,
      plot.grouplabel = plot.grouplabel,
      plot.figsize = plot.figsize,
      plot.figunit = plot.figunit,
      plot.filetype = plot.filetype,
      plot.palette = plot.palette,
      plot.palette.end = plot.palette.end,
      limm.nComps = limm.nComps,
      limm.limit = limm.limit,
      plot.loadinggroupcolumn = plot.loadinggroupcolumn,
      plot.loadinggrouplabel = plot.loadinggrouplabel,
      plot.myTheme = plot.myTheme,
      explanatorylimit = explanatorylimit,
      limitsCI = limitsCI,
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
      variablelist = unique(df$variable),
      timelist = levels(df$time),
      grouplist = levels(df$group),
      ALASCA.version = printVer(get = "version"),
      ALASCA.version.date = printVer(get = "date")
    )
    object$stratificationVector = object$df[, get(object$stratificationColumn)]
    printVer()
  }
  class(object) <- "ALASCA"
  if (object$doDebug) {
    cat(".. Has initialized the ALASCA model. Next step is to clean it and check input\n")
    currentTs <- Sys.time()
  }

  # Clean input
  object <- sanitizeObject(object)

  if (object$doDebug) {
    cat("..* sanitizeObject:", Sys.time() - currentTs, "s\n")
    cat(".. Has cleaned the ALASCA model. Next step is building it\n")
    currentTs <- Sys.time()
  }
  # Build the ALASCA model
  object <- buildModel(object)
  if (object$doDebug) cat("..* buildModel:", Sys.time() - currentTs, "s\n")

  if (object$minimizeObject) {
    # To save space, we remove unnecessary embedded data
    if (object$doDebug) currentTs <- Sys.time()
    object <- removeEmbedded(object)
    if (object$doDebug) cat("..* removeEmbedded:", Sys.time() - currentTs, "s\n")
  }

  # Validate the model
  if (object$validate) {
    if (object$doDebug) {
      cat(".. You chose to validate the model. Starting validation\n")
      if (object$doDebug) currentTs <- Sys.time()
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
  ALASCA.version <- "0.0.0.111"
  ALASCA.version.date <- "2022-02-22"
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
    object$valCol <- as.character(object$formula)[2]

    # Check formula from user
    object$formulaTerms <- colnames(attr(terms.formula(object$formula), "factors"))
    object$allFormulaTerms <- unlist(strsplit(c(object$formulaTerms, object$participantColumn, object$stratificationColumn, object$keepColumn, "group"), split = "\\:|\\+|\\||\\*"))
    object$allFormulaTerms <- gsub(" ", "", object$allFormulaTerms)
    object$allFormulaTerms <- unique(object$allFormulaTerms[object$allFormulaTerms != "1"])
    
    ## We need to keep original IDs to have a unique identifier later on
    if (object$validationMethod == "bootstrap") object$allFormulaTerms <- unique(c(object$allFormulaTerms,"originalIDbeforeBootstrap"))
    if (!"originalIDbeforeBootstrap" %in% colnames(object$df)) object$df[, originalIDbeforeBootstrap := -1]
    if (object$validationMethod == "bootstrap") object$allFormulaTerms <- unique(c(object$allFormulaTerms,"uniqueIDforBootstrap"))
    if (!"uniqueIDforBootstrap" %in% colnames(object$df)) object$df[, uniqueIDforBootstrap := -1]
    
    # Remove surplus data for efficiency
    if (is.na(object$plot.loadinggroupcolumn)) {
      object$df <- object$df[, .SD, .SDcols = c(object$allFormulaTerms, "variable", "value")]
    } else {
      object$df <- object$df[, .SD, .SDcols = c(object$allFormulaTerms, "variable", "value", object$plot.loadinggroupcolumn)]
    }
    
    
    if (!is.na(object$method)) {
      # The user has specified a method to use
      if (object$method == "LMM") {
        if (!any(grepl("\\|", object$formulaTerms))) {
          stop("The model must contain at least one random effect. Sure you wanted linear mixed models?")
        }
      } else if (object$method == "LM") {
        if (any(grepl("\\|", object$formulaTerms))) {
          stop("The model contains at least one random effect. Sure you not wanted linear mixed models instead?")
        }
      } else if (object$method %in% c("KM", "KMM", "Limm","Lim")) {

      } else {
        stop("You entered an undefined method. Use `LMM` or `LM`!")
      }
      if (object$useRfast) {
        cat("Will use Rfast!\n")

        # Validation of regression only works for LMs at the moment
        object$validateRegression <- FALSE
      }
    } else {
      # Find which default method to use
      if (any(grepl("\\|", object$formulaTerms))) {
        object$method <- "LMM"
        cat("Will use linear mixed models!\n")
      } else {
        object$method <- "LM"
        cat("Will use linear models!\n")
      }
    }

    # Check that the input is as expected
    if (!("time" %in% colnames(object$df))) {
      stop("The dataframe must contain a column names 'time'")
    }
    if (!("group" %in% colnames(object$df))) {
      stop("The dataframe must contain a column names 'group'")
    }
    if (!("variable" %in% colnames(object$df))) {
      stop("The dataframe must contain a column names 'variable'")
    }

    if (all(is.na(object$stratificationVector))) {
      cat("Using group for stratification.\n")
      object$stratificationVector <- object$df$group
    }
    if (object$method %in% c("KM", "KMM")) {
      if (!is.data.frame(object$lowerLimit) | any(!(unique(object$df$variable) %in% unique(object$lowerLimit$variable)))) {
        stop("Some lower limits are not specified!")
      }
      for (i in unique(object$df$variable)) {
        maxValue <- max(object$df$value[object$df$variable == i])
        object$df$value[object$df$variable == i] <- maxValue - object$df$value[object$df$variable == i]
        object$lowerLimit$value[object$lowerLimit$variable == i] <- maxValue - object$lowerLimit$value[object$lowerLimit$variable == i]
        object$df$belowLowerLimit[object$df$variable == i] <- object$df$value[object$df$variable == i] > object$lowerLimit$value[object$lowerLimit$variable == i]
        object$df$value[object$df$variable == i & object$df$belowLowerLimit[object$df$variable == i]] <- object$lowerLimit$value[object$lowerLimit$variable == i]
      }
    }

    # Change value column if necessary
    if (as.character(object$formula)[2] != "value") {
      cat("Changing", as.character(object$formula)[2], "to `value`.\n")
      object$df$value <- object$df[, get(object$valCol)]
      object$formula <- formula(paste(
        "value ~",
        as.character(object$formula)[3]
      ))
    }

    if (object$method %in% c("LMM", "Limm")) {
      if (object$participantColumn != "ID") {
        object$df$ID <- object$df[, get(object$participantColumn)]
        tmp <- object$formulaTerms[!grepl("\\|", object$formulaTerms)]
        object$formula <- formula(paste(
          "value ~",
          paste(tmp, collapse = "+")
        ))
      } else if (any(grepl("\\|", object$formulaTerms))) {
        if (sum(grepl("\\|", object$formulaTerms)) > 1) {
          stop("Multiple random effects, couldn't determine participant-id. Please specify `participantColumn`")
        } else {
          tmp <- object$formulaTerms[grepl("\\|", object$formulaTerms)]
          tmp <- gsub(" ", "", tmp)
          tmp <- strsplit(tmp, "\\|")
          object$participantColumn <- tmp[[1]][2]
          object$df$ID <- object$df[, get(object$participantColumn)]
          tmp <- object$formulaTerms[!grepl("\\|", object$formulaTerms)]
          object$formula <- formula(paste(
            "value ~",
            paste(tmp, collapse = "+")
          ))
        }
      } else if (object$participantColumn == "ID") {
        if (!("ID" %in% colnames(object$df))) {
          stop("Please specify participant-id in `participantColumn`")
        }
      }

      if (object$method %in% c("LMM", "Limm")) {
        if (object$useRfast) {
          # Using Rfast
          rterms <- object$formulaTerms[!grepl("\\|", object$formulaTerms)]
          object$newformula <- formula(paste("value ~ ", paste(rterms, collapse = "+")))
        } else {
          # Using lme4
          rterms <- object$formulaTerms[grepl("\\|", object$formulaTerms)]
          rterms <- paste0("(", rterms, ")")
          object$newformula <- formula(paste("value ~ modmat+", paste(rterms, collapse = "+")))
        }
      }
    } else if (object$method %in% c("LM", "Lim")) {
      object$newformula <- value ~ modmat
    }

    if (object$savetodisk) {
      object$db.driver <- RSQLite::dbDriver("SQLite")
      object$db.filename <- getFilename(object, prefix = "validation/", filetype = "db")
      object$db.con <- DBI::dbConnect(object$db.driver, dbname = object$db.filename)
      # object$numvariablelist <- 1:length(object$variablelist)
      # names(object$numvariablelist) <- object$variablelist
      # DBI::dbWriteTable(object$db.con, "covars", data.frame(id = object$numvariablelist, covar = object$variablelist))

      # object$timelist <- unique(object$df$time)
      # object$numtimelist <- 1:length(object$timelist)
      # names(object$numtimelist) <- object$timelist
      # DBI::dbWriteTable(object$db.con, "time", data.frame(id = object$numtimelist, time = object$timelist))

      # object$grouplist <- unique(object$df$group)
      # object$numgrouplist <- 1:length(object$grouplist)
      # names(object$numgrouplist) <- object$grouplist
      # DBI::dbWriteTable(object$db.con, "group", data.frame(id = object$numgrouplist, group = object$grouplist))
    }
  }

  if (object$doDebug) {
    cat(".... Making factors for time, group and variable\n")
  }
  object$df[, time := factor(time, levels = object$timelist), ]
  object$df[, group := factor(group, levels = object$grouplist), ]
  object$df[, variable := factor(variable, levels = object$variablelist), ]

  if (object$minimizeObject) {
    # This is usually a validation object
    object$df <- object$df[object$validationParticipants]
    object$df$ID <- factor(object$df$ID)
  }

  if (object$doDebug) {
    cat(".... Checking for missing information\n")
  }
  if (any(is.na(object$df$time) | is.na(object$df$group))) {
    warning("\n\n!!! -> Oh dear, at least on of your rows is missing either time or group. I have removed it/them for now, but you should check if this is a problem...\n\n")
    object$df <- object$df[!is.na(time) & !is.na(group)]
  }

  # Use sum coding?
  if (object$useSumCoding) {
    if (object$doDebug) {
      cat(".... Use sum coding\n")
    }
    contrasts(object$df$group) <- contr.sum(length(unique(object$df$group)))
  }

  # Keep a copy of unscaled data
  object$dfRaw <- object$df

  # The user provided a custom function
  if (is.function(object$scaleFun)) {
    if (!object$minimizeObject) {
      cat("Scaling data with custom function...\n")
    }
    object$df <- object$scaleFun(object$df)
    object$df$value <- object$df[, get(object$valCol)]

    # The user do not want to scale
  } else if (object$scaleFun == "none") {
    if (!object$minimizeObject) {
      warning("Not scaling data...\n")
    }

    # Use a deafult scaling
  } else if (is.character(object$scaleFun)) {
    object$scaleFun <- getScaleFun(object$scaleFun)
    if (!object$minimizeObject) {
      cat("Scaling data...\n")
    }
    object$df <- object$scaleFun(object$df)
    object$df$value <- object$df[, get(object$valCol)]
  } else {
    stop("Unknown scaling function")
  }
  if (!object$minimizeObject) {
    # Check what terms that is present in formula
    object$hasGroupTerm <- ifelse(any(object$formulaTerms == "group"), TRUE, FALSE)
    object$hasInteractionTerm <- ifelse(any(object$formulaTerms == "group:time" | object$formulaTerms == "time:group"), TRUE, FALSE)
    object$covars <- object$formulaTerms[!(object$formulaTerms %in% c("time", "group", "group:time", "time:group", object$keepTerms))]
    if (object$doDebug) {
      cat(".... Group term in formula? ", object$hasGroupTerm, "\n")
      cat(".... Interaction term in formula? ", object$hasInteractionTerm, "\n")
      cat(".... Identified the following covariates in addition to time and troup: ", object$covars, "\n")
    }
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
#' @return An ALASCA object
getScaleFun <- function(scaleFun_string) {
  if (scaleFun_string == "sdall") {
    scaleFun <- function(df) {
      # Scale by the SD of all rows
      df[, value := as.double(value)][, value := value / sd(value), by = variable]
    }
  } else if (scaleFun_string == "sdref") {
    scaleFun <- function(df) {
      # Scale by the SD of all rows in the refence group
      df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1]]), by = variable]
    }
  } else if (scaleFun_string == "sdt1") {
    scaleFun <- function(df) {
      # Scale by the SD of all baseline rows
      df[, value := as.double(value)][, value := value / sd(value[time == levels(time)[1]]), by = variable]
    }
  } else if (scaleFun_string == "sdreft1") {
    scaleFun <- function(df) {
      # Scale by the SD of all baseline rows in the reference group
      df[, value := as.double(value)][, value := value / sd(value[group == levels(group)[1] & time == levels(time)[1]]), by = variable]
    }
  } else {
    stop("Unknown scaling method. Please use of on the following: `none`, `sdall`, `sdref`, `sdreft1`, `sdt1`")
  }
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
  if (!is.null(object$limm.nComps)) {
    cat("Kept",object$limm.nComps,"components from initial PCA, explaining",100*cumsum(object$limm.explanatory_power)[object$limm.nComps],"% of variation\n", file = file, append = TRUE)
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
