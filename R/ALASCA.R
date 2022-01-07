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
#' @param plot.figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
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
#' model <- ALASCA(df = df, formula = value~time*group + (1|ID))
#' 
#' @export
ALASCA <- function(df,
                   formula,
                   separateTimeAndGroup = TRUE,
                   pAdjustMethod = NA,
                   participantColumn = "ID",
                   validate = FALSE,
                   scaleFun = TRUE,
                   forceEqualBaseline = FALSE,
                   useSumCoding = FALSE,
                   method = NA,
                   useRfast = TRUE,
                   stratificationVector = NA,
                   minimizeObject = FALSE,
                   plot.xlabel = "Time",
                   plot.grouplabel = "Group",
                   plot.figsize = c(12, 8, 300),
                   plot.figunit = "mm",
                   plot.filetype = "png",
                   plot.palette = NA,
                   plot.loadinggroupcolumn = NA,
                   doDebug = FALSE,
                   nValFold = 7,
                   nValRuns = 50,
                   keepTerms = c(""),
                   save = FALSE,
                   filename = NA,
                   filepath = NA,
                   lowerLimit = NA,
                   optimizeScore = TRUE,
                   validateRegression = FALSE,
                   validationMethod = "loo",
                   validationObject = NA,
                   validationParticipants = NA){

  if(!is.na(validationObject[1])){
    # This is a validation run
    
    object <- validationObject
    # Overwrite some data
    
    ## Unscaled values
    object$df <- validationObject$dfRaw
    
    ## Avoid recursion
    object$validate <- FALSE
    
    ## Save space
    object$minimizeObject <- TRUE
    
    ## Selected participants for this run
    object$validationParticipants <- validationParticipants
    
    ## Keep original object?
    object$validationObject <- NULL #validationObject
    
    object$doDebug <- FALSE
  }else{
    object <- list(df = data.table::as.data.table(df),
                   formula = formula,
                   separateTimeAndGroup = separateTimeAndGroup,
                   pAdjustMethod = pAdjustMethod,
                   participantColumn = participantColumn,
                   validate = validate,
                   scaleFun = scaleFun,
                   forceEqualBaseline = forceEqualBaseline,
                   useSumCoding = FALSE,
                   method = method,
                   plot.xlabel = plot.xlabel,
                   plot.grouplabel = plot.grouplabel,
                   plot.figsize = plot.figsize,
                   plot.figunit = plot.figunit,
                   plot.filetype = plot.filetype,
                   plot.palette = plot.palette,
                   plot.loadinggroupcolumn = plot.loadinggroupcolumn,
                   minimizeObject = minimizeObject,
                   doDebug = doDebug,
                   nValFold = nValFold,
                   nValRuns = nValRuns,
                   useRfast = useRfast,
                   keepTerms = keepTerms,
                   initTime = Sys.time(),
                   save = save,
                   lowerLimit = lowerLimit,
                   filename = filename,
                   filepath = filepath,
                   optimizeScore = optimizeScore,
                   stratificationVector = stratificationVector,
                   keepValidationObjects = TRUE,
                   validateRegression = ifelse(validate,validateRegression,FALSE),
                   validationMethod = validationMethod,
                   validationObject = validationObject,
                   validationParticipants = validationParticipants,
                   ALASCA.version = printVer(get = "version"),
                   ALASCA.version.date = printVer(get = "date")
    )
    printVer()
  }
  class(object) <- "ALASCA"
  if(object$doDebug){
    cat(".. Has initialized the ALASCA model. Next step is to clean it and check input\n")
  }

  object <- sanitizeObject(object)
  
  if(object$doDebug){
    cat(".. Has cleaned the ALASCA model. Next step is building it\n")
  }
  object <- buildModel(object)
  
  if(object$minimizeObject){
    # To save space, we remove unnecessary embedded data
    object <- removeEmbedded(object)
  }
  if(object$validate){
    if(object$doDebug){
      cat(".. You chose to validate the model. Starting validation\n")
    }
    object <- validate(object)
  }
  if(object$save){
    saveALASCAModel(object)
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
RMASCA <- function(...){
  object <- ALASCA(..., method = "LMM")
  return(object)
}

#' Print ALASCA version
#'
#' @return String
#' @export
printVer <- function(object = FALSE, get = NA, print = TRUE){
  ALASCA.version <- "0.0.0.98"
  ALASCA.version.date <- "2022-01-04"
  if(is.list(object)){
    ALASCA.version <- object$ALASCA.version
    ALASCA.version.date <- object$ALASCA.version.date
  }
  if(is.na(get)){
    if(print){
      cat("\n\n====== ALASCA ======\n\n")
      cat(paste0(ALASCA.version, " (", ALASCA.version.date, ")\n\n"))
    }
  }else{
    if(get == "date"){
      return(ALASCA.version.date)
    }else{
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
SMASCA <- function(...){
  object <- ALASCA(..., method = "LM")
  return(object)
}

#' Sanitize an ALASCA object
#'
#' This function checks that the input to an ALASCA object is as expected
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
sanitizeObject <- function(object){
  
  if(!object$minimizeObject){
    object$valCol <- as.character(object$formula)[2]
    # Check formula from user
    formulaTerms <- colnames(attr(terms.formula(object$formula),"factors"))
    object$formulaTerms <- formulaTerms
    if(!is.na(object$method)){ 
      # The user has specified a method to use
      if(object$method == "LMM"){
        if(!any(grepl("\\|",formulaTerms))){
          stop("The model must contain at least one random effect. Sure you wanted linear mixed models?")
        }
      }else if(object$method == "LM"){
        if(any(grepl("\\|",formulaTerms))){
          stop("The model contains at least one random effect. Sure you not wanted linear mixed models instead?")
        }
      }else if(object$method %in% c("KM", "KMM")){
        
      }else{
        stop("You entered an undefined method. Use `LMM` or `LM`!")
      }
      if(object$useRfast){
        cat("Will use Rfast!\n")
        
        # Validation of regression only works for LMs at the moment
        object$validateRegression <- FALSE
      }
    }else{
      # Find which default method to use
      if(any(grepl("\\|",formulaTerms))){
        object$method <- "LMM"
        cat("Will use linear mixed models!\n")
      }else{
        object$method <- "LM"
        cat("Will use linear models!\n")
      }
    }
    
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
  
    if(all(is.na(object$stratificationVector))){
      cat("Using group for stratification.\n")
      object$stratificationVector <- object$df$group
    }
    if(object$method %in% c("KM", "KMM")){
      if(!is.data.frame(object$lowerLimit) | any(!(unique(object$df$variable) %in% unique(object$lowerLimit$variable) ))){
        stop("Some lower limits are not specified!")
      }
      for(i in unique(object$df$variable)){
        maxValue <- max(object$df$value[object$df$variable == i])
        object$df$value[object$df$variable == i] <- maxValue-object$df$value[object$df$variable == i]
        object$lowerLimit$value[object$lowerLimit$variable == i] <- maxValue-object$lowerLimit$value[object$lowerLimit$variable == i]
        object$df$belowLowerLimit[object$df$variable == i] <- object$df$value[object$df$variable == i] > object$lowerLimit$value[object$lowerLimit$variable == i]
        object$df$value[object$df$variable == i & object$df$belowLowerLimit[object$df$variable == i]] <- object$lowerLimit$value[object$lowerLimit$variable == i]
      }
    }
    
    # Change value column if necessary
    if(as.character(object$formula)[2] != "value"){
      cat("Changing",as.character(object$formula)[2],"to `value`.\n")
      object$df$value <- object$df[, get(object$valCol)]
      object$formula <- formula(paste("value ~",
                                      as.character(object$formula)[3]))
    }
  
    if(object$method %in% c("LMM")){
      if(object$participantColumn != "ID"){
        object$df$ID <- object$df[, get(object$participantColumn)]
        tmp <- formulaTerms[!grepl("\\|",formulaTerms)]
        object$formula <- formula(paste("value ~",
                                        paste(tmp, collapse = "+")))
      }else if(any(grepl("\\|",formulaTerms))){
        if(sum(grepl("\\|",formulaTerms)) > 1){
          stop("Multiple random effects, couldn't determine participant-id. Please specify `participantColumn`")
        }else{
          tmp <- formulaTerms[grepl("\\|",formulaTerms)]
          tmp <- gsub(" ", "",tmp)
          tmp <- strsplit(tmp, "\\|")
          object$participantColumn <- tmp[[1]][2]
          object$df$ID <- object$df[, get(object$participantColumn)]
          tmp <- formulaTerms[!grepl("\\|",formulaTerms)]
          object$formula <- formula(paste("value ~",
                                          paste(tmp, collapse = "+")))
        }
      }else if(object$participantColumn == "ID"){
        if(!("ID" %in% colnames(object$df))){
          stop("Please specify participant-id in `participantColumn`")
        }
      }
      
      if(object$method == "LMM"){
        if(object$useRfast){
          # Using Rfast
          rterms <- formulaTerms[!grepl("\\|",formulaTerms)]
          object$newformula <- formula(paste("value ~ ", paste(rterms, collapse = "+")))
        }else{
          #Using lme4
          rterms <- formulaTerms[grepl("\\|",formulaTerms)]
          rterms <- paste0("(",rterms, ")")
          object$newformula <- formula(paste("value ~ modmat+", paste(rterms, collapse = "+")))
        }
      }
    }else if(object$method %in% c("LM")){
      object$newformula <- value ~ modmat
    }
  }
  
  
  if(object$minimizeObject){
    # This is usually a validation object
    object$df <- object$df[object$validationParticipants]
    object$df$ID <- factor(object$df$ID)
  }
  
  if(any(is.na(object$filepath))){
    object$filepath <- paste0("ALASCA/",strftime(object$initTime, format = "%Y%m%d_%H%M%S"),"/")
  }
  
  if(any(is.na(object$filename))){
    object$filename <- "ALASCA"
  }
  
  if(object$doDebug){
    cat(".... Making factors for time, group and variable\n")
  }
  if(!is.factor(object$df$time)){
    object$df$time <- factor(object$df$time)
  }
  if(!is.factor(object$df$group)){
    object$df$group <- factor(object$df$group)
  }
  if(!is.factor(object$df$variable)){
    object$df$variable <- factor(object$df$variable)
  }
  
  if(object$doDebug){
    cat(".... Checking for missing information\n")
  }
  if(any(is.na(object$df$time) | is.na(object$df$group))){
    warning("\n\n!!! -> Oh dear, at least on of your rows is missing either time or group. I have removed it/them for now, but you should check if this is a problem...\n\n")
    object$df <- object$df[!is.na(time) & !is.na(group)]
  }
  
  # Use sum coding?
  if(object$useSumCoding){
    if(object$doDebug){
      cat(".... Use sum coding\n")
    }
    contrasts(object$df$group) <- contr.sum(length(unique(object$df$group)))
  }
  
  # Keep a copy of unscaled data
  object$dfRaw <- object$df
  
  if(is.function(object$scaleFun)){
    if(!object$minimizeObject){
      cat("Scaling data with custom function...\n")
    }
    object$df <- object$scaleFun(object$df)
    object$df$value <- object$df[, get(object$valCol)]
  }else if(object$scaleFun == TRUE){
    if(!object$minimizeObject){
      cat("Scaling data...\n")
    }
    object$df[,value:=as.double(value)][, value := scale(value), by = variable]
  }else{
    if(!object$minimizeObject){
      warning("Not scaling data...\n")
    }
  }
  if(!object$minimizeObject){
    # Check what terms that is present in formula
    object$hasGroupTerm <- ifelse(any(formulaTerms == "group"), TRUE, FALSE)
    object$hasInteractionTerm <- ifelse(any(formulaTerms == "group:time" | formulaTerms == "time:group"), TRUE, FALSE)
    object$covars <- formulaTerms[!(formulaTerms %in% c("time","group","group:time","time:group",object$keepTerms))]
    if(object$doDebug){
      cat(".... Group term in formula? ",object$hasGroupTerm,"\n")
      cat(".... Interaction term in formula? ",object$hasInteractionTerm,"\n")
      cat(".... Identified the following covariates in addition to time and troup: ",object$covars,"\n")
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
removeEmbedded <- function(object){
  object$df <- NULL
  object$dfRaw <- NULL
  object$parts <- NULL
  object$partsWithVariable <- NULL
  object$validationObject <- NULL
  object$regr.model <- NULL
  object$RegressionCoefficients <- NULL
  object$effect.matrix <- NULL
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
flipIt <- function(object, component = NA, effect = "both"){
  if(any(is.na(component))){
    component <- unique(object$ALASCA$score$time$PC)
  }
  for(i in component){
    if(effect %in% c("both", "time")){
      object$ALASCA$score$time$score[object$ALASCA$score$time$PC == i] <- -object$ALASCA$score$time$score[object$ALASCA$score$time$PC == i]
      object$ALASCA$loading$time$loading[object$ALASCA$loading$time$PC == i] <- -object$ALASCA$loading$time$loading[object$ALASCA$loading$time$PC == i]
    }
    if(object$separateTimeAndGroup & effect %in% c("both", "group")){
      object$ALASCA$score$group$score[object$ALASCA$score$group$PC == i] <- -object$ALASCA$score$group$score[object$ALASCA$score$group$PC == i]
      object$ALASCA$loading$group$loading[object$ALASCA$loading$group$PC == i] <- -object$ALASCA$loading$group$loading[object$ALASCA$loading$group$PC == i]
    }
    if(object$validate){
      if(effect %in% c("both", "time")){
        object$ALASCA$score$time$high[object$ALASCA$score$time$PC == i] <- -object$ALASCA$score$time$high[object$ALASCA$score$time$PC == i]
        object$ALASCA$loading$time$high[object$ALASCA$loading$time$PC == i] <- -object$ALASCA$loading$time$high[object$ALASCA$loading$time$PC == i]
        object$ALASCA$score$time$low[object$ALASCA$score$time$PC == i] <- -object$ALASCA$score$time$low[object$ALASCA$score$time$PC == i]
        object$ALASCA$loading$time$low[object$ALASCA$loading$time$PC == i] <- -object$ALASCA$loading$time$low[object$ALASCA$loading$time$PC == i]
      }
      if(object$separateTimeAndGroup & effect %in% c("both", "group")){
        object$ALASCA$score$group$high[object$ALASCA$score$group$PC == i] <- -object$ALASCA$score$group$high[object$ALASCA$score$group$PC == i]
        object$ALASCA$loading$group$high[object$ALASCA$loading$group$PC == i] <- -object$ALASCA$loading$group$high[object$ALASCA$loading$group$PC == i]
        object$ALASCA$score$group$low[object$ALASCA$score$group$PC == i] <- -object$ALASCA$score$group$low[object$ALASCA$score$group$PC == i]
        object$ALASCA$loading$group$low[object$ALASCA$loading$group$PC == i] <- -object$ALASCA$loading$group$low[object$ALASCA$loading$group$PC == i]
      }
    }
  }
  
  if(object$validate == TRUE){
    for(i in seq_along(object$validation$temp_objects)){
      object$validation$temp_objects[[i]] <- flipIt(object$validation$temp_objects[[i]], component = component, effect = effect)
    }
  }
  return(object)
}

#' Summary
#'
#' Gives some general information
#'
#' @inheritParams ALASCA
#' @export
summary.ALASCA <- function(object){
  cat("================ ALASCA ================\n")
  cat("Model initialized ", as.character(object$initTime), " using ",object$method," on ",length(unique(mod$RegressionCoefficients$covar))," variables. ", sep = "")
  if(object$validate){
    cat("The model been validated.\n")
  }else{
    cat("The model has *not* been validated yet.\n")
  }
  cat("Terms in model:\n   * ",paste(names(lme4::fixef(object$regr.model[[1]])), collapse = "\n   * "),"\n", sep = "")
  cat("\nPCs explaining at least 5% of variation:\n   Time: ", 
      paste(getRelevantPCs(object$pca$score$explained$time), collapse = ", "), " (",
      paste(round(100*object$pca$score$explained$time[getRelevantPCs(object$pca$score$explained$time)], 2), collapse = "%, "), "%)", 
         ifelse(object$separateTimeAndGroup,paste0("\n   Group: ",
                                                   paste(getRelevantPCs(object$pca$score$explained$group), collapse = ", "), " (",
                                                   paste(round(100*object$pca$score$explained$group[getRelevantPCs(object$pca$score$explained$group)], 2), collapse = "%, "), "%)"),"\n"), sep = "")
  cat("\nNumber of data points (based on ",names(object$regr.model)[[1]]," measurements",ifelse(object$missingMeasurements,", some variables have more/fewer observations",""),"):\n", sep = "")
  aggregate(data = subset(object$df, variable == names(object$regr.model)[[1]]), as.formula(paste(as.character(object$formula)[[2]], " ~ group + time")), FUN = length)
  if(is.function(object$scaleFun)){
    cat("\nScaling function:\n")
    object$scaleFun
  }else{
    cat("\nNo scaling performed.\n")
  }
}

#' Save ALASCA object
#'
#' @param object An ALASCA object
#' @return An ALASCA object
#' @export
saveALASCAModel <- function(object){
  if(!dir.exists(object$filepath)){
    dir.create(object$filepath, recursive = TRUE)
  }
  fname <- paste0(object$filepath,object$filename,".Rdata")
  cnt <- 1
  while(file.exists(fname)){
    fname <- paste0(object$filepath,object$filename,"_",cnt,".Rdata")
    cnt <- cnt + 1
  }
  save(object, file=fname)
  cat(paste0("- Saved model to ", fname,"\n"))
}

#' Append a new model to an existing model
#'
#' 
#'
#' @param target The existing model
#' @param object The new model
#' @param method `rotate` or `project` (default)
#' @return An ALASCA object
#' @export
appendModel <- function(target, object, method = "project"){
  if(method == "rotate"){
    target$ALASCA$score$time$model <- "Model 1"
    object$ALASCA$score$time$model <- "Model 2"
    object <- ALASCA:::rotateMatrix(object = object, target = target)
    target$ALASCA$score$time <- rbind(target$ALASCA$score$time, object$ALASCA$score$time)
    target$ALASCA$loading$time$model <- "Model 1"
    object$ALASCA$loading$time$model <- "Model 2"
    target$ALASCA$loading$time <- rbind(target$ALASCA$loading$time, object$ALASCA$loading$time)
    if(target$separateTimeAndGroup){
      target$ALASCA$score$group$model <- "Model 1"
      object$ALASCA$score$group$model <- "Model 2"
      target$ALASCA$score$group <- rbind(target$ALASCA$score$group, object$ALASCA$score$group)
      target$ALASCA$loading$group$model <- "Model 1"
      object$ALASCA$loading$group$model <- "Model 2"
      target$ALASCA$loading$group <- rbind(target$ALASCA$loading$group, object$ALASCA$loading$group)
    }
  }else{
    tmp_pca <- as.data.frame(
      scale(object$effect.matrix[object$effect.matrix$comp == "TIME",1:(ncol(object$effect.matrix)-1)],
            target$pca$time$center,
            target$pca$time$scale) %*% target$pca$time$rotation
      )
    tmp_pca$group <- object$parts$group
    tmp_pca$time <- object$parts$time
    tmp_pca <- tmp_pca[!duplicated(tmp_pca),]
    tmp_pca <- reshape2::melt(tmp_pca, id.vars = c("time","group"))
    colnames(tmp_pca) <- c("time", "group", "PC", "score")
    tmp_pca$PC <- as.numeric(gsub("PC","",tmp_pca$PC))
    target$ALASCA$score$time <- rbind(target$ALASCA$score$time, tmp_pca)
    
    # temp_df <- merge(object$mod.pred, target$ALASCA$loading$time, by.y = "covars", by.x = "variable", all = TRUE)
    # temp_df$score <- temp_df$loading * temp_df$pred
    # temp_df <- aggregate(data = temp_df, score ~ time + group + PC, FUN = sum)
    # target$ALASCA$score$time <- rbind(target$ALASCA$score$time, temp_df)
    # if(target$separateTimeAndGroup){
    #   temp_df <- merge(object$mod.pred, target$ALASCA$loading$group, by.y = "covars", by.x = "variable", all = TRUE)
    #   temp_df$score <- temp_df$loading * temp_df$pred
    #   temp_df <- aggregate(data = temp_df, score ~ time + group + PC, FUN = sum)
    #   target$ALASCA$score$group <- rbind(target$ALASCA$score$group, temp_df)
    # }
  }
  
  return(target)
}