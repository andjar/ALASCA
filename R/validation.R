#' Validate the ALASCA model LMM
#'
#' This function performs leave-one-out robustness testing of your ALASCA model. If you didn't specify the number of runs `nValRuns` when initializing the model (see \code{\link{ALASCA}}), you can do it by running for example `model$nValRuns <- 100` prior to calling `validate`. Your dataset is divided into `nValFold` partitions, keeping group proportions, and one of these are left out. `nValFold` is set the same way as  `nValRuns`.
#'
#' @param object An ALASCA object
#' @param participantColumn The name of the column containing participant identifier. Needed if not set during initialization of the model.
#' @param validateRegression Whether to validate regression models
#' @return An ALASCA object
#' 
#' @examples
#' load("PE.Rdata")
#' model$nValRuns <- 10
#' model.val <- validate(model, participantColumn = "ID")
#' 
#' @export
validate <- function(object, participantColumn = FALSE, validateRegression = FALSE){
  if(object$validate){
    #stop("The object has already been validated")
  }
  object$validate <- TRUE
  if(validateRegression){
    object$validateRegression <- TRUE
  }

  if(object$method == "LMM"){
    if(participantColumn == FALSE){
      if(object$participantColumn == FALSE){
        stop("You need to specify the column containing participant id in `participantColumn`")
      }
    }else{
      object$participantColumn <- participantColumn
    }
  }
  
  cat("Running validation...\n")

  start.time.all <- Sys.time()

  if(any(is.na(object$validationIDs))){
    # Generate random samples
    if(object$savetodisk){
      limPC_time <- getRelevantPCs(object = object,  object$ALASCA$loading$explained$time)
      if(object$separateTimeAndGroup){
        limPC_group <- getRelevantPCs(object = object,  object$ALASCA$loading$explained$group)
      }
      temp_object <- lapply(1:object$nValRuns, FUN = function(ii){
        cat("- Run ",ii," of ",object$nValRuns,"\n")
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- prepareValidationRun(object)
        
        # Rotate new loadings/scores to the original model
        if(object$optimizeScore){
          temp_object <- rotateMatrixOptimizeScore(object = temp_object, target = object)
        }else{
          temp_object <- rotateMatrix(object = temp_object, target = object)
        }
        
        temp_object <- cleanALASCA(temp_object)
        
        # Save to disk
        fname <- paste0("val_",ii)
        loading <- data.frame(temp_object$ALASCA$loading$time[temp_object$ALASCA$loading$time$PC %in% limPC_time,], model = ii)
        #loading$covars <- object$numvariablelist[loading$covars]
        DBI::dbWriteTable(object$db.con, "time.loading", loading, append= T, overwrite=F)
        scores <- data.frame(temp_object$ALASCA$score$time[temp_object$ALASCA$score$time$PC %in% limPC_time,], model = ii)
        #scores$time <- object$numtimelist[scores$time]
        #if(!object$separateTimeAndGroup){
        #  scores$group <- object$numgrouplist[scores$group]
        #}
        DBI::dbWriteTable(object$db.con, "time.score", scores, append= T, overwrite=F)
        
        if(object$separateTimeAndGroup){
          loading <- data.frame(temp_object$ALASCA$loading$group[temp_object$ALASCA$loading$group$PC %in% limPC_group,], model = ii)
          #loading$covars <- object$numvariablelist[loading$covars]
          DBI::dbWriteTable(object$db.con, "group.loading", loading, append= T, overwrite=F)
          
          scores <- data.frame(temp_object$ALASCA$score$group[temp_object$ALASCA$score$group$PC %in% limPC_group,], model = ii)
          #scores$time <- object$numtimelist[scores$time]
          #scores$group <- object$numgrouplist[scores$group]
          DBI::dbWriteTable(object$db.con, "group.score", scores, append= T, overwrite=F)
        }
        if(object$validateRegression){
          DBI::dbWriteTable(object$db.con, "mod.pred", data.frame(temp_object$mod.pred, model = ii), append= T, overwrite=F)
        }
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs"))/ii
        cat("--- Used ",round(difftime(Sys.time(), start.time.this, units = c("secs")),2)," seconds. Est. time remaining: ",round((object$nValRuns-ii)*time_all,2)," seconds \n")
        fname
      })
    }else{
      temp_object <- lapply(1:object$nValRuns, FUN = function(ii){
        cat("- Run ",ii," of ",object$nValRuns,"\n")
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- prepareValidationRun(object)
        
        # Rotate new loadings/scores to the original model
        if(object$optimizeScore){
          temp_object <- rotateMatrixOptimizeScore(object = temp_object, target = object)
        }else{
          temp_object <- rotateMatrix(object = temp_object, target = object)
        }
        temp_object <- cleanALASCA(temp_object)
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs"))/ii
        cat("--- Used ",round(difftime(Sys.time(), start.time.this, units = c("secs")),2)," seconds. Est. time remaining: ",round((object$nValRuns-ii)*time_all,2)," seconds \n")
        temp_object
      })
    }
  }else{
    # Reuse previous samples
    cat("Using predefined samples...\n")
    
    if(object$savetodisk){
      limPC_time <- getRelevantPCs(object = object,  object$ALASCA$loading$explained$time)
      if(object$separateTimeAndGroup){
        limPC_group <- getRelevantPCs(object = object,  object$ALASCA$loading$explained$group)
      }
      temp_object <- lapply(1:object$nValRuns, FUN = function(ii){
        cat("- Run ",ii," of ",object$nValRuns,"\n")
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- prepareValidationRun(object, runN = ii)
        
        # Rotate new loadings/scores to the original model
        if(object$optimizeScore){
          temp_object <- rotateMatrixOptimizeScore(object = temp_object, target = object)
        }else{
          temp_object <- rotateMatrix(object = temp_object, target = object)
        }
        
        temp_object <- cleanALASCA(temp_object)
        
        # Save to disk
        fname <- paste0("val_",ii)
        loading <- data.frame(temp_object$ALASCA$loading$time[temp_object$ALASCA$loading$time$PC %in% limPC_time,], model = ii)
        #loading$covars <- object$numvariablelist[loading$covars]
        DBI::dbWriteTable(object$db.con, "time.loading", loading, append= T, overwrite=F)
        scores <- data.frame(temp_object$ALASCA$score$time[temp_object$ALASCA$score$time$PC %in% limPC_time,], model = ii)
        #scores$time <- object$numtimelist[scores$time]
        #if(!object$separateTimeAndGroup){
        #  scores$group <- object$numgrouplist[scores$group]
        #}
        DBI::dbWriteTable(object$db.con, "time.score", scores, append= T, overwrite=F)
        
        if(object$separateTimeAndGroup){
          loading <- data.frame(temp_object$ALASCA$loading$group[temp_object$ALASCA$loading$group$PC %in% limPC_group,], model = ii)
          #loading$covars <- object$numvariablelist[loading$covars]
          DBI::dbWriteTable(object$db.con, "group.loading", loading, append= T, overwrite=F)
          
          scores <- data.frame(temp_object$ALASCA$score$group[temp_object$ALASCA$score$group$PC %in% limPC_group,], model = ii)
          #scores$time <- object$numtimelist[scores$time]
          #scores$group <- object$numgrouplist[scores$group]
          DBI::dbWriteTable(object$db.con, "group.score", scores, append= T, overwrite=F)
        }
        if(object$validateRegression){
          DBI::dbWriteTable(object$db.con, "mod.pred", data.frame(temp_object$mod.pred, model = ii), append= T, overwrite=F)
        }
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs"))/ii
        cat("--- Used ",round(difftime(Sys.time(), start.time.this, units = c("secs")),2)," seconds. Est. time remaining: ",round((object$nValRuns-ii)*time_all,2)," seconds \n")
        fname
      })
    }else{
      temp_object <- lapply(1:object$nValRuns, FUN = function(ii){
        cat("- Run ",ii," of ",object$nValRuns,"\n")
        start.time.this <- Sys.time()
        
        # Make resampled model
        temp_object <- prepareValidationRun(object, runN = ii)
        
        # Rotate new loadings/scores to the original model
        if(object$optimizeScore){
          temp_object <- rotateMatrixOptimizeScore(object = temp_object, target = object)
        }else{
          temp_object <- rotateMatrix(object = temp_object, target = object)
        }
        temp_object <- cleanALASCA(temp_object)
        
        time_all <- difftime(Sys.time(), start.time.all, units = c("secs"))/ii
        cat("--- Used ",round(difftime(Sys.time(), start.time.this, units = c("secs")),2)," seconds. Est. time remaining: ",round((object$nValRuns-ii)*time_all,2)," seconds \n")
        temp_object
      })
    }
  }
  
  if(grepl("permutation",object$validationMethod)){
    cat("- Calculates P values...\n")
    object <- getPermutationPValues(object, objectlist = temp_object)
  }
  
  cat("- Calculates percentiles for score and loading...\n")
  object <- getValidationPercentiles(object, objectlist = temp_object)
  
  if(object$keepValidationObjects){
    object$validation$temp_objects <- temp_object
  }
  return(object)
}

.procrustes <- function(loadings, target){
  s= t(loadings)%*%target
  w1 = s %*% t(s)
  v1 = t(s) %*% s
  w <- eigen(w1) $vectors
  ew <- diag(eigen(w1) $values)
  v <- eigen(v1) $vectors
  ev <- diag(eigen(v1) $values)
  o = t(w) %*% s %*% v
  k = diag(  ((diag(o)) / abs(diag(o))) , nrow = nrow(o), ncol =nrow(o))
  ww = w %*% k
  out <- list()
  out$t1 = ww %*% t(v) # Rotation matrix
  out$procrust = loadings %*% out$t1 # Rotated loadings
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
rotateMatrixOptimizeScore <- function(object, target){
  # We are only looking at components explaining more than a predefined value
  PCloading <- getRelevantPCs(object = target, target$ALASCA$loading$explained$time)
  
  # PCA can give loadings with either sign, so we have to check whether this improves the rotation
  N   <- length(PCloading)
  vec <- c(-1, 1)
  lst <- lapply(numeric(N), function(x) vec)
  signMatrix <- as.matrix(expand.grid(lst))
  signVar <- Reduce(cbind,lapply(1:nrow(signMatrix), function(i){
    c <- .procrustes(loadings= as.matrix(t(t(object$pca$loading$time[,PCloading]) * signMatrix[i,])),
                    target = as.matrix(target$pca$loading$time[,PCloading]))
    (target$pca$score$time[,PCloading] - as.matrix(t(t(object$pca$score$time[,PCloading]) * signMatrix[i,])) %*% solve(c$t1) )^2
  }))
  minSignVar <- which(signVar == min(signVar))[1]
  object$pca$loading$time[,PCloading] <- t(t(object$pca$loading$time[,PCloading]) * signMatrix[minSignVar,])
  object$pca$score$time[,PCloading] <- t(t(object$pca$score$time[,PCloading]) * signMatrix[minSignVar,])

  c <- .procrustes(loadings= as.matrix(object$pca$loading$time[,PCloading]),
                  target = as.matrix(target$pca$loading$time[,PCloading]))
  object$pca$loading$time[,PCloading] <- c$procrust
  object$pca$score$time[,PCloading] <- as.matrix(object$pca$score$time[,PCloading]) %*% solve(c$t1)
  
  if(object$separateTimeAndGroup){
    # We are only looking at components explaining more than 5% of variation
    PCloading <- getRelevantPCs(object = target, target$ALASCA$loading$explained$group)
    
    # PCA can give loadings with either sign, so we have to check whether this improves the rotation
    N   <- length(PCloading)
    vec <- c(-1, 1)
    lst <- lapply(numeric(N), function(x) vec)
    signMatrix <- as.matrix(expand.grid(lst))
    signVar <- Reduce(cbind,lapply(1:nrow(signMatrix), function(i){
      c <- .procrustes(loadings= as.matrix(t(t(object$pca$loading$group[,PCloading]) * signMatrix[i,])),
                      target = as.matrix(target$pca$loading$group[,PCloading]))
      (target$pca$score$group[,PCloading] - as.matrix(t(t(object$pca$score$group[,PCloading]) * signMatrix[i,])) %*% solve(c$t1) )^2
    }))
    minSignVar <- which(signVar == min(signVar))[1]
    object$pca$loading$group[,PCloading] <- t(t(object$pca$loading$group[,PCloading]) * signMatrix[minSignVar,])
    object$pca$score$group[,PCloading] <- t(t(object$pca$score$group[,PCloading]) * signMatrix[minSignVar,])

    c <- .procrustes(loadings= as.matrix(object$pca$loading$group[,PCloading]),
                    target = as.matrix(target$pca$loading$group[,PCloading]))
    object$pca$loading$group[,PCloading] <- c$procrust
    object$pca$score$group[,PCloading] <- as.matrix(object$pca$score$group[,PCloading]) %*% solve(c$t1)
  }
  
  return(object)
}

#' Rotate PCA
#'
#' This function rotates loadings and scores during validation
#'
#' @param object ALASCA object to be rotated (and returned)
#' @param target ALASCA object acting as target
#' @return An ALASCA object
rotateMatrix <- function(object, target){
  PCloading <- getRelevantPCs(object = target, target$ALASCA$loading$explained$time)
  c <- .procrustes(loadings= as.matrix(object$pca$loading$time[,PCloading]),
                  target = as.matrix(target$pca$loading$time[,PCloading]))
  
  object$pca$loading$time[,PCloading] <- c$procrust
  object$pca$score$time[,PCloading] <- as.matrix(object$pca$score$time[,PCloading]) %*% solve(c$t1)
  
   if(object$separateTimeAndGroup){
    PCloading <- getRelevantPCs(object = target, target$ALASCA$loading$explained$group)
    c <- .procrustes(loadings= as.matrix(object$pca$loading$group[,PCloading]),
                    target = as.matrix(target$pca$loading$group[,PCloading]))
    
    object$pca$loading$group[,PCloading] <- c$procrust
    object$pca$score$group[,PCloading] <- as.matrix(object$pca$score$group[,PCloading]) %*% solve(c$t1)
  }
  
  return(object)
}

#' Extract percentiles
#'
#' This function extract percentiles during validation
#'
#' @param object ALASCA object
#' @param objectlist List of ALASCA objects
#' @return An ALASCA object
getValidationPercentiles <- function(object, objectlist){

  if("low" %in% colnames(object$ALASCA$loading$time)){
    object$ALASCA$loading$time$low <- NULL
    object$ALASCA$loading$time$high <- NULL
    object$ALASCA$score$time$low <- NULL
    object$ALASCA$score$time$high <- NULL
    if(object$separateTimeAndGroup){
      object$ALASCA$loading$group$low <- NULL
      object$ALASCA$loading$group$high <- NULL
      object$ALASCA$score$group$low <- NULL
      object$ALASCA$score$group$high <- NULL
    }
  }
  
  object <- getValidationPercentilesLoading(object, objectlist)
  object <- getValidationPercentilesScore(object, objectlist)
  if(object$validateRegression){
    object <- getValidationPercentilesRegression(object, objectlist)
  }
  if(nrow(getCovars(object)) > 0){
    object <- getValidationPercentilesCovars(object, objectlist)
  }
  
  return(object)
}

#' Extract percentiles for regressions
#'
#' This function extract percentiles for validation of regression
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesRegression <- function(object, objectlist){
  if(object$savetodisk){
    res <- DBI::dbSendQuery(object$db.con, "SELECT * FROM 'mod.pred'")
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  }else{
    df <- rbindlist(lapply(objectlist, function(x) x$mod.pred))
  }
  df <- df[, as.list(quantile(pred, probs = object$limitsCI, type = 1)),by = .(group, time, variable)]
  colnames(df) = c("group", "time", "variable", "low", "high")
  
  object$mod.pred <- merge(object$mod.pred, df)
  return(object)
}

#' Extract percentiles for covariates
#'
#' This function extract percentiles for validation of covariates
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesCovars <- function(object, objectlist){
  if(object$savetodisk){
    stop("Not working yet")
  }else{
    df <- rbindlist(lapply(objectlist, function(x) getCovars(x)))
  }
  df <- df[, as.list(quantile(estimate, probs = object$limitsCI, type = 1)),by = .(covar, variable)]
  colnames(df) = c("covar", "variable", "low", "high")
  
  object$CovarCoefficients <- merge(object$CovarCoefficients, df)
  return(object)
}

#' Extract percentiles for loading
#'
#' This function extract percentiles during validation of loadings
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesLoading <- function(object, objectlist){
  PC_time <- getRelevantPCs(object = object,  object$ALASCA$loading$explained$time)
  if(object$savetodisk){
    res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.loading' WHERE PC IN(",paste(PC_time, collapse = ', '),")"))
    df_time <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  }else{
    df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$loading$time[x$ALASCA$loading$time$PC %in% PC_time,]), fill = TRUE)
  }
  
  object$validation$time$loading <- df_time[, as.list(quantile(loading, probs = object$limitsCI, type = 1)),by = .(PC, covars)]
  colnames(object$validation$time$loading) = c("PC", "covars", "low", "high")
  object$ALASCA$loading$time <- merge(object$ALASCA$loading$time, object$validation$time$loading, all.x = TRUE)
  
  if(object$separateTimeAndGroup){
    PC_group <- getRelevantPCs(object = object,  object$ALASCA$loading$explained$group)
    if(object$savetodisk){
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'group.loading' WHERE PC IN(",paste(PC_group, collapse = ', '),")"))
      df_group <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    }else{
      df_group <- rbindlist(lapply(objectlist, function(x) x$ALASCA$loading$group[x$ALASCA$loading$group$PC %in% PC_group,]), fill = TRUE)
    }
    
    object$validation$group$loading <- df_group[, as.list(quantile(loading, probs = object$limitsCI, type = 1)),by = .(PC, covars)]
    colnames(object$validation$group$loading) = c("PC", "covars", "low", "high")
    object$ALASCA$loading$group <- merge(object$ALASCA$loading$group, object$validation$group$loading, all.x = TRUE)
  }
  return(object)
}

#' Extract percentiles for score
#'
#' This function extract percentiles during validation of scores
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesScore <- function(object, objectlist){
  if(object$separateTimeAndGroup){
    # Separate time and group effects
    
    PC_time <- getRelevantPCs(object = object,  object$ALASCA$score$explained$time)
    if(object$savetodisk){
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(",paste(PC_time, collapse = ', '),")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    }else{
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$time[x$ALASCA$score$time$PC %in% PC_time,]), fill = TRUE)
    }
    
    object$validation$time$score <- df_time[, as.list(quantile(score, probs = object$limitsCI, type = 1)),by = .(PC, time)]
    colnames(object$validation$time$score) = c("PC", "time", "low", "high")
    object$ALASCA$score$time <- merge(object$ALASCA$score$time, object$validation$time$score, all.x = TRUE)
    
    PC_group <- getRelevantPCs(object = object,  object$ALASCA$score$explained$group)
    if(object$savetodisk){
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'group.score' WHERE PC IN(",paste(PC_group, collapse = ', '),")"))
      df_group <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    }else{
      df_group <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$group[x$ALASCA$score$group$PC %in% PC_group,]), fill = TRUE)
    }
    
    object$validation$group$score <- df_group[, as.list(quantile(score, probs = object$limitsCI, type = 1)),by = .(PC, time, group)]
    colnames(object$validation$group$score) = c("PC", "time", "group", "low", "high")
    object$ALASCA$score$group <- merge(object$ALASCA$score$group, object$validation$group$score, all.x = TRUE)
  }else{
    # Pooled time and groups effects
    PC_time <- getRelevantPCs(object = object,  object$ALASCA$score$explained$time)
    if(object$savetodisk){
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(",paste(PC_time, collapse = ', '),")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    }else{
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$time[x$ALASCA$score$time$PC %in% PC_time,]), fill = TRUE)
    }
    
    object$validation$time$score <- df_time[, as.list(quantile(score, probs = object$limitsCI, type = 1)),by = .(PC, time, group)]
    colnames(object$validation$time$score) = c("PC", "time", "group", "low", "high")
    object$ALASCA$score$time <- merge(object$ALASCA$score$time, object$validation$time$score, all.x = TRUE)
  }
  return(object)
}

#' Get relevant PCs
#'
#' This function extract percentiles during validation
#'
#' @param x Explanatory power of PC
#' @return A vector with relevant PCs
getRelevantPCs <- function(object, x){
  PC <- x >= object$exploratorylimit
  PC[1:2] <- TRUE
  return(which(PC))
}

#' Validate underlying regression models
#'
#' This function calcuates predictions from each regression model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
getRegressionPredictions <- function(object){
  if(!object$minimizeObject){
    # This is not a validation run
    cat("Calculating predictions from regression models...\n")
  }
    regCoeffAll <- reshape2::dcast(data = object[["RegressionCoefficients"]], covar ~ variable, value.var = "estimate")
    regModel <- unique(model.matrix(object$formula, data = object$df))
    if(object$forceEqualBaseline){
      regModel <- regModel[,!grepl(paste0("time",levels(object$df$time)[1]), colnames(regModel))]
    }
    newdata <- data.table::rbindlist(
      lapply(object$variablelist, function(x){
        regCoeff <- as.matrix(regCoeffAll[regCoeffAll$covar == x, -1])
        data.frame(
          time = object$df$time[as.numeric(rownames(regModel))],
          group = object$df$group[as.numeric(rownames(regModel))],
          pred = colSums(regCoeff[,colnames(regModel)] * t(regModel)),
          variable = x
        )
      })
    )
    object$mod.pred <- newdata[, .(pred = mean(pred)), by = c("variable", "time", "group")]
  
  if(!object$minimizeObject){
    # This is not a validation run
    cat("Finished calculating predictions from regression models!\n")
  }
  return(object)
}

#' Make a single validation run
#'
#' This function ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
prepareValidationRun <- function(object, runN = NA){
  if(object$validationMethod %in% c("loo", "jack-knife", "jackknife")){
    # Use leave-one-out validation
    selectedParts <- data.frame()
    
    if(object$method %in% c("LMM", "Rfast")){
      if(any(is.na(object$validationIDs))){
        # For each group, divide the participants into nValFold groups, and select nValFold-1 of them
        selectedParts <- lapply(unique(object$stratificationVector), function(gr){
          selectedParts_temp_all <- unique(object$df[object$stratificationVector == gr,ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$nValFold
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
        temp_object <- ALASCA(validationObject = object,
                              validationParticipants = object$df[,ID] %in% unlist(selectedParts))
      }else{
        temp_object <- ALASCA(validationObject = object,
                              validationParticipants = object$df[,ID] %in% object$validationIDs[runN,x])
      }
      
    }else if(object$method == "LM"){
      object$df$ID <- c(1:nrow(object$df))
      if(any(is.na(object$validationIDs))){
        # For each group, divide the participants into nValFold groups, and select nValFold-1 of them
        selectedParts <- lapply(unique(object$stratificationVector), function(gr){
          selectedParts_temp_all <- unique(object$df[object$stratificationVector == gr,ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$nValFold
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
      
      temp_object <- ALASCA(validationObject = object,
                            validationParticipants = object$df[,ID] %in% unlist(selectedParts))
    }else{
      temp_object <- ALASCA(validationObject = object,
                            validationParticipants = object$df[,ID] %in% object$validationIDs[runN,x])
    }
    }
  }else if(object$validationMethod == "bootstrap"){
    # Use bootstrap validation
    # When using bootstrapping, we resample participants with replacement
    bootobject <- object
    bootdf_temp <- object$dfRaw
    bootdf <- data.frame()
    cc_id <- 0 # Will become the new participant ID
    
    if(object$method %in% c("LMM", "Rfast")){
      # Loop through all the groups and create a new dataframe with resampled values
      bootobject$newIDs <- c()
      bootobject$originalIDs <- c()
      if(any(is.na(object$validationIDs))){
        for(i in unique(object$stratificationVector)){
          # Get ID of all members of stratification group
          selectedParts_temp_all <- unique(object$df[object$stratificationVector == i,ID])
          
          # Resample participants
          selectedParts_temp_selected <- sample(selectedParts_temp_all, length(selectedParts_temp_all), replace = TRUE)
          newIDs <- seq(cc_id+1, cc_id+length(selectedParts_temp_selected))
          cc_id <- max(newIDs)
          bootobject$originalIDs <- c(bootobject$originalIDs, selectedParts_temp_selected)
          bootobject$newIDs <- c(bootobject$newIDs, newIDs)
          
          # Create data frame from resampled participants
          bootdf <- rbind(bootdf,
                          data.table::rbindlist(
                            lapply(seq_along(selectedParts_temp_selected), function(x){
                              seldf <- bootdf_temp[bootdf_temp$ID == selectedParts_temp_selected[x], ]
                              seldf$originalIDbeforeBootstrap <- seldf$ID
                              if(object$validationAssignNewID) seldf$ID <- newIDs[x] # Replace ID
                              seldf
                            })
                          )
          )
        }
        
      }else{
        newIDs <- seq(1, length(object$validationIDs[runN,]))
        bootdf <- rbind(bootdf,
                        data.table::rbindlist(
                          lapply(seq_along(object$validationIDs[runN,]), function(x){
                            seldf <- bootdf_temp[bootdf_temp$ID == object$validationIDs[runN,x], ]
                            seldf$originalIDbeforeBootstrap <- seldf$ID
                            if(object$validationAssignNewID) seldf$ID <- newIDs[x] # Replace ID
                            seldf
                          })
                        )
        )
      }
      bootobject$dfRaw <- bootdf
      
      temp_object <- ALASCA(validationObject = bootobject,
                            validationParticipants = rep(TRUE,nrow(bootobject$dfRaw)))
    }else if(df$method == "LM"){
      stop("Bootstrapping not implemented for LMs yet")
    }
  }
  
  return(temp_object)
}