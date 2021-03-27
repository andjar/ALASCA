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
  temp_object <- list()
  time_mean <- c()

  for(ii in 1:object$nValRuns){
    cat("- Run ",ii," of ",object$nValRuns,"\n")
    start.time <- Sys.time()
    
    # Make resampled model
    temp_object[[ii]] <- prepareValidationRun(object)
      
    
    # Rotate new loadings/scores to the original model
    temp_object[[ii]] <- rotateMatrix(object = temp_object[[ii]], target = object)
    temp_object[[ii]] <- cleanALASCA(temp_object[[ii]])
    
    end.time <- Sys.time()
    time_mean[ii] <- end.time - start.time
    cat("--- Used ",round(time_mean[ii],2)," seconds. Est. time remaining: ",round((object$nValRuns-ii)*mean(time_mean),2)," seconds \n")
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

#' Rotate PCA
#'
#' This function rotates loadings and scores during validation
#'
#' @param object ALASCA object to be rotated (and returned)
#' @param target ALASCA object acting as target
#' @return An ALASCA object
rotateMatrix <- function(object, target){
  # We are only looking at components explaining more than 5% of variation
  PCloading <- target$ALASCA$loading$explained$time > 0.05
  PCloading[1:2] <- TRUE
  PCloading <- which(PCloading)
  
  
  a_l <- object$pca$loading
  b_l <- target$pca$loading
  a_s <- object$pca$score
  b_s <- target$pca$score
  
  procrustes <- function(loadings, target){
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
  
  # PCA can give loadings with either sign, so we have to check whether this improves the rotation
  N   <- length(PCloading)
  vec <- c(-1, 1)
  lst <- lapply(numeric(N), function(x) vec)
  signMatrix <- as.matrix(expand.grid(lst))
  signVar <- Reduce(cbind,lapply(1:nrow(signMatrix), function(i){
    c <- procrustes(loadings= as.matrix(t(t(a_l$time[,PCloading]) * signMatrix[i,])),
                    target = as.matrix(b_l$time[,PCloading]))
    ifelse(object$optimizeScore,
           sum((b_s$time[,PCloading] - as.matrix(t(t(a_s$time[,PCloading]) * signMatrix[i,])) %*% solve(c$t1) )^2),
           sum((b_l$time[,PCloading] - c$procrust)^2))
  }))
  minSignVar <- which(signVar == min(signVar))[1]
  object$pca$loading$time[,PCloading] <- t(t(object$pca$loading$time[,PCloading]) * signMatrix[minSignVar,])
  object$pca$score$time[,PCloading] <- t(t(object$pca$score$time[,PCloading]) * signMatrix[minSignVar,])
  a_l <- object$pca$loading

  c <- procrustes(loadings= as.matrix(a_l$time[,PCloading]),
                  target = as.matrix(b_l$time[,PCloading]))
  object$pca$loading$time[,PCloading] <- c$procrust
  object$pca$score$time[,PCloading] <- as.matrix(object$pca$score$time[,PCloading]) %*% solve(c$t1)
  
  if(object$separateTimeAndGroup){
    # We are only looking at components explaining more than 5% of variation
    PCloading <- target$pca$loading$explained$group > 0.05
    PCloading[1:2] <- TRUE
    PCloading <- which(PCloading)
    
    # PCA can give loadings with either sign, so we have to check whether this improves the rotation
    N   <- length(PCloading)
    vec <- c(-1, 1)
    lst <- lapply(numeric(N), function(x) vec)
    signMatrix <- as.matrix(expand.grid(lst))
    signVar <- Reduce(cbind,lapply(1:nrow(signMatrix), function(i){
      c <- procrustes(loadings= as.matrix(t(t(a_l$group[,PCloading]) * signMatrix[i,])),
                      target = as.matrix(b_l$group[,PCloading]))
      ifelse(object$optimizeScore,
             sum((b_s$group[,PCloading] - as.matrix(t(t(a_s$group[,PCloading]) * signMatrix[i,])) %*% solve(c$t1) )^2),
             sum((b_l$group[,PCloading] - c$procrust)^2))
    }))
    minSignVar <- which(signVar == min(signVar))[1]
    object$pca$loading$group[,PCloading] <- t(t(object$pca$loading$group[,PCloading]) * signMatrix[minSignVar,])
    object$pca$score$group[,PCloading] <- t(t(object$pca$score$group[,PCloading]) * signMatrix[minSignVar,])
    a_l <- object$pca$loading

    c <- procrustes(loadings= as.matrix(a_l$group[,PCloading]),
                    target = as.matrix(b_l$group[,PCloading]))
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

  object <- getValidationPercentilesLoading(object, objectlist)
  object <- getValidationPercentilesScore(object, objectlist)
  if(object$validateRegression){
    object <- getValidationPercentilesRegression(object, objectlist)
  }
  
  return(object)
}

#' Extract percentiles for loading
#'
#' This function extract percentiles during validation
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesRegression <- function(object, objectlist){
  df <- Reduce(rbind,lapply(objectlist, function(x) x$mod.pred))
  df <- aggregate(data = df, pred ~ group + time + variable, FUN = function(x) quantile(x , probs = c(0.025, 0.975) ))
  df$low <- df$pred[,1]
  df$high <- df$pred[,2]
  df$pred <- NULL
  colnames(df) <- c("group", "time", "variable", "low", "high")
  object$mod.pred <- merge(object$mod.pred, df)
  return(object)
}

#' Extract percentiles for loading
#'
#' This function extract percentiles during validation
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesLoading <- function(object, objectlist){
  df_time <- Reduce(rbind,lapply(objectlist, function(x) x$pca$loading$time))
  PC_time <- getRelevantPCs(object$pca$loading$explained$time)
  perc_time <- aggregate(data = df_time, . ~ covars, FUN = function(x) quantile(x , probs = c(0.025, 0.975) ))
  perc_time <- Reduce(rbind,lapply(2:ncol(perc_time), function(x) data.frame(low = perc_time[[x]][,1], 
                                                                             high = perc_time[[x]][,2],
                                                                             PC = as.numeric(substr(names(perc_time)[x], 3, nchar(names(perc_time)[x]))), 
                                                                             covars = perc_time$covars)))
  
  object$validation$time$loading <- subset(perc_time, PC %in% PC_time)
  names(object$validation$time$loading)[names(object$validation$time$loading) == 'value'] <- 'loading'
  object$ALASCA$loading$time <- merge(object$ALASCA$loading$time, object$validation$time$loading, all.x = TRUE)
  
  if(object$separateTimeAndGroup){
    df_group <- Reduce(rbind,lapply(objectlist, function(x) x$pca$loading$group))
    PC_group <- getRelevantPCs(object$pca$loading$explained$group)
    perc_group <- aggregate(data = df_group, . ~ covars, FUN = function(x) quantile(x , probs = c(0.025, 0.975) ))
    perc_group <- Reduce(rbind,lapply(2:ncol(perc_group), function(x) data.frame(low = perc_group[[x]][,1], 
                                                                                 high = perc_group[[x]][,2],
                                                                                 PC = as.numeric(substr(names(perc_group)[x], 3, nchar(names(perc_group)[x]))), 
                                                                                 covars = perc_group$covars)))
    object$validation$group$loading <- subset(perc_group, PC %in% PC_group)
    names(object$validation$group$loading)[names(object$validation$group$loading) == 'value'] <- 'loading'
    object$ALASCA$loading$group <- merge(object$ALASCA$loading$group, object$validation$group$loading, all.x = TRUE)
  }
  return(object)
}

#' Extract percentiles for score
#'
#' This function extract percentiles during validation
#'
#' @inheritParams getValidationPercentiles
#' @return An ALASCA object
getValidationPercentilesScore <- function(object, objectlist){
  if(object$separateTimeAndGroup){
    # Separate time and group effects
    
    df_time <- Reduce(rbind,lapply(objectlist, function(x) x$pca$score$time))
    PC_time <- getRelevantPCs(object$pca$score$explained$time)
    perc_time <- aggregate(data = df_time, . ~ time, FUN = function(x) quantile(x , probs = c(0.025, 0.975) ))
    perc_time <- Reduce(rbind,lapply(2:ncol(perc_time), function(x) data.frame(low = perc_time[[x]][,1], 
                                                                               high = perc_time[[x]][,2],
                                                                               PC = as.numeric(substr(names(perc_time)[x], 3, nchar(names(perc_time)[x]))), 
                                                                               time = perc_time$time)))
    #perc_time <- switchSign(reshape2::melt(object$pca$score$time, id.vars = "time"), perc_time)
    object$validation$time$score <- subset(perc_time, PC %in% PC_time)
    names(object$validation$time$score)[names(object$validation$time$score) == 'value'] <- 'score'
    object$ALASCA$score$time <- merge(object$ALASCA$score$time, object$validation$time$score, all.x = TRUE)
    
    df_group <- Reduce(rbind,lapply(objectlist, function(x) x$pca$score$group))
    PC_group <- getRelevantPCs(object$pca$score$explained$group)
    perc_group <- aggregate(data = df_group, . ~ group + time, FUN = function(x) quantile(x , probs = c(0.025, 0.975) ))
    perc_group <- Reduce(rbind,lapply(3:ncol(perc_group), function(x) data.frame(low = perc_group[[x]][,1], 
                                                                                 high = perc_group[[x]][,2],
                                                                                 PC = as.numeric(substr(names(perc_group)[x], 3, nchar(names(perc_group)[x]))), 
                                                                                 time = perc_group$time,
                                                                                 group = perc_group$group)))
    #perc_group <- switchSign(reshape2::melt(object$pca$score$group, id.vars = c("time", "group")), perc_group, PC_group)
    object$validation$group$score <- subset(perc_group, PC %in% PC_group)
    names(object$validation$group$score)[names(object$validation$group$score) == 'value'] <- 'score'
    object$ALASCA$score$group <- merge(object$ALASCA$score$group, object$validation$group$score, all.x = TRUE)
  }else{
    # Pool time and groups effects
    
    df_time <- Reduce(rbind,lapply(objectlist, function(x) x$pca$score$time))
    PC_time <- getRelevantPCs(object$pca$score$explained$time)
    perc_time <- aggregate(data = df_time, . ~ time + group, FUN = function(x) quantile(x , probs = c(0.025, 0.975) ))
    perc_time <- Reduce(rbind,lapply(3:ncol(perc_time), function(x) data.frame(low = perc_time[[x]][,1], 
                                                                               high = perc_time[[x]][,2],
                                                                               PC = as.numeric(substr(names(perc_time)[x], 3, nchar(names(perc_time)[x]))), 
                                                                               time = perc_time$time,
                                                                               group = perc_time$group)))
    object$validation$time$score <- subset(perc_time, PC %in% PC_time)
    names(object$validation$time$score)[names(object$validation$time$score) == 'value'] <- 'score'
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
getRelevantPCs <- function(x){
  PC <- x > 0.05
  PC[1:2] <- TRUE
  return(which(PC))
}

#' Switch sign
#'
#' This function extract percentiles during validation
#'
#' @param value Point estimate
#' @param perc Uncertainty
#' @param PCs Components
#' @return An ALASCA object
switchSign <- function(value, perc, PCs){
  value$variable <- NULL
  perc <- merge(perc, value, all.x = TRUE)
  perc <- subset(perc, PC %in% PCs)
  for(p in unique(perc$PC)){
    if(any(perc$value[perc$PC == p] < perc$low[perc$PC == p]) | any(perc$value[perc$PC == p] > perc$high[perc$PC == p])){
      perc$high[perc$PC == p] <- perc$high[perc$PC == p]*(-1)
      perc$low[perc$PC == p] <- perc$low[perc$PC == p]*(-1)
    }
  }
  return(perc)
}

#' Validate underlying regression models
#'
#' This function extract percentiles during validation
#'
#' @param object An ALASCA object
#' @param variable Variable names for which models to validate (default `NA` for all)
#' @return An ALASCA object
getRegressionPredictions <- function(object){
  if(!object$minimizeObject){
    # This is not a validation run
    cat("Calculating predictions from regression models...\n")
  }
  newdata <- Reduce(rbind,lapply(unique(object$df$variable), function(x){
    xi <- which(x == names(object$regr.model))
    model <- object$regr.model[[xi]]
    if(object$method == "LMM"){
      data.frame(
        pred = lme4:::predict.merMod(model, re.form=NA),
        time = object$df[variable == x,time],
        group = object$df[variable == x,group],
        variable = x
      )
    }else if(object$method == "LM"){
      data.frame(
        pred = predict(model),
        time = object$df[variable == x,time],
        group = object$df[variable == x,group],
        variable = x
      )
    }
  }))
  object$mod.pred <- aggregate(data = newdata, pred~time+group + variable, FUN = "mean")
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
prepareValidationRun <- function(object){
  if(object$validationMethod == "loo"){
    # Use leave-one-out validation
    selectedParts <- data.frame()
    
    if(object$method %in% c("LMM", "Rfast")){
      # For each group, divide the participants into nValFold groups, and select nValFold-1 of them
      selectedParts <- lapply(unique(object$stratificationVector), function(gr){
        selectedParts_temp_all <- unique(object$df[object$stratificationVector == gr,ID])
        selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$nValFold
        selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
        selectedParts_temp_all[selectedParts_temp_ticket != 1]
      })
      
      temp_object <- ALASCA(validationObject = object,
                            validationParticipants = object$df[,ID] %in% unlist(selectedParts))
    }else if(df$method == "LM"){
      object$df$ID <- c(1:nrow(object$df))
      # For each group, divide the participants into nValFold groups, and select nValFold-1 of them
      selectedParts <- lapply(unique(object$stratificationVector), function(gr){
        selectedParts_temp_all <- unique(object$df[object$stratificationVector == gr,ID])
        selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$nValFold
        selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
        selectedParts_temp_all[selectedParts_temp_ticket != 1]
      })
      
      temp_object <- ALASCA(validationObject = object,
                            validationParticipants = object$df[,ID] %in% unlist(selectedParts))
    }
  }else if(object$validationMethod == "permutation2"){
    # Validation by permutation
    # Randomize samples
    parts <- data.frame(
      time = object$df$time,
      group = object$df$group,
      ID = object$df[,ID]
    )
    parts <- parts[!duplicated(parts),]
    times <- parts$time
    groups <- parts$group
    
    # sample id
    parts_orig <- paste0(object$df[,ID], object$df$time, object$df$group)
    u_parts_orig <- unique(parts_orig)
    
    temp_object <- object
    for(i in 1:length(u_parts_orig)){
      rIDTime <- sample(seq_along(times), 1)
      rIDGroup <- sample(seq_along(groups), 1)
      temp_object$dfRaw$time[parts_orig == u_parts_orig[i]] <- times[rIDTime]
      times <- times[-rIDTime]
      temp_object$dfRaw$group[parts_orig == u_parts_orig[i]] <- groups[rIDGroup]
      groups <- groups[-rIDGroup]
    }
    temp_object <- ALASCA(validationObject = temp_object,
                          validationParticipants = !is.na(object$df[,ID]))
  }else if(object$validationMethod == "permutation"){
    # Validation by permutation
    # Randomize individuals across groups
    # Randomize time within individual
    parts_g <- data.frame(
      group = object$df$group,
      ID = object$df[,ID]
    )
    parts_g <- parts_g[!duplicated(parts_g),]
    groups <- parts_g$group
    
    # sample id
    parts_orig <- unique(paste0(object$df[,ID]))
    
    temp_object <- object
    for(i in unique(object$df[,ID])){
      rIDGroup <- sample(seq_along(groups), 1)
      temp_object$dfRaw$group[temp_object$df[,ID] == i] <- groups[rIDGroup]
      groups <- groups[-rIDGroup]
      times <- unique(temp_object$dfRaw$time[temp_object$df[,ID] == i])
      for(j in unique(temp_object$dfRaw$time[temp_object$df[,ID] == i])){
        rIDTime <- sample(seq_along(times), 1)
        temp_object$dfRaw$time[temp_object$df[,ID] == i & temp_object$df$time == j] <- times[rIDTime]
        times <- times[-rIDTime]
      }
    }
    temp_object <- ALASCA(validationObject = temp_object,
                          validationParticipants = !is.na(object$df[,ID]))
  }
  
  return(temp_object)
}

#' Get P values
#'
#' This function calculates P values
#'
#' @param target The full-model ALASCA object
#' @param objectlist List of permutated ALASCA objects
#' @return An ALASCA object
getPermutationPValues <- function(target, objectlist){
  
  if(target$separateTimeAndGroup){
    # Calculate the Frobenius sum squares of the first l principal components
    SS_f_time <- Reduce(rbind,lapply(seq_along(objectlist), function(i){
      tmp <- reshape2::melt(objectlist[[i]]$pca$score$time, id.vars = c("time"))
      tmp$variable <- as.numeric(gsub("PC","",tmp$variable))
      tmp2 <- aggregate(data=tmp[tmp$variable %in% getRelevantPCs(objectlist[[i]]$pca$loading$explained$time),], value~time, FUN = function(x) sum(x^2))
      tmp2$model <- i
      tmp2
    }))
    SS_time <- Reduce(rbind,lapply(1, function(i){
      tmp <- reshape2::melt(target$pca$score$time, id.vars = c("time"))
      tmp$variable <- as.numeric(gsub("PC","",tmp$variable))
      tmp2 <- aggregate(data=tmp[tmp$variable %in% getRelevantPCs(target$pca$loading$explained$time),], value~time, FUN = function(x) sum(x^2))
      tmp2
    }))
    SS_f_group <- Reduce(rbind,lapply(seq_along(objectlist), function(i){
      tmp <- reshape2::melt(objectlist[[i]]$pca$score$group, id.vars = c("time","group"))
      tmp$variable <- as.numeric(gsub("PC","",tmp$variable))
      tmp2 <- aggregate(data=tmp[tmp$variable %in% getRelevantPCs(objectlist[[i]]$pca$loading$explained$group),], value~group+time, FUN = function(x) sum(x^2))
      tmp2$model <- i
      tmp2
    }))
    SS_group <- Reduce(rbind,lapply(1, function(i){
      tmp <- reshape2::melt(target$pca$score$group, id.vars = c("time","group"))
      tmp$variable <- as.numeric(gsub("PC","",tmp$variable))
      tmp2 <- aggregate(data=tmp[tmp$variable %in% getRelevantPCs(target$pca$loading$explained$group),], value~group+time, FUN = function(x) sum(x^2))
      tmp2
    }))
    
    pvals_time <- Reduce(rbind,lapply(unique(SS_time$time), function(x){
      data.frame(
        p.value = sum( SS_f_time$value[SS_f_time$time == x] >= SS_time$value[SS_time$time == x] ) / sum(SS_f_time$time == x),
        effect = x,
        nRuns = target$nValRuns
      )
    }))
    
    pvals_group <- Reduce(rbind,lapply(unique(paste(SS_group$time,SS_group$group)), function(x){
      data.frame(
        p.value = sum( SS_f_group$value[paste(SS_f_group$time,SS_f_group$group) == x] >= SS_group$value[paste(SS_group$time,SS_group$group) == x] ) / sum(paste(SS_f_group$time,SS_f_group$group) == x),
        effect = x,
        nRuns = target$nValRuns
      )
    }))
    
    pvals <- rbind(pvals_time,pvals_group)
    
  }else{
    # Calculate the Frobenius sum squares of the first l principal components
    SS_f <- Reduce(rbind,lapply(seq_along(objectlist), function(i){
      tmp <- reshape2::melt(objectlist[[i]]$pca$score$time, id.vars = c("time","group"))
      tmp$variable <- as.numeric(gsub("PC","",tmp$variable))
      tmp2 <- aggregate(data=tmp[tmp$variable %in% getRelevantPCs(objectlist[[i]]$pca$loading$explained$time),], value~group+time, FUN = function(x) sum(x^2))
      tmp2$model <- i
      tmp2
    }))
    SS <- Reduce(rbind,lapply(1, function(i){
      tmp <- reshape2::melt(target$pca$score$time, id.vars = c("time","group"))
      tmp$variable <- as.numeric(gsub("PC","",tmp$variable))
      tmp2 <- aggregate(data=tmp[tmp$variable %in% getRelevantPCs(target$pca$loading$explained$time),], value~group+time, FUN = function(x) sum(x^2))
      tmp2
    }))
    pvals <- Reduce(rbind,lapply(unique(paste(SS$time,SS$group)), function(x){
      data.frame(
        p.value = sum( SS_f$value[paste(SS_f$time,SS_f$group) == x] >= SS$value[paste(SS$time,SS$group) == x] ) / sum(paste(SS_f$time,SS_f$group) == x),
        effect = x,
        nRuns = target$nValRuns
      )
    }))
  }
  
  min_p_value <- 1/pvals$nRuns[1]
  pvals$p.value[pvals$p.value < min_p_value] <- min_p_value
  target$pvals <- pvals
  return(target)
}