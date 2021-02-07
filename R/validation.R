#' Validate the ALASCA model LMM
#'
#' This function performs leave-one-out robustness testing of your ALASCA model. If you didn't specify the number of runs `nValRuns` when initializing the model (see \code{\link{ALASCA}}), you can do it by running for example `model$nValRuns <- 100` prior to calling `validate`. Your dataset is divided into `nValFold` partitions, keeping group proportions, and one of these are left out. `nValFold` is set the same way as  `nValRuns`.
#'
#' @param object An ALASCA object
#' @param participantColumn The name of the column containing participant identifier. Needed if not set during initialization of the model.
#' @return An ALASCA object
#' 
#' @examples
#' load("PE.Rdata")
#' model$nValRuns <- 10
#' model.val <- validate(model, participantColumn = "ID")
#' 
#' @export
validate <- function(object, participantColumn = FALSE){
  object$validate <- TRUE

  if(participantColumn == FALSE){
    if(object$participantColumn == FALSE){
      stop("You need to specify the column containing participant id in `participantColumn`")
    }
  }else{
    object$participantColumn <- participantColumn
  }
  cat("Running validation...\n")
  partColumn <- which(colnames(object$df) == object$participantColumn)
  temp_object <- list()
  time_mean <- c()

  for(ii in 1:object$nValRuns){
    cat("- Run ",ii," of ",object$nValRuns,"\n")
    start.time <- Sys.time()
    
    # For each group, divide the participants into nValFold groups, and select nValFold-1 of them
    selectedParts <- Reduce(c,lapply(unique(object$df$group), function(gr){
      selectedParts_temp_all <- unique(object$df[object$df$group == gr,partColumn])
      selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$nValFold
      selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
      selectedParts_temp <- selectedParts_temp_all[selectedParts_temp_ticket != 1]
      selectedParts_temp
    }))
    
    # Make ALASCA model from the selected subset
    temp_object[[ii]] <- ALASCA(validationObject = object,
                                validationParticipants = object$df[,partColumn] %in% selectedParts)
    
    # Rotate new loadings/scores to the original model
    temp_object[[ii]] <- rotateMatrix(object = temp_object[[ii]], target = object)
    
    end.time <- Sys.time()
    time_mean[ii] <- end.time - start.time
    cat("--- Used ",round(time_mean[ii],2)," seconds. Est. time remaining: ",round((object$nValRuns-ii)*mean(time_mean),2)," seconds \n")
  }
  
  # Tried parallelization, but it didn't improve performance significantly
  # doMC::registerDoMC(cores=parallel::detectCores()-1)
  # temp_object <- foreach::"%dopar%"(foreach::foreach(1:object$nValRuns, .inorder = FALSE, .packages=c("ALASCA","foreach")),{
  #   selectedParts <- c()
  #   for(gr in unique(object$df$group)){
  #     selectedParts_temp_all <- unique(object$df[object$df$group == gr,partColumn])
  #     selectedParts_temp_ticket <- sample(1:object$nValFold, length(selectedParts_temp_all), replace = TRUE)
  #     selectedParts_temp <- selectedParts_temp_all[selectedParts_temp_ticket != 1]
  #     selectedParts <- c(selectedParts, selectedParts_temp)
  #   }
  #   temp_object <- ALASCA(df = NA,
  #                               validate = FALSE, # to avoid recursion
  #                               minimizeObject = TRUE,
  #                               validationObject = object,
  #                               validationParticipants = object$df[,partColumn] %in% selectedParts)
  #   temp_object <- rotateMatrix(object = temp_object, target = object)
  # })
  object <- getValidationPercentiles(object, objectlist = temp_object)
  object$validation$temp_objects <- temp_object
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
  
  
  a <- object$pca$loading
  b <- target$pca$loading
  
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
    c <- procrustes(loadings= as.matrix(a$time[,PCloading]* signMatrix[i,]),
                    target = as.matrix(b$time[,PCloading]))
    sum((b$time[,PCloading] - c$procrust)^2)
  }))
  minSignVar <- which(signVar == min(signVar))[1]
  object$pca$loading$time[,PCloading] <- object$pca$loading$time[,PCloading] * signMatrix[minSignVar,]
  object$pca$score$time[,PCloading] <- object$pca$score$time[,PCloading] * signMatrix[minSignVar,]
  a <- object$pca$loading
  
  c <- procrustes(loadings= as.matrix(a$time[,PCloading]),
                  target = as.matrix(b$time[,PCloading]))
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
      c <- procrustes(loadings= as.matrix(a$group[,PCloading]* signMatrix[i,]),
                      target = as.matrix(b$group[,PCloading]))
      sum((b$group[,PCloading] - c$procrust)^2)
    }))
    minSignVar <- which(signVar == min(signVar))[1]
    object$pca$loading$group[,PCloading] <- object$pca$loading$group[,PCloading] * signMatrix[minSignVar,]
    object$pca$score$group[,PCloading] <- object$pca$score$group[,PCloading] * signMatrix[minSignVar,]
    a <- object$pca$loading
    
    c <- procrustes(loadings= as.matrix(a$group[,PCloading]),
                    target = as.matrix(b$group[,PCloading]))
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
    perc_time <- switchSign(reshape2::melt(object$pca$score$time, id.vars = "time"), perc_time)
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
    perc_group <- switchSign(reshape2::melt(object$pca$score$group, id.vars = c("time", "group")), perc_group)
    object$validation$group$score <- subset(perc_group, PC %in% PC_group)
    names(object$validation$group$score)[names(object$validation$group$score) == 'value'] <- 'score'
    object$ALASCA$score$group <- merge(object$ALASCA$score$group, object$validation$group$score, all.x = TRUE)
  }else{
    # Pool time and groups effects
    
    df_time <- Reduce(rbind,lapply(objectlist, function(x) x$pca$score$time))
    PC_time <- getRelevantPCs(object$pca$score$explained$time > 0.05)
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
#' @param value
#' @param perc
#' @param PCs
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
