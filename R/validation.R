

#' Validate the RMASCA model LMM
#'
#' This function calculates the LMM coefficients for the RMASCA-model
#'
#' @param object An RMASCA object
#' @return An RMASCA object
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

  for(ii in 1:object$nValRuns){
    cat("- Run ",ii," of ",object$nValRuns,"\n")
    selectedParts <- c()
    for(gr in unique(object$df$group)){
      selectedParts_temp_all <- unique(object$df[object$df$group == gr,partColumn])
      selectedParts_temp_ticket <- sample(1:object$nValFold, length(selectedParts_temp_all), replace = TRUE)
      selectedParts_temp <- selectedParts_temp_all[selectedParts_temp_ticket != 1]
      selectedParts <- c(selectedParts, selectedParts_temp)
    }
    temp_object[[ii]] <- RMASCA(df = NA,
                                validate = FALSE, # to avoid recursion
                                minimizeObject = TRUE,
                                validationObject = object,
                                validationParticipants = object$df[,partColumn] %in% selectedParts)
    temp_object[[ii]] <- rotateMatrix(object = temp_object[[ii]], target = object)
  }
  object <- getValidationPercentiles(object, objectlist = temp_object)

  return(object)
}

#' Rotate PCA
#'
#' This function calculates the LMM coefficients for the RMASCA-model
#'
#' @param object RMASCA object to be rotated (and returned)
#' @param target RMASCA object acting as target
#' @return An RMASCA object
rotateMatrix <- function(object, target){
  a <- object$RMASCA$loading
  b <- target$RMASCA$loading
  c <- GPArotation::targetT(L = as.matrix(a$time[,1:2]), Target = as.matrix(b$time[,1:2]))
  object$RMASCA$loading$time[,1:2] <- c$loadings
  object$RMASCA$score$time[,1:2] <- as.matrix(object$RMASCA$score$time[,1:2]) %*% solve(t(c$Th))
  if(object$separateTimeAndGroup){
    c <- GPArotation::targetT(L = as.matrix(a$group[,1:2]), Target = as.matrix(b$group[,1:2]))
    object$RMASCA$loading$group[,1:2] <- c$loadings
    object$RMASCA$score$group[,1:2] <- as.matrix(object$RMASCA$score$group[,1:2]) %*% solve(t(c$Th))
  }
  return(object)
}

#' Extract percentiles
#'
#' This function calculates the LMM coefficients for the RMASCA-model
#'
#' @param object RMASCA object
#' @param objectlist List of RMASCA objects
#' @return An RMASCA object
getValidationPercentiles <- function(object, objectlist){
  df_loading_time <- data.frame()
  df_score_time <- data.frame()
  df_loading_group <- data.frame()
  df_score_group <- data.frame()
  if(object$separateTimeAndGroup){
    for(i in 1:length(objectlist)){
      df_loading_time_temp <- getLoadings(objectlist[[i]])$time
      df_loading_time <- rbind(df_loading_time, df_loading_time_temp)

      df_score_time_temp <- getScores(objectlist[[i]])$time
      df_score_time <- rbind(df_score_time, df_score_time_temp)

      df_loading_group_temp <- getLoadings(objectlist[[i]])$group
      df_loading_group <- rbind(df_loading_group, df_loading_group_temp)

      df_score_group_temp <- getScores(objectlist[[i]])$group
      df_score_group <- rbind(df_score_group, df_score_group_temp)
    }
    perc_loading_time <- data.frame()
    perc_score_time <- data.frame()
    perc_loading_group <- data.frame()
    perc_score_group <- data.frame()
    for(p in 1:2){
      for(t in unique(df_loading_time$covars)){
        perc_loading_time_temp <- data.frame(
          a = quantile(df_loading_time[df_loading_time$covars == t,p], probs = c(0.025)),
          b = quantile(df_loading_time[df_loading_time$covars == t,p], probs = c(0.975)),
          c = t,
          d = p
        )
        perc_loading_time <- rbind(perc_loading_time, perc_loading_time_temp)
      }
      for(t in unique(df_score_time$time)){
        perc_score_time_temp <- data.frame(
          a = quantile(df_score_time[df_score_time$time == t,p], probs = c(0.025)),
          b = quantile(df_score_time[df_score_time$time == t,p], probs = c(0.975)),
          c = t,
          d = p
        )
        perc_score_time <- rbind(perc_score_time, perc_score_time_temp)
      }
      for(t in unique(df_loading_group$covars)){
        perc_loading_group_temp <- data.frame(
          a = quantile(df_loading_group[df_loading_group$covars == t,p], probs = c(0.025)),
          b = quantile(df_loading_group[df_loading_group$covars == t,p], probs = c(0.975)),
          c = t,
          d = p
        )
        perc_loading_group <- rbind(perc_loading_group, perc_loading_group_temp)
      }
      for(t in unique(df_score_group$time)){
        for(g in unique(df_score_group$group)){
          perc_score_group_temp <- data.frame(
            a = quantile(df_score_group[df_score_group$time == t & df_score_group$group == g,p], probs = c(0.025)),
            b = quantile(df_score_group[df_score_group$time == t & df_score_group$group == g,p], probs = c(0.975)),
            c = t,
            g = g,
            d = p
          )
          perc_score_group <- rbind(perc_score_group, perc_score_group_temp)
        }
      }
    }
    colnames(perc_loading_time) <- c("low", "high", "covars","PC")
    rownames(perc_loading_time) <- NULL
    colnames(perc_score_time) <- c("low", "high", "time","PC")
    rownames(perc_score_time) <- NULL
    colnames(perc_loading_group) <- c("low", "high", "covars","PC")
    rownames(perc_loading_group) <- NULL
    colnames(perc_score_group) <- c("low", "high", "time","group","PC")
    rownames(perc_score_group) <- NULL

    # Check if we need to switch sign
    loadings <- getLoadings(object)$time
    for(comp in 1:2){
      temp_limits <- perc_loading_time[perc_loading_time$PC == comp,]
      temp_limits <- merge(temp_limits, loadings[,c(comp, ncol(loadings))], by = "covars")
      if( temp_limits$low > temp_limits[,3] | temp_limits[,3] > temp_limits$high ){
        perc_loading_time[perc_loading_time$PC == comp,] <- perc_loading_time[perc_loading_time$PC == comp,]*-1
        perc_score_time[perc_score_time$PC == comp,] <- perc_score_time[perc_score_time$PC == comp,]*-1
      }
    }

    # Check if we need to switch sign
    loadings <- getLoadings(object)$group
    for(comp in 1:2){
      temp_limits <- perc_loading_group[perc_loading_group$PC == comp,]
      temp_limits <- merge(temp_limits, loadings[,c(comp, ncol(loadings))], by = "covars")
      if( temp_limits$low > temp_limits[,3] | temp_limits[,3] > temp_limits$high ){
        perc_loading_group[perc_loading_group$PC == comp,] <- perc_loading_group[perc_loading_group$PC == comp,]*-1
        perc_score_group[perc_score_group$PC == comp,] <- perc_score_group[perc_score_group$PC == comp,]*-1
      }
    }

    object$validation$time$loading  <- perc_loading_time
    object$validation$time$score  <- perc_score_time
    object$validation$group$loading  <- perc_loading_group
    object$validation$group$score  <- perc_score_group

  }else{
    for(i in 1:length(objectlist)){
      df_loading_time_temp <- objectlist[[i]]$RMASCA$loading$time
      df_loading_time <- rbind(df_loading_time, df_loading_time_temp)

      df_score_time_temp <- objectlist[[i]]$RMASCA$score$time
      df_score_time <- rbind(df_score_time, df_score_time_temp)
    }
    perc_loading_time <- data.frame()
    perc_score_time <- data.frame()
    for(p in 1:2){
      for(t in unique(df_loading_time$covars)){
        perc_loading_time_temp <- data.frame(
          a = quantile(df_loading_time[df_loading_time$covars == t,p], probs = c(0.025)),
          b = quantile(df_loading_time[df_loading_time$covars == t,p], probs = c(0.975)),
          c = t,
          d = p
        )
        perc_loading_time <- rbind(perc_loading_time, perc_loading_time_temp)
      }
      for(t in unique(df_score_time$time)){
        for(g in unique(df_score_time$group)){
          perc_score_time_temp <- data.frame(
            a = quantile(df_score_time[df_score_time$time == t & df_score_time$group == g,p], probs = c(0.025)),
            b = quantile(df_score_time[df_score_time$time == t & df_score_time$group == g,p], probs = c(0.975)),
            c = t,
            g = g,
            d = p
          )
          perc_score_time <- rbind(perc_score_time, perc_score_time_temp)
        }
      }
    }
    colnames(perc_loading_time) <- c("low", "high", "covars","PC")
    rownames(perc_loading_time) <- NULL
    colnames(perc_score_time) <- c("low", "high", "time","group","PC")
    rownames(perc_score_time) <- NULL

    # Check if we need to switch sign
    loadings <- getLoadings(object)$time
    for(comp in 1:2){
      temp_limits <- perc_loading_time[perc_loading_time$PC == comp,]
      temp_limits <- merge(temp_limits, loadings[,c(comp, ncol(loadings))], by = "covars")
      if( temp_limits$low > temp_limits[,3] | temp_limits[,3] > temp_limits$high ){
        perc_loading_time[perc_loading_time$PC == comp,] <- perc_loading_time[perc_loading_time$PC == comp,]*-1
        perc_score_time[perc_score_time$PC == comp,] <- perc_score_time[perc_score_time$PC == comp,]*-1
      }
    }

    object$validation$time$loading  <- perc_loading_time
    object$validation$time$score  <- perc_score_time
  }


  return(object)
}
