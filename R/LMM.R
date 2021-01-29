#' Get LMM coefficients
#'
#' This function calculates the LMM coefficients for the RMASCA-model
#'
#' @param object An RMASCA object to be sanitized
#' @return An RMASCA object
getLMECoefficients <- function(object){
  if(!object$minimizeObject){
    cat("Calculating LMM coefficients...\n")
  }
  lmer.model <- list()
  cc <- 1
  for(i in unique(df$variable)){
    if(!object$minimizeObject){
      cat("- ",i,"\n")
    }
    lmer.model[[cc]] <- lmerTest::lmer(object$formula, data = subset(object$df, variable == i))
    cc <- cc + 1
  }
  if(!object$minimizeObject){
    object$lmer.models <- lmer.model
    cat("Finished calculating LMM coefficients!\n")
  }

  fdf <- data.frame()
  cyts <- unique(object$df$variable)
  cc <- 1
  for(i in 1:length(lmer.model)){
    tmp_ef <- lme4::fixef(lmer.model[[i]])
    for(j in 1:length(tmp_ef)){
      fdf[cc,1] <- cyts[i]
      fdf[cc,2] <- summary(lmer.model[[i]])[["coefficients"]][j,1]
      fdf[cc,3] <- summary(lmer.model[[i]])[["coefficients"]][j,5]
      fdf[cc,4] <- names(tmp_ef)[j]
      cc <- cc + 1
    }
  }

  colnames(fdf) <- c("covar", "estimate", "pvalue", "variable")
  if(!is.na(object$pAdjustMethod)){
    object$LMM.coefficients$pvalue_adj <- p.adjust(object$LMM.coefficients$pvalue, method = object$pAdjustMethod)
  }
  object$LMM.coefficients <- fdf
  return(object)
}

#' Remove unwanted covariables
#'
#' This function removes coefficients that we do not want in our PCA
#'
#' @param object An RMASCA object to be sanitized
#' @return An RMASCA object
removeCovars <- function(object){
  for(i in unique(object$covars)){
    object$LMM.coefficients <- subset(object$LMM.coefficients, substr(variable, 1, nchar(i)) != i)
  }
  return(object)
}

#' Separate time and group effects
#'
#' This function separates time and group variables if separateTimeAndGroup = TRUE
#'
#' @param object An RMASCA object to be sanitized
#' @return An RMASCA object
separateLMECoefficients <- function(object){
  object$LMM.coefficients$comp <- "TIME"
  if(object$separateTimeAndGroup){
    object$LMM.coefficients$comp[!(object$LMM.coefficients$variable == "(Intercept)" |
                                  (substr(object$LMM.coefficients$variable, 1, 4) == "time" & !grepl(":",object$LMM.coefficients$variable)))
                                 ] <- "GROUP"
  }
  return(object)
}

getEffectMatrix <- function(object){
  if(!object$minimizeObject){
    cat("Calculating effect matrix (may take some time...)\n")
  }
  parts <- subset(object$df, variable == object$df$variable[1]) # just to get time and group for each participant

  covars <- unique(object$LMM.coefficients$covar)
  refTime <- levels(object$df$time)[1]
  refGroup <- levels(object$df$group)[1]
  fdf_pred <- data.frame()

  # if time and group effects are separated, we need to do that here:
  if(object$separateTimeAndGroup){
    fdf_time_pred <- data.frame()
    fdf_group_pred <- data.frame()

    for(j in 1:length(covars)){
      for(i in 1:nrow(parts)){
        fdf_time_pred[i,j] <- object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == "(Intercept)"]
        if(parts$time[i] != refTime){
          fdf_time_pred[i,j] <- fdf_time_pred[i,j] +object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == paste0("time",parts$time[i])]
        }
        fdf_group_pred[i,j] <- 0
        if(parts$group[i] != refGroup){
          if(object$hasGroupTerm){
            fdf_group_pred[i,j] <- object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == paste0("group",parts$group[i])]
          }
          if(parts$time[i] != refTime & object$hasInteractionTerm){
            fdf_group_pred[i,j] <- fdf_group_pred[i,j] + object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == paste0("time",parts$time[i],":group",parts$group[i])]
          }
        }
      }
    }
    fdf_time_pred$comp <- "TIME"
    fdf_group_pred$comp <- "GROUP"
    fdf_pred <- rbind(fdf_time_pred, fdf_group_pred)
  }else{
    for(i in 1:nrow(parts)){
      for(j in 1:length(covars)){

        fdf_pred[i,j] <- object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == "(Intercept)"]
        if(parts$time[i] != refTime){
          fdf_pred[i,j] <- fdf_pred[i,j] + object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == paste0("time",parts$time[i])]
        }
        if(parts$group[i] != refGroup){
          if(object$hasGroupTerm){
            fdf_pred[i,j] <- fdf_pred[i,j] + object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == paste0("group",parts$group[i])]
          }
          if(parts$time[i] != refTime & object$hasInteractionTerm){
            fdf_pred[i,j] <- fdf_pred[i,j] + object$LMM.coefficients$estimate[object$LMM.coefficients$covar == covars[j] & object$LMM.coefficients$variable == paste0("time",parts$time[i],":group",parts$group[i])]
          }
        }

      }
    }
    fdf_pred$comp <- "TIME"
  }
  colnames(fdf_pred) <- c(as.character(covars),"comp")
  object$effect.matrix <- fdf_pred
  object$parts$time <- parts$time
  object$parts$group <- parts$group
  if(object$minimizeObject){
    cat("Finished calculating effect matrix!\n")
  }
  return(object)
}
