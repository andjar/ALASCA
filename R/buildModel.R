#' Organize the ALASCA model construction
#'
#' This function builds the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
buildModel <- function(object){
  if(!object$minimizeObject){
    # This is not a validation run
    cat("Calculating ",object$method," coefficients...\n")
  }
  object <- runRegression(object)
  if(object$method != "Rfast"){
    # With Rfast, we've already got the coefficients
    object <- getRegressionCoefficients(object)
  }
  
  if(!object$minimizeObject){
    # This is not a validation run
    cat("Finished calculating regression coefficients!\n")
  }
  if(object$method != "Rfast"){
    if(object$minimizeObject){
      if(object$validateRegression){
        object <- getRegressionPredictions(object)
      }
    }else{
      object <- getRegressionPredictions(object)
    }
  }
  object <- removeCovars(object)
  object <- separateLMECoefficients(object)
  object <- getEffectMatrix(object)
  object <- doPCA(object)
  object <- cleanPCA(object)
  object <- cleanALASCA(object)
  return(object)
}

#' Run regressions
#'
#' This function runs the underlying regression models
#'
#' @param object An ALASCA object
#' @return An ALASCA object
runRegression <- function(object){
  
  if(object$method == "Rfast"){
    object$RegressionCoefficients <- Reduce(rbind,lapply(unique(object$df$variable), function(x){
      #start.time <- Sys.time()
      df <- object$df[variable == x]
      modmat <- model.matrix(object$newformula, data = df)
      if(object$forceEqualBaseline){
        modmat <- modmat[,!grepl(paste0("time",levels(object$df$time)[1]), colnames(modmat))]
      }
      mod <- data.frame(
        estimate = Rfast::rint.reg(y = df[,value], 
                                   x = modmat[,2:ncol(modmat)], 
                                   id = as.numeric(factor(df[,ID])), 
                                   ranef = FALSE)$be,
        pvalue = 1,
        covar = as.character(x),
        variable = colnames(modmat)
      )
      #end.time <- Sys.time()
      #cat("\n\n",end.time - start.time,"\n")
      mod
    }))
    return(object)
  }else if(object$method == "LM"){
        object$regr.model <- lapply(unique(object$df$variable), function(x){
          modmat <- model.matrix(object$formula, data = object$df[variable == x])
          modmat <- modmat[,-1]
          if(object$forceEqualBaseline){
            modmat <- modmat[,!grepl(paste0("time",unique(object$df$time)[1]), colnames(modmat))]
          }
          environment(object$newformula) <- environment()
          regr.model <- lm(object$newformula, data = object$df[variable == x])
          attr(regr.model, "name") <- x
          regr.model
        })
  }else{
    object$regr.model <- lapply(unique(object$df$variable), function(x){
      modmat <- model.matrix(object$formula, data = object$df[variable == x])
      modmat <- modmat[,-1]
      if(object$forceEqualBaseline){
        modmat <- modmat[,!grepl(paste0("time",levels(object$df$time)[1]), colnames(modmat), fixed = TRUE)]
      }
      #modmat <- modmat[,ncol(modmat):1]
      environment(object$newformula) <- environment()
      regr.model <- lmerTest::lmer(object$newformula, data = object$df[variable == x])
      attr(regr.model, "name") <- x
      regr.model
    })
  }
  names(object$regr.model) <- unique(object$df$variable)
  return(object)
}

#' Get regression coefficients
#'
#' This function extract the regression coefficients for the ALASCA model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
getRegressionCoefficients <- function(object){
    fdf <- Reduce(rbind,lapply(object$regr.model, function(y){
      if(object$method == "LM"){
        tmp_ef <- coef(y)
        a <- as.data.frame(summary(y)[["coefficients"]][,c(1,4)])
      }else{
        tmp_ef <- lme4::fixef(y)
        a <- as.data.frame(summary(y)[["coefficients"]][,c(1,5)])
      }
      a$covar <- attr(y, "name")
      a$variable <- rownames(a)
      rownames(a) <- NULL
      a
    }))
    fdf$variable <- gsub("modmat","",fdf$variable, fixed = TRUE)
    
    colnames(fdf) <- c("estimate", "pvalue", "covar", "variable")
    
    if(!is.na(object$pAdjustMethod)){
      cat("Adjusting p values...\n")
      object$RegressionCoefficients$pvalue_adj <- NA
      for(i in unique(object$RegressionCoefficients$covar)){
        object$RegressionCoefficients$pvalue_adj[object$RegressionCoefficients$covar == i, ] <- p.adjust(object$RegressionCoefficients$pvalue[object$RegressionCoefficients$covar == i, ], method = object$pAdjustMethod)
      }
    }
    object$RegressionCoefficients <- fdf
  return(object)
}

#' Remove unwanted covariables
#'
#' This function removes coefficients that we do not want in our PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
removeCovars <- function(object){
  object$CovarCoefficients <- data.frame()
  for(i in unique(object$covars)){
    object$CovarCoefficients <- rbind(object$CovarCoefficients, subset(object$RegressionCoefficients, substr(variable, 1, nchar(i)) == i))
    object$RegressionCoefficients <- subset(object$RegressionCoefficients, substr(variable, 1, nchar(i)) != i)
  }
  return(object)
}

#' Separate time and group effects
#'
#' This function separates time and group variables if separateTimeAndGroup = TRUE
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
separateLMECoefficients <- function(object){
  object$RegressionCoefficients$comp <- "TIME"
  if(object$separateTimeAndGroup){
    object$RegressionCoefficients$comp[!(object$RegressionCoefficients$variable == "(Intercept)" |
                                  (substr(object$RegressionCoefficients$variable, 1, 4) == "time" & !grepl(":",object$RegressionCoefficients$variable, fixed = "TRUE")))
                                 ] <- "GROUP"
  }
  return(object)
}

#' Get effect matrix
#'
#' This function separates time and group variables if separateTimeAndGroup = TRUE
#'
#' @param object An ALASCA object
#' @return An ALASCA object
getEffectMatrix <- function(object){
  if(!object$minimizeObject){
    cat("Calculating effect matrix\n")
  }
  parts <- object$df[variable == object$df$variable[1]]
  #parts <- object$df[!duplicated(cbind(object$df$ID, object$df$time))]
  Dmatrix <- model.matrix(object$formula, data = object$df[variable == object$df$variable[1]])
  #Dmatrix <- Dmatrix[,ncol(Dmatrix):1]

  if(object$separateTimeAndGroup){
    BmatrixTime <- object$RegressionCoefficients[object$RegressionCoefficients$comp == "TIME",c("covar","estimate","variable")]
    BmatrixTime <- reshape2::dcast(BmatrixTime, formula = variable~covar, value.var = "estimate")
    selectDColumnsTime <- colnames(Dmatrix) %in% BmatrixTime$variable
    rowOrder <- c()
    for(i in 1:nrow(BmatrixTime)){
      rowOrder[i] <- which(BmatrixTime$variable == colnames(Dmatrix[,selectDColumnsTime])[i])
    }
    BmatrixTime <- BmatrixTime[rowOrder,]
    if(any(colnames(Dmatrix[,selectDColumnsTime]) != BmatrixTime$variable)){
      stop("Column mismatch for time in getEffectMatrix")
    }
    AmatrixTime <- as.data.frame(as.matrix(Dmatrix[,selectDColumnsTime])%*%as.matrix(BmatrixTime[,2:ncol(BmatrixTime)]))
    AmatrixTime$comp <- "TIME"
    
    BmatrixGroup <- object$RegressionCoefficients[object$RegressionCoefficients$comp == "GROUP",c("covar","estimate","variable")]
    BmatrixGroup <- reshape2::dcast(BmatrixGroup, formula = variable~covar, value.var = "estimate")
    selectDColumnsGroup <- colnames(Dmatrix) %in% BmatrixGroup$variable
    rowOrder <- c()
    for(i in 1:nrow(BmatrixGroup)){
      rowOrder[i] <- which(BmatrixGroup$variable == colnames(Dmatrix[,selectDColumnsGroup])[i])
    }
    BmatrixGroup <- BmatrixGroup[rowOrder,]
    if(any(colnames(Dmatrix[,selectDColumnsGroup]) != BmatrixGroup$variable)){
      stop("Column mismatch for group in getEffectMatrix")
    }
    AmatrixGroup <- as.data.frame(as.matrix(Dmatrix[,selectDColumnsGroup])%*%as.matrix(BmatrixGroup[,2:ncol(BmatrixGroup)]))
    AmatrixGroup$comp <- "GROUP"
    object$effect.matrix <- rbind(AmatrixTime, AmatrixGroup)
    
  }else{
    BmatrixTime <- object$RegressionCoefficients[object$RegressionCoefficients$comp == "TIME",c("covar","estimate","variable")]
    BmatrixTime <- reshape2::dcast(BmatrixTime, formula = variable~covar, value.var = "estimate")
    selectDColumnsTime <- colnames(Dmatrix) %in% BmatrixTime$variable
    rowOrder <- c()
    for(i in 1:nrow(BmatrixTime)){
      rowOrder[i] <- which(BmatrixTime$variable == colnames(Dmatrix[,selectDColumnsTime])[i])
    }
    BmatrixTime <- BmatrixTime[rowOrder,]
    if(any(colnames(Dmatrix[,selectDColumnsTime]) != BmatrixTime$variable)){
      stop("Column mismatch for time in getEffectMatrix")
    }
    AmatrixTime <- as.data.frame(as.matrix(Dmatrix[,selectDColumnsTime])%*%as.matrix(BmatrixTime[,2:ncol(BmatrixTime)]))
    AmatrixTime$comp <- "TIME"
    object$effect.matrix <- AmatrixTime
  }
  
  object$parts$time <- parts$time
  if(object$keepTerms != ""){
    keepTerms <- c("group", object$keepTerms)
    object$parts$group <- apply( parts[ , ..keepTerms ] , 1 , paste , collapse = " - " )
  }else{
    object$parts$group <- parts$group
  }
  if(!object$minimizeObject){
    cat("Finished calculating effect matrix!\n")
  }
  return(object)
}