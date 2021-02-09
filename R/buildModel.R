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
  object <- getRegressionCoefficients(object)
  if(!object$minimizeObject){
    # This is a validation run
    cat("Finished calculating regression coefficients!\n")
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
  
  regr.model <- list()
  fdf <- data.frame()
  cyts <- unique(object$df$variable)
  cc <- 1
  ccc <- 1
  
  if(object$forceEqualBaseline){
    if(object$doDebug){
      cat(".... Starts removing interaction between groups and first time point\n")
    }
    # Remove interaction between group and first time point
    if(object$method == "LM"){
      regr.model <- lm(object$formula, data = subset(object$df, variable == unique(object$df$variable)[1]))
    }else{
      regr.model <- lmerTest::lmer(object$formula, data = subset(object$df, variable == unique(object$df$variable)[1]))
    }
    X <- model.matrix(regr.model)
    baselineLabel <- paste0("time", unique(object$df$time)[1])
    if(object$doDebug){
      cat(".... Baseline found to be: ",baselineLabel,"\n")
      cat(".... Column names of X: ",colnames(X),"\n")
    }
    X <- X[, substr(colnames(X),1,nchar(baselineLabel)) != baselineLabel]
    newFormula <- as.character(object$formula)
    newFormulaPred <- strsplit(as.character(newFormula[3]), "\\+")[[1]]
    newFormulaPred <- newFormulaPred[Reduce(cbind,lapply(newFormulaPred, function (x) grepl("\\(",x)))]
    newFormula <- paste(newFormula[2],"~","X +",newFormulaPred, collapse = " ")
    object$newFormula <- formula(newFormula)
    if(object$doDebug){
      cat(".... Old formula (given as ~ `outcome` `predictors`): ",as.character(object$formula),"\n")
      cat(".... New formula (given as ~ `outcome` `predictors`): ",as.character(object$newFormula),"\n")
      cat(".... Column names of X: ",colnames(X),"\n")
      cat(".... Size of X (rowsxcols): ",nrow(X),"x",ncol(X),"\n")
    }
    object$X <- X
  }
  
  if(object$doDebug){
    cat(".. Beginning regression\n")
  }
  object$regr.model <- lapply(unique(object$df$variable), function(i){
    if(object$doDebug){
      cat(".... Variable: ",i," (",levels(object$df$variable)[i],")\n")
      cat(".... Number of rows in subset: ",nrow(subset(object$df, variable == i)),"\n")
    }
    if(object$forceEqualBaseline){
      if(object$doDebug){
        cat("!!! --- Making special exception --- !!!\n")
        p1 <- paste0(subset(object$df, variable == unique(object$df$variable)[1])$studyno, subset(object$df, variable == unique(object$df$variable)[1])$time)
        p2 <- paste0(subset(object$df, variable == i)$studyno, subset(object$df, variable == i)$time)
        cat(".... Length p1: ",length(p1),"\n")
        cat(".... Length p2: ",length(p2),"\n")
        X <- subset(object$X, p1 %in% p2)
      }
      if(object$method == "LM"){
        regr.model <- lm(object$newFormula, data = subset(object$df, variable == i))
      }else{
        regr.model <- lmerTest::lmer(object$newFormula, data = subset(object$df, variable == i))
      }
    }else{
      if(object$method == "LM"){
        regr.model <- lm(object$formula, data = subset(object$df, variable == i))
      }else{
        regr.model <- lmerTest::lmer(object$formula, data = subset(object$df, variable == i))
      }
    }
    attr(regr.model, "name") <- i
    regr.model
  }
  )
  
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
  if(object$forceEqualBaseline){
    fdf$variable[substr(fdf$variable,1,1) == "X"] <- substr(fdf$variable[substr(fdf$variable,1,1) == "X"],2,nchar(fdf$variable[substr(fdf$variable,1,1) == "X"]))
  }
  
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
  for(i in unique(object$covars)){
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
                                  (substr(object$RegressionCoefficients$variable, 1, 4) == "time" & !grepl(":",object$RegressionCoefficients$variable)))
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
  parts <- subset(object$df, variable == object$df$variable[1])
  if(object$method == "LM"){
    Dmatrix <- model.matrix(object$regr.model[[1]],"X")
  }else{
    Dmatrix <- lme4::getME(object$regr.model[[1]],"X")
  }
  
  if(object$forceEqualBaseline){
    colnames(Dmatrix)[substr(colnames(Dmatrix),1,1) == "X"] <- substr(colnames(Dmatrix)[substr(colnames(Dmatrix),1,1) == "X"],2,nchar(colnames(Dmatrix)[substr(colnames(Dmatrix),1,1) == "X"]))
  }
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
  object$parts$group <- parts$group
  if(!object$minimizeObject){
    cat("Finished calculating effect matrix!\n")
  }
  return(object)
}