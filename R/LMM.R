#' Get LMM coefficients
#'
#' This function calculates the LMM coefficients for the ALASCA-model
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
getLMECoefficients <- function(object){
  if(!object$minimizeObject){
    cat("Calculating LMM coefficients...\n")
  }
  lmer.model <- list()
  fdf <- data.frame()
  cyts <- unique(object$df$variable)
  cc <- 1
  ccc <- 1
  #start.time <- Sys.time()
  object$lmer.models <- lapply(unique(object$df$variable), function(i){
      lmer.model <- lmerTest::lmer(object$formula, data = subset(object$df, variable == i))
      if(object$forceEqualBaseline){
        X <- model.matrix(lmer.model)
        baselineLabel <- paste0("time", unique(object$df$time)[1])
        X <- X[, substr(colnames(X),1,nchar(baselineLabel)) != baselineLabel]
        newFormula <- as.character(object$formula)
        newFormulaPred <- strsplit(as.character(newFormula[3]), "\\+")[[1]]
        newFormulaPred <- newFormulaPred[Reduce(cbind,lapply(newFormulaPred, function (x) grepl("\\(",x)))]
        newFormula <- paste(newFormula[2],"~","X +",newFormulaPred, collapse = " ")
        object$newFormula <- formula(newFormula)
        lmer.model <- lmerTest::lmer(object$newFormula, data = subset(object$df, variable == i))
      }
      attr(lmer.model, "name") <- i
      lmer.model
    }
  )
  #end.time <- Sys.time()
  #cat("Time 1: ", end.time - start.time, "\n")
  names(object$lmer.models) <- unique(object$df$variable)
  #start.time <- Sys.time()
  fdf <- Reduce(rbind,lapply(object$lmer.models, function(y){
    tmp_ef <- lme4::fixef(y)
    a <- as.data.frame(summary(y)[["coefficients"]][,c(1,5)])
    a$covar <- attr(y, "name")
    a$variable <- rownames(a)
    rownames(a) <- NULL
    a
  }))
  if(object$forceEqualBaseline){
    fdf$variable[substr(fdf$variable,1,1) == "X"] <- substr(fdf$variable[substr(fdf$variable,1,1) == "X"],2,nchar(fdf$variable[substr(fdf$variable,1,1) == "X"]))
  }
  # fdf <- Reduce(rbind,lapply(seq_along(object$lmer.models), function(y, n, i){
  #       tmp_ef <- lme4::fixef(y[[i]])
  #       a <- as.data.frame(summary(y[[i]])[["coefficients"]][,c(1,5)])
  #       a$covar <- n[[i]]
  #       a$variable <- rownames(a)
  #       rownames(a) <- NULL
  #       a
  #     }, y = object$lmer.models, n = names(object$lmer.models)
  #   ))
  #end.time <- Sys.time()
  #cat("Time 2: ", end.time - start.time, "\n")
  if(!object$minimizeObject){
    cat("Finished calculating LMM coefficients!\n")
  }
  
  colnames(fdf) <- c("estimate", "pvalue", "covar", "variable")
  
  if(!is.na(object$pAdjustMethod)){
    cat("Adjusting p values...\n")
    object$LMM.coefficients$pvalue_adj <- NA
    for(i in unique(object$LMM.coefficients$covar)){
      object$LMM.coefficients$pvalue_adj[object$LMM.coefficients$covar == i, ] <- p.adjust(object$LMM.coefficients$pvalue[object$LMM.coefficients$covar == i, ], method = object$pAdjustMethod)
    }
  }
  object$LMM.coefficients <- fdf

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
    object$LMM.coefficients <- subset(object$LMM.coefficients, substr(variable, 1, nchar(i)) != i)
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
  object$LMM.coefficients$comp <- "TIME"
  if(object$separateTimeAndGroup){
    object$LMM.coefficients$comp[!(object$LMM.coefficients$variable == "(Intercept)" |
                                  (substr(object$LMM.coefficients$variable, 1, 4) == "time" & !grepl(":",object$LMM.coefficients$variable)))
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
  Dmatrix <- lme4::getME(object$lmer.model[[1]],"X")
  if(object$forceEqualBaseline){
    colnames(Dmatrix)[substr(colnames(Dmatrix),1,1) == "X"] <- substr(colnames(Dmatrix)[substr(colnames(Dmatrix),1,1) == "X"],2,nchar(colnames(Dmatrix)[substr(colnames(Dmatrix),1,1) == "X"]))
  }
  if(object$separateTimeAndGroup){
    BmatrixTime <- object$LMM.coefficients[object$LMM.coefficients$comp == "TIME",c("covar","estimate","variable")]
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
    
    BmatrixGroup <- object$LMM.coefficients[object$LMM.coefficients$comp == "GROUP",c("covar","estimate","variable")]
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
    BmatrixTime <- object$LMM.coefficients[object$LMM.coefficients$comp == "TIME",c("covar","estimate","variable")]
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