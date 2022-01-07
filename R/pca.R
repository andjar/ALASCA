#' Perform PCA
#'
#' This function performs PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
doPCA <- function(object){
  object$pca$time <- prcomp(object$effect.matrix[object$effect.matrix$comp == "TIME",1:(ncol(object$effect.matrix)-1)], scale = FALSE, center = TRUE)
  if(object$separateTimeAndGroup){
    object$pca$group <- prcomp(object$effect.matrix[object$effect.matrix$comp == "GROUP",1:(ncol(object$effect.matrix)-1)], scale = FALSE, center = TRUE)
  }
  return(object)
}

#' Clean the PCA data
#'
#' This function makes the pca output more useful
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
cleanPCA <- function(object){

  loading_time <- as.data.frame(object$pca$time$rotation[,])
  loading_time$covars <- rownames(loading_time)
  object$pca$loading$time <- loading_time

  PC_time <- as.data.frame(object$pca$time$x[,])
  PC_time$time <- object$parts$time
  if(!object$separateTimeAndGroup){
    PC_time$group <- object$parts$group
  }
  object$pca$score$time <- PC_time
  object$pca$score$explained$time <- object$pca$time$sdev^2 / sum(object$pca$time$sdev^2)
  object$pca$loading$explained$time <- object$pca$score$explained$time

  if(object$separateTimeAndGroup){
    loading_group <- as.data.frame(object$pca$group$rotation[,])
    loading_group$covars <- rownames(loading_group)
    object$pca$loading$group <- loading_group
    PC_group <- as.data.frame(object$pca$group$x[,])
    PC_group$time <- rownames(PC_group)
    PC_group$group <- object$parts$group
    PC_group$time <- object$parts$time
    object$pca$score$group <- PC_group
    object$pca$score$explained$group <- object$pca$group$sdev^2 / sum(object$pca$group$sdev^2)
    object$pca$loading$explained$group <- object$pca$score$explained$group
  }
  
  return(object)
}

cleanALASCA<- function(object){
  # Clean up a copy
  if(object$separateTimeAndGroup){
    # Time effect
    object$ALASCA$loading$time <- object$pca$loading$time[!duplicated(object$pca$loading$time),]
    object$ALASCA$loading$time <- reshape2::melt(object$ALASCA$loading$time, id.vars = "covars")
    colnames(object$ALASCA$loading$time) <- c("covars","PC","loading")
    object$ALASCA$loading$time$PC <- as.numeric(gsub("PC","", object$ALASCA$loading$time$PC))
    
    object$ALASCA$score$time <- object$pca$score$time[!duplicated(object$pca$score$time),]
    object$ALASCA$score$time <- reshape2::melt(object$ALASCA$score$time, id.vars = c("time"))
    colnames(object$ALASCA$score$time) <- c("time","PC","score")
    object$ALASCA$score$time$PC <- as.numeric(gsub("PC","", object$ALASCA$score$time$PC))
    
    object$ALASCA$score$explained$time <- object$pca$score$explained$time
    object$ALASCA$loading$explained$time <- object$pca$loading$explained$time
    
    # Group effect
    object$ALASCA$loading$group <- object$pca$loading$group[!duplicated(object$pca$loading$group),]
    object$ALASCA$loading$group <- reshape2::melt(object$ALASCA$loading$group, id.vars = "covars")
    colnames(object$ALASCA$loading$group) <- c("covars","PC","loading")
    object$ALASCA$loading$group$PC <- as.numeric(gsub("PC","", object$ALASCA$loading$group$PC))
    
    object$ALASCA$score$group <- object$pca$score$group[!duplicated(object$pca$score$group),]
    object$ALASCA$score$group <- reshape2::melt(object$ALASCA$score$group, id.vars = c("time","group"))
    colnames(object$ALASCA$score$group) <- c("time","group","PC","score")
    object$ALASCA$score$group$PC <- as.numeric(gsub("PC","", object$ALASCA$score$group$PC))
    
    object$ALASCA$score$explained$group <- object$pca$score$explained$group
    object$ALASCA$loading$explained$group <- object$pca$loading$explained$group
    
  }else{
    object$ALASCA$loading$time <- object$pca$loading$time[!duplicated(object$pca$loading$time),]
    object$ALASCA$loading$time <- reshape2::melt(object$ALASCA$loading$time, id.vars = "covars")
    colnames(object$ALASCA$loading$time) <- c("covars","PC","loading")
    object$ALASCA$loading$time$PC <- as.numeric(gsub("PC","", object$ALASCA$loading$time$PC))
    
    object$ALASCA$score$time <- object$pca$score$time[!duplicated(object$pca$score$time),]
    object$ALASCA$score$time <- reshape2::melt(object$ALASCA$score$time, id.vars = c("time","group"))
    colnames(object$ALASCA$score$time) <- c("time","group","PC","score")
    object$ALASCA$score$time$PC <- as.numeric(gsub("PC","", object$ALASCA$score$time$PC))
    
    object$ALASCA$score$explained$time <- object$pca$score$explained$time
    object$ALASCA$loading$explained$time <- object$pca$loading$explained$time
  }
  return(object)
}