#' Perform PCA
#'
#' This function performs PCA
#'
#' @param object An RMASCA object to be sanitized
#' @return An RMASCA object
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
#' @param object An RMASCA object to be sanitized
#' @return An RMASCA object
cleanPCA <- function(object){

  loading_time <- as.data.frame(object$pca$time$rotation[,])
  loading_time$covars <- rownames(loading_time)
  object$RMASCA$loading$time <- loading_time

  PC_time <- as.data.frame(object$pca$time$x[,])
  PC_time$time <- object$parts$time
  if(!object$separateTimeAndGroup){
    PC_time$group <- object$parts$group
  }
  object$RMASCA$score$time <- PC_time
  object$RMASCA$score$explained$time <- object$pca$time$sdev^2 / sum(object$pca$time$sdev^2)
  object$RMASCA$loading$explained$time <- object$RMASCA$score$explained$time

  if(object$separateTimeAndGroup){
    loading_group <- as.data.frame(object$pca$group$rotation[,])
    loading_group$covars <- rownames(loading_group)
    object$RMASCA$loading$group <- loading_group
    PC_group <- as.data.frame(object$pca$group$x[,])
    PC_group$time <- rownames(PC_group)
    PC_group$group <- object$parts$group
    PC_group$time <- object$parts$time
    object$RMASCA$score$group <- PC_group
    object$RMASCA$score$explained$group <- object$pca$group$sdev^2 / sum(object$pca$group$sdev^2)
    object$RMASCA$loading$explained$group <- object$RMASCA$score$explained$group
  }
  return(object)
}
