#' Perform PCA
#'
#' This function performs PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
doPCA <- function(object) {
  object$pca$time <- prcomp(
    object$effect.matrix[object$effect.matrix$comp == "TIME",
                         seq_len(ncol(object$effect.matrix) - 1)],
    scale = FALSE,
    center = TRUE)
  if (object$separateTimeAndGroup) {
    object$pca$group <- prcomp(
      object$effect.matrix[object$effect.matrix$comp == "GROUP",
                           seq_len(ncol(object$effect.matrix) - 1)],
      scale = FALSE,
      center = TRUE)
  }
  return(object)
}

#' Perform PCA for Limm-PCA
#'
#' This function performs PCA before "the real PCA"
#'
#' @param object An ALASCA object
#' @return An ALASCA object
doLimmPCA <- function(object){
  wide_data <- dcast(data = object$df, as.formula(paste(paste(object$allFormulaTerms, collapse = " + "), "~ variable")))
  if (object$doDebug) currentTs <- Sys.time()
  temp_pca_values <- prcomp(
    wide_data[, (length(object$allFormulaTerms)+1):ncol(wide_data)],
    scale = FALSE,
    center = TRUE)
  if (object$doDebug) cat("* First PCA:", Sys.time() - currentTs, "s\n")
  
  # Remove surplus columns
  if(ncol(temp_pca_values$rotation) > object$limm.nComps){
    temp_pca_values$rotation <- temp_pca_values$rotation[,-c((object$limm.nComps+1):ncol(temp_pca_values$rotation))]
    temp_pca_values$x <- temp_pca_values$x[, -c((object$limm.nComps+1):ncol(temp_pca_values$x))]
  }
  
  object$Limm$pca <- temp_pca_values
  
  # Check if the pca model needs reflection to better fit the main model
  if (object$doDebug) currentTs <- Sys.time()
  for (i in seq_len(ncol(temp_pca_values$rotation))) {
    V1 <- sum((temp_pca_values$rotation[,i] - object$Limm$main$pca$rotation[,i])^2)
    V2 <- sum((-temp_pca_values$rotation[,i] - object$Limm$main$pca$rotation[,i])^2)
    if(V2 < V1){
      temp_pca_values$rotation[,i] = -temp_pca_values$rotation[,i]
      temp_pca_values$x[,i] = -temp_pca_values$x[,i]
    }
  }
  if (object$doDebug) cat("* Remove surplus PCs:", Sys.time() - currentTs, "s\n")
  
  object$Limm$loadings <- temp_pca_values$rotation
  object$Limm$df <- object$df
  object$df <- melt(data = cbind(wide_data[, .SD, .SDcols = object$allFormulaTerms], temp_pca_values$x), id.vars = object$allFormulaTerms, variable.factor = FALSE)
  object$variablelist <- unique(object$df$variable)
  object$stratificationVector <- object$df[, get(object$stratificationColumn)]
  return(object)
}

#' Clean the PCA data
#'
#' This function makes the pca output more useful
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
cleanPCA <- function(object) {
  loading_time <- as.data.frame(object$pca$time$rotation)
  loading_time$covars <- rownames(loading_time)
  object$pca$loading$time <- setDT(loading_time)
  setkey(object$pca$loading$time, covars)

  PC_time <- as.data.frame(object$pca$time$x)
  PC_time$time <- object$parts$time
  if (object$separateTimeAndGroup) {
    PC_time$group <- object$grouplist[1]
  }else{
    PC_time$group <- object$parts$group
  }
  object$pca$score$time <- setDT(PC_time[!duplicated(paste(PC_time$time, PC_time$group)),])
  setkey(object$pca$score$time, time, group)
  object$pca$score$explained$time <- object$pca$time$sdev^2 / sum(object$pca$time$sdev^2)
  object$pca$loading$explained$time <- object$pca$score$explained$time

  PCloading <- getRelevantPCs(object = object, effect = "time")
  for (i in PCloading) {
    # Ensure that the highest loading has positive sign
    nVar <- which(abs(object$pca$loading$time[, ..i]) == max(abs(object$pca$loading$time[, ..i])))[1]
    mVar <- object$pca$loading$time$covars[[nVar]]
    sVar <- sign(object$pca$loading$time[[nVar, i]])
    object$pca$loading$time[, (i) := sVar * .SD, .SDcols = i ]
    object$pca$score$time[, (i) := sVar * .SD, .SDcols = i ]
  }

  if (object$separateTimeAndGroup) {
    loading_group <- as.data.frame(object$pca$group$rotation)
    loading_group$covars <- rownames(loading_group)
    object$pca$loading$group <- setDT(loading_group)
    setkey(object$pca$loading$group, covars)
    
    PC_group <- as.data.frame(object$pca$group$x)
    PC_group$time <- rownames(PC_group)
    PC_group$group <- object$parts$group
    PC_group$time <- object$parts$time
    object$pca$score$group <- setDT(PC_group[!duplicated(paste(PC_group$time, PC_group$group)),])
    setkey(object$pca$score$group, time, group)
    object$pca$score$explained$group <- object$pca$group$sdev^2 / sum(object$pca$group$sdev^2)
    object$pca$loading$explained$group <- object$pca$score$explained$group

    PCloading <- getRelevantPCs(object = object, effect = "group")
    for (i in PCloading) {
      # Ensure that the highest loading has positive sign
      nVar <- which(abs(object$pca$loading$group[, ..i]) == max(abs(object$pca$loading$group[, ..i])))[1]
      mVar <- object$pca$loading$group$covars[[nVar]]
      sVar <- sign(object$pca$loading$group[[nVar, i]])
      object$pca$loading$group[, (i) := sVar * .SD, .SDcols = i ]
      object$pca$score$group[, (i) := sVar * .SD, .SDcols = i ]
    }
  }

  return(object)
}

cleanALASCA <- function(object) {
  # We need to create new group names for the combined group and keepTerms
  if (object$keepTerms != "") {
    if (object$separateTimeAndGroup) {
      object$grouplist <- unique(object$pca$score$group$group)
    } else {
      object$grouplist <- unique(object$pca$score$time$group)
    }
  }

  # Clean up a copy
  # Time effect
  object$ALASCA$loading$time <- object$pca$loading$time
  object$ALASCA$loading$time <- melt(object$ALASCA$loading$time, id.vars = "covars")
  colnames(object$ALASCA$loading$time) <- c("covars", "PC", "loading")
  if(object$method %in% c("Limm", "Lim")){
    object$ALASCA$loading$time <- rbindlist(lapply(unique(object$ALASCA$loading$time$PC), function(selPC){
      ref <- object$ALASCA$loading$time[PC == selPC,]
      data.frame(
        covars = rownames(object$Limm$loadings),
        PC = selPC,
        loading = rowSums(object$Limm$loadings*ref$loading[match(colnames(object$Limm$loadings), ref$covars)][col(object$Limm$loadings)])
      )
    }))
    object$variablelist <- unique(object$ALASCA$loading$time$covars)
  }
  object$ALASCA$loading$time[, PC := as.numeric(gsub("PC", "", PC)), ]

  object$ALASCA$score$time <- object$pca$score$time
  if (object$separateTimeAndGroup) {
    object$ALASCA$score$time$group <- object$grouplist[1]
  }
  object$ALASCA$score$time <- melt(object$ALASCA$score$time, id.vars = c("time", "group"))
  colnames(object$ALASCA$score$time) <- c("time", "group", "PC", "score")
  object$ALASCA$score$time[, time := factor(time, levels = object$timelist), ]
  object$ALASCA$score$time[, group := factor(group, levels = object$grouplist), ]
  object$ALASCA$score$time[, PC := as.numeric(gsub("PC", "", PC)), ]

  object$ALASCA$score$explained$time <- object$pca$score$explained$time
  object$ALASCA$loading$explained$time <- object$pca$loading$explained$time

  if (object$separateTimeAndGroup) {
    # Group effect
    object$ALASCA$loading$group <- object$pca$loading$group
    object$ALASCA$loading$group <- melt(object$ALASCA$loading$group, id.vars = "covars")
    colnames(object$ALASCA$loading$group) <- c("covars", "PC", "loading")
    object$ALASCA$loading$group[, PC := as.numeric(gsub("PC", "", PC)), ]

    object$ALASCA$score$group <- object$pca$score$group
    object$ALASCA$score$group <- melt(object$ALASCA$score$group, id.vars = c("time", "group"))
    colnames(object$ALASCA$score$group) <- c("time", "group", "PC", "score")
    object$ALASCA$score$group[, time := factor(time, levels = object$timelist), ]
    object$ALASCA$score$group[, group := factor(group, levels = object$grouplist), ]
    if(object$method %in% c("Limm", "Lim")){
      object$ALASCA$loading$group <- rbindlist(lapply(unique(object$ALASCA$loading$group$PC), function(selPC){
        ref <- object$ALASCA$loading$group[PC == selPC,]
        data.frame(
          covars = rownames(object$Limm$loadings),
          PC = selPC,
          loading = rowSums(object$Limm$loadings*ref$loading[match(colnames(object$Limm$loadings), ref$covars)][col(object$Limm$loadings)])
        )
      }))
    }
    object$ALASCA$score$group[, PC := as.numeric(gsub("PC", "", PC)), ]

    object$ALASCA$score$explained$group <- object$pca$score$explained$group
    object$ALASCA$loading$explained$group <- object$pca$loading$explained$group
  }

  return(object)
}
