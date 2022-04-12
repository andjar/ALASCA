#' Perform PCA
#'
#' This function performs PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
do_pca <- function(object) {
  object[["pca"]][["time"]] <- prcomp(
    object$effect_matrix[object$effect_matrix$comp == "TIME",
                         seq_len(ncol(object$effect_matrix) - 1)],
    scale = FALSE,
    center = !object$scale_function.center
    )
  if (object$separate_time_and_group) {
    object[["pca"]][["group"]] <- prcomp(
      object$effect_matrix[object$effect_matrix$comp == "GROUP",
                           seq_len(ncol(object$effect_matrix) - 1)],
      scale = FALSE,
      center = !object$scale_function.center
    )
  }
  return(object)
}

#' Perform PCA for Limm-PCA
#'
#' This function performs PCA before "the real PCA"
#'
#' @param object An ALASCA object
#' @return An ALASCA object
reduce_dimensions <- function(object){
  
  log4r::debug(object$log, "Starting reduce_dimensions")

  wide_data <- dcast(data = object$df, ...~variable
                     )
  if (object$do_debug) currentTs <- Sys.time()
  temp_pca_values <- object$function.pca(
    wide_data[, .SD, .SDcols = -object$all_formula_terms],
    center = !object$scale_function.center
  )
  
  object$reduce_dimensions.explanatory_power <- temp_pca_values$sdev^2 / sum(temp_pca_values$sdev^2)
  
  # Remove surplus columns
  if (is.null(object$reduce_dimensions.nComps)) {
    object$reduce_dimensions.nComps <- which(cumsum(object$reduce_dimensions.explanatory_power) >= object$reduce_dimensions.limit)[1]
    log4r::info(object$log, paste("Keeping",
                                  object$reduce_dimensions.nComps,
                                  "components from initial PCA, explaining",
                                  100*cumsum(object$reduce_dimensions.explanatory_power)[object$reduce_dimensions.nComps],
                                  "% of variation"))
  }
  if(ncol(temp_pca_values$rotation) > object$reduce_dimensions.nComps){
    temp_pca_values$rotation <- temp_pca_values$rotation[,-c((object$reduce_dimensions.nComps+1):ncol(temp_pca_values$rotation))]
    temp_pca_values$x <- temp_pca_values$x[, -c((object$reduce_dimensions.nComps+1):ncol(temp_pca_values$x))]
  }
  
  # Check if the pca model needs reflection to better fit the main model
  if (object$do_debug) currentTs <- Sys.time()
  for (i in seq_len(ncol(temp_pca_values$rotation))) {
    V1 <- sum((temp_pca_values$rotation[,i] - object$Limm$main$pca$rotation[,i])^2)
    V2 <- sum((-temp_pca_values$rotation[,i] - object$Limm$main$pca$rotation[,i])^2)
    if(V2 < V1){
      temp_pca_values$rotation[,i] = -temp_pca_values$rotation[,i]
      temp_pca_values$x[,i] = -temp_pca_values$x[,i]
    }
  }
  log4r::info(object$log, "Removing surplus PCs")
  
  object$Limm$loadings <- temp_pca_values$rotation
  object$Limm$pca <- temp_pca_values
  object$Limm$df <- object$df
  object$df <- melt(data = cbind(wide_data[, .SD, .SDcols = object$all_formula_terms], temp_pca_values$x), id.vars = object$all_formula_terms, variable.factor = FALSE)
  object$variablelist <- unique(object$df$variable)
  object$stratification_vector <- object$df[, get(object$stratification_column)]
  log4r::info(object$log, "Completed reduce_dimensions")
  return(object)
}

#' Clean the PCA data
#'
#' This function makes the pca output more useful
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
clean_pca <- function(object) {
  log4r::debug(object$log, "Starting clean_pca")
  
  # Clean scores ----
  PC_time <- as.data.frame(object$pca$time$x)
  PC_time$time <- object$parts$time
  if (object$separate_time_and_group) {
    PC_time$group <- object$grouplist[1]
  }else{
    PC_time$group <- object$parts$group
  }
  object$pca$score$time <- setDT(PC_time[!duplicated(paste(PC_time$time, PC_time$group)),])
  setkey(object$pca$score$time, time, group)
  object$pca$score$explained$time <- object$pca$time$sdev^2 / sum(object$pca$time$sdev^2)
  
  # Clean loadings ----
  object$pca$loading$time <- setDT(as.data.frame(object$pca$time$rotation), keep.rownames="covars")
  object$pca$loading$explained$time <- object$pca$score$explained$time
  
  if(object$reduce_dimensions){
    # Loadings must be back-transformed
    object$pca$loading$time <- setDT(
      as.data.frame(
        as.matrix(object$Limm$loadings) %*% as.matrix(object$pca$loading$time[order(as.numeric(substr(covars, 3, nchar(covars)))), !"covars"])
        ),
      keep.rownames="covars")
    object$variablelist <- unique(object$pca$loading$time$covars)
    object$regression_coefficients <- rbindlist(
      lapply(unique(object$regression_coefficients$variable), function(x){
        setDT(
          data.frame(
            variable = x,
            pvalue = NA,
            estimate = as.matrix(object$Limm$loadings) %*% as.matrix(object$regression_coefficients[variable == x & order(as.numeric(substr(covar, 3, nchar(covar)))), "estimate"])
          ),
        keep.rownames="covar")
      })
    )
  }
  
  setkey(object$pca$loading$time, covars)
  
  PCloading <- paste0("PC", get_relevant_pcs(object = object, effect = "time"))
  for (i in PCloading) {
    # Ensure that the highest loading has positive sign
    nVar <- object$pca$loading$time[, .I[which.max(abs(get(i)))]]
    sVar <- object$pca$loading$time[nVar, sign(get(i))]
    set(object$pca$loading$time, j = (i), value = object$pca$loading$time[, get(i)] * sVar)
    set(object$pca$score$time, j = (i), value = object$pca$score$time[, get(i)] * sVar)
  }

  if (object$separate_time_and_group) {
    # Clean scores ----
    PC_group <- as.data.frame(object$pca$group$x)
    PC_group$time <- rownames(PC_group)
    PC_group$group <- object$parts$group
    PC_group$time <- object$parts$time
    object$pca$score$group <- setDT(PC_group[!duplicated(paste(PC_group$time, PC_group$group)),])
    setkey(object$pca$score$group, time, group)
    object$pca$score$explained$group <- object$pca$group$sdev^2 / sum(object$pca$group$sdev^2)
    
    # Clean loadings ----
    object$pca$loading$group <- setDT(as.data.frame(object$pca$group$rotation), keep.rownames="covars")
    object$pca$loading$explained$group <- object$pca$score$explained$group
    
    if(object$reduce_dimensions){
      # Loadings must be back-transformed
      object$pca$loading$group <- setDT(
        as.data.frame(
          as.matrix(object$Limm$loadings) %*% as.matrix(object$pca$loading$group[order(as.numeric(substr(covars, 3, nchar(covars)))), !"covars"])
          ),
        keep.rownames="covars")
    }
    
    setkey(object$pca$loading$group, covars)

    PCloading <- paste0("PC", get_relevant_pcs(object = object, effect = "group"))
    for (i in PCloading) {
      # Ensure that the highest loading has positive sign
      nVar <- object$pca$loading$group[, .I[which.max(abs(get(i)))]]
      sVar <- object$pca$loading$time[nVar, sign(get(i))]
      set(object$pca$loading$group, j = (i), value = object$pca$loading$group[, get(i)] * sVar)
      set(object$pca$score$group, j = (i), value = object$pca$score$group[, get(i)] * sVar)
    }
  }

  log4r::debug(object$log, "Completed clean_pca")
  
  return(object)
}

clean_alasca <- function(object) {

  log4r::debug(object$log, "Starting clean_alasca")
  # We need to create new group names for the combined group and keep_terms
  if (object$keep_terms != "") {
    if (object$separate_time_and_group) {
      object$grouplist <- unique(object$pca$score$group$group)
    } else {
      object$grouplist <- unique(object$pca$score$time$group)
    }
  }

  # Clean up a copy
  # Time effect
  object$ALASCA$loading$time <- melt(object$pca$loading$time, id.vars = "covars", variable.factor = FALSE)
  colnames(object$ALASCA$loading$time) <- c("covars", "PC", "loading")
  object$ALASCA$loading$time[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

  object$ALASCA$score$time <- object$pca$score$time
  if (object$separate_time_and_group) {
    object$ALASCA$score$time$group <- object$grouplist[1]
  }
  object$ALASCA$score$time <- melt(object$ALASCA$score$time, id.vars = c("time", "group"), variable.factor = FALSE)
  colnames(object$ALASCA$score$time) <- c("time", "group", "PC", "score")
  object$ALASCA$score$time[, time := factor(time, levels = object$timelist), ]
  object$ALASCA$score$time[, group := factor(group, levels = object$grouplist), ]
  object$ALASCA$score$time[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

  object$ALASCA$score$explained$time <- object$pca$score$explained$time
  object$ALASCA$loading$explained$time <- object$pca$loading$explained$time

  if (object$separate_time_and_group) {
    # Group effect
    object$ALASCA$loading$group <- melt(object$pca$loading$group, id.vars = "covars", variable.factor = FALSE)
    colnames(object$ALASCA$loading$group) <- c("covars", "PC", "loading")
    object$ALASCA$loading$group[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

    object$ALASCA$score$group <- object$pca$score$group
    object$ALASCA$score$group <- melt(object$ALASCA$score$group, id.vars = c("time", "group"), variable.factor = FALSE)
    colnames(object$ALASCA$score$group) <- c("time", "group", "PC", "score")
    object$ALASCA$score$group[, time := factor(time, levels = object$timelist), ]
    object$ALASCA$score$group[, group := factor(group, levels = object$grouplist), ]
    object$ALASCA$score$group[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

    object$ALASCA$score$explained$group <- object$pca$score$explained$group
    object$ALASCA$loading$explained$group <- object$pca$loading$explained$group
  }
  
  log4r::debug(object$log, "Finished clean_alasca")

  return(object)
}
