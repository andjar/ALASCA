#' Perform PCA
#'
#' This function performs PCA
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
do_pca <- function() {
  self[["pca"]][["time"]] <- prcomp(
    self$effect_matrix[self$effect_matrix$comp == "TIME",
                         seq_len(ncol(self$effect_matrix) - 1)],
    scale = FALSE,
    center = !self$scale_function.center
    )
  if (self$separate_time_and_group) {
    self[["pca"]][["group"]] <- prcomp(
      self$effect_matrix[self$effect_matrix$comp == "GROUP",
                           seq_len(ncol(self$effect_matrix) - 1)],
      scale = FALSE,
      center = !self$scale_function.center
    )
  }
  #invisible(self)
}

#' Perform PCA for Limm-PCA
#'
#' This function performs PCA before "the real PCA"
#'
#' @param object An ALASCA object
#' @return An ALASCA object
do_reduce_dimensions <- function(){
  
  log4r::debug(self$log, "Starting reduce_dimensions")

  wide_data <- dcast(data = self$df, ...~variable
                     )
  if (self$do_debug) currentTs <- Sys.time()
  temp_pca_values <- self$function.pca(
    wide_data[, .SD, .SDcols = -self$all_formula_terms],
    center = !self$scale_function.center
  )
  
  self$reduce_dimensions.explanatory_power <- temp_pca_values$sdev^2 / sum(temp_pca_values$sdev^2)
  
  # Remove surplus columns
  if (is.null(self$reduce_dimensions.nComps)) {
    self$reduce_dimensions.nComps <- which(cumsum(self$reduce_dimensions.explanatory_power) >= self$reduce_dimensions.limit)[1]
    log4r::info(self$log, paste("Keeping",
                                  self$reduce_dimensions.nComps,
                                  "components from initial PCA, explaining",
                                  100*cumsum(self$reduce_dimensions.explanatory_power)[self$reduce_dimensions.nComps],
                                  "% of variation"))
  }
  if(ncol(temp_pca_values$rotation) > self$reduce_dimensions.nComps){
    temp_pca_values$rotation <- temp_pca_values$rotation[,-c((self$reduce_dimensions.nComps+1):ncol(temp_pca_values$rotation))]
    temp_pca_values$x <- temp_pca_values$x[, -c((self$reduce_dimensions.nComps+1):ncol(temp_pca_values$x))]
  }
  
  # Check if the pca model needs reflection to better fit the main model
  if (self$do_debug) currentTs <- Sys.time()
  for (i in seq_len(ncol(temp_pca_values$rotation))) {
    V1 <- sum((temp_pca_values$rotation[,i] - self$Limm$main$pca$rotation[,i])^2)
    V2 <- sum((-temp_pca_values$rotation[,i] - self$Limm$main$pca$rotation[,i])^2)
    if(V2 < V1){
      temp_pca_values$rotation[,i] = -temp_pca_values$rotation[,i]
      temp_pca_values$x[,i] = -temp_pca_values$x[,i]
    }
  }
  log4r::info(self$log, "Removing surplus PCs")
  
  self$Limm$loadings <- temp_pca_values$rotation
  self$Limm$pca <- temp_pca_values
  self$Limm$df <- self$df
  self$df <- melt(data = cbind(wide_data[, .SD, .SDcols = self$all_formula_terms], temp_pca_values$x), id.vars = self$all_formula_terms, variable.factor = FALSE)
  self$variablelist <- unique(self$df$variable)
  self$stratification_vector <- self$df[, get(self$stratification_column)]
  log4r::info(self$log, "Completed reduce_dimensions")
  #invisible(self)
}

#' Clean the PCA data
#'
#' This function makes the pca output more useful
#'
#' @param object An ALASCA object to be sanitized
#' @return An ALASCA object
clean_pca <- function() {
  log4r::debug(self$log, "Starting clean_pca")
  
  # Clean scores ----
  PC_time <- as.data.frame(self$pca$time$x)
  PC_time$time <- self$parts$time
  if (self$separate_time_and_group) {
    PC_time$group <- self$grouplist[1]
  }else{
    PC_time$group <- self$parts$group
  }
  self$pca$score$time <- setDT(PC_time[!duplicated(paste(PC_time$time, PC_time$group)),])
  setkey(self$pca$score$time, time, group)
  self$pca$score$explained$time <- self$pca$time$sdev^2 / sum(self$pca$time$sdev^2)
  
  # Clean loadings ----
  self$pca$loading$time <- setDT(as.data.frame(self$pca$time$rotation), keep.rownames="covars")
  self$pca$loading$explained$time <- self$pca$score$explained$time
  
  if(self$reduce_dimensions){
    # Loadings must be back-transformed
    self$pca$loading$time <- setDT(
      as.data.frame(
        as.matrix(self$Limm$loadings) %*% as.matrix(self$pca$loading$time[order(as.numeric(substr(covars, 3, nchar(covars)))), !"covars"])
        ),
      keep.rownames="covars")
    self$variablelist <- unique(self$pca$loading$time$covars)
    self$regression_coefficients <- rbindlist(
      lapply(unique(self$regression_coefficients$variable), function(x){
        setDT(
          data.frame(
            variable = x,
            pvalue = NA,
            estimate = as.matrix(self$Limm$loadings) %*% as.matrix(self$regression_coefficients[variable == x & order(as.numeric(substr(covar, 3, nchar(covar)))), "estimate"])
          ),
        keep.rownames="covar")
      })
    )
  }
  
  setkey(self$pca$loading$time, covars)
  
  PCloading <- paste0("PC", get_relevant_pcs(self = self, effect = "time"))
  for (i in PCloading) {
    # Ensure that the highest loading has positive sign
    nVar <- self$pca$loading$time[, .I[which.max(abs(get(i)))]]
    sVar <- self$pca$loading$time[nVar, sign(get(i))]
    set(self$pca$loading$time, j = (i), value = self$pca$loading$time[, get(i)] * sVar)
    set(self$pca$score$time, j = (i), value = self$pca$score$time[, get(i)] * sVar)
  }

  if (self$separate_time_and_group) {
    # Clean scores ----
    PC_group <- as.data.frame(self$pca$group$x)
    PC_group$time <- rownames(PC_group)
    PC_group$group <- self$parts$group
    PC_group$time <- self$parts$time
    self$pca$score$group <- setDT(PC_group[!duplicated(paste(PC_group$time, PC_group$group)),])
    setkey(self$pca$score$group, time, group)
    self$pca$score$explained$group <- self$pca$group$sdev^2 / sum(self$pca$group$sdev^2)
    
    # Clean loadings ----
    self$pca$loading$group <- setDT(as.data.frame(self$pca$group$rotation), keep.rownames="covars")
    self$pca$loading$explained$group <- self$pca$score$explained$group
    
    if(self$reduce_dimensions){
      # Loadings must be back-transformed
      self$pca$loading$group <- setDT(
        as.data.frame(
          as.matrix(self$Limm$loadings) %*% as.matrix(self$pca$loading$group[order(as.numeric(substr(covars, 3, nchar(covars)))), !"covars"])
          ),
        keep.rownames="covars")
    }
    
    setkey(self$pca$loading$group, covars)

    PCloading <- paste0("PC", get_relevant_pcs(self = self, effect = "group"))
    for (i in PCloading) {
      # Ensure that the highest loading has positive sign
      nVar <- self$pca$loading$group[, .I[which.max(abs(get(i)))]]
      sVar <- self$pca$loading$time[nVar, sign(get(i))]
      set(self$pca$loading$group, j = (i), value = self$pca$loading$group[, get(i)] * sVar)
      set(self$pca$score$group, j = (i), value = self$pca$score$group[, get(i)] * sVar)
    }
  }

  log4r::debug(self$log, "Completed clean_pca")
  
  #invisible(self)
}

clean_alasca <- function() {

  log4r::debug(self$log, "Starting clean_alasca")
  # We need to create new group names for the combined group and keep_terms
  if (self$keep_terms != "") {
    if (self$separate_time_and_group) {
      self$grouplist <- unique(self$pca$score$group$group)
    } else {
      self$grouplist <- unique(self$pca$score$time$group)
    }
  }

  # Clean up a copy
  # Time effect
  self$ALASCA$loading$time <- melt(self$pca$loading$time, id.vars = "covars", variable.factor = FALSE)
  colnames(self$ALASCA$loading$time) <- c("covars", "PC", "loading")
  self$ALASCA$loading$time[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

  self$ALASCA$score$time <- self$pca$score$time
  if (self$separate_time_and_group) {
    self$ALASCA$score$time$group <- self$grouplist[1]
  }
  self$ALASCA$score$time <- melt(self$ALASCA$score$time, id.vars = c("time", "group"), variable.factor = FALSE)
  colnames(self$ALASCA$score$time) <- c("time", "group", "PC", "score")
  self$ALASCA$score$time[, time := factor(time, levels = self$timelist), ]
  self$ALASCA$score$time[, group := factor(group, levels = self$grouplist), ]
  self$ALASCA$score$time[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

  self$ALASCA$score$explained$time <- self$pca$score$explained$time
  self$ALASCA$loading$explained$time <- self$pca$loading$explained$time

  if (self$separate_time_and_group) {
    # Group effect
    self$ALASCA$loading$group <- melt(self$pca$loading$group, id.vars = "covars", variable.factor = FALSE)
    colnames(self$ALASCA$loading$group) <- c("covars", "PC", "loading")
    self$ALASCA$loading$group[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

    self$ALASCA$score$group <- self$pca$score$group
    self$ALASCA$score$group <- melt(self$ALASCA$score$group, id.vars = c("time", "group"), variable.factor = FALSE)
    colnames(self$ALASCA$score$group) <- c("time", "group", "PC", "score")
    self$ALASCA$score$group[, time := factor(time, levels = self$timelist), ]
    self$ALASCA$score$group[, group := factor(group, levels = self$grouplist), ]
    self$ALASCA$score$group[, PC := as.integer(substr(PC, 3, nchar(PC))), ]

    self$ALASCA$score$explained$group <- self$pca$score$explained$group
    self$ALASCA$loading$explained$group <- self$pca$loading$explained$group
  }
  
  log4r::debug(self$log, "Finished clean_alasca")

  #invisible(self)
}
