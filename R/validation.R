#' Validate the ALASCA model LMM
#'
#' This function performs leave-one-out robustness testing of your ALASCA model. If you didn't specify the number of runs `n_validation_runs` when initializing the model (see \code{\link{ALASCA}}), you can do it by running for example `model$n_validation_runs <- 100` prior to calling `validate`. Your dataset is divided into `n_validation_folds` partitions, keeping group proportions, and one of these are left out. `n_validation_folds` is set the same way as  `n_validation_runs`.
#'
#' @param object An ALASCA object
#' @param participant_column The name of the column containing participant identifier. Needed if not set during initialization of the model.
#' @param validate_regression Whether to validate regression models
#' @return An ALASCA object
#'
#' @examples
#' load("PE.Rdata")
#' model$n_validation_runs <- 10
#' model.val <- validate(model, participant_column = "ID")
#' @export
validate <- function(object, participant_column = FALSE, validate_regression = FALSE) {
  if (object$validate) {
    # stop("The object has already been validated")
  }
  object$validate <- TRUE
  if (validate_regression) {
    object$validate_regression <- TRUE
  }

  if (object$method == "LMM") {
    if (participant_column == FALSE) {
      if (object$participant_column == FALSE) {
        log4r::error(object$log, "You need to specify the column containing participant id in `participant_column`")
        stop()
      }
    } else {
      object$participant_column <- participant_column
    }
  }

  log4r::info(object$log, "Starting validation")

  start.time.all <- Sys.time()

  if (any(is.na(object$validation_ids))) {
    # Generate random samples
    if (object$save_to_disk) {
      limPC_time <- get_relevant_pcs(object = object, effect = "time")
      if (object$separate_time_and_group) {
        limPC_group <- get_relevant_pcs(object = object, effect = "group")
      }
      temp_object <- lapply(seq_len(object$n_validation_runs), FUN = function(ii) {
        log4r::info(object$log, paste0("- Run ", ii, " of ", object$n_validation_runs))
        start.time.this <- Sys.time()

        # Make resampled model
        temp_object <- prepare_validation_run(object)

        # Rotate new loadings/scores to the original model
        if (object$optimize_score) {
          temp_object <- rotate_matrix_optimize_score(object = temp_object, target = object)
        } else {
          temp_object <- rotate_matrix(object = temp_object, target = object)
        }

        temp_object <- clean_alasca(temp_object)

        # Save to disk
        fname <- paste0("val_", ii)
        loading <- data.frame(temp_object$ALASCA$loading$time[temp_object$ALASCA$loading$time$PC %in% limPC_time, ], model = ii)
        # loading$covars <- object$numvariablelist[loading$covars]
        DBI::dbWriteTable(object$db.con, "time.loading", loading, append = T, overwrite = F)
        scores <- data.frame(temp_object$ALASCA$score$time[temp_object$ALASCA$score$time$PC %in% limPC_time, ], model = ii)
        # scores$time <- object$numtimelist[scores$time]
        # if(!object$separate_time_and_group){
        #  scores$group <- object$numgrouplist[scores$group]
        # }
        DBI::dbWriteTable(object$db.con, "time.score", scores, append = T, overwrite = F)
        if (nrow(get_covars(object)) > 0) DBI::dbWriteTable(object$db.con, "covars", data.frame(get_covars(temp_object), model = ii), append = T, overwrite = F)

        if (object$separate_time_and_group) {
          loading <- data.frame(temp_object$ALASCA$loading$group[temp_object$ALASCA$loading$group$PC %in% limPC_group, ], model = ii)
          # loading$covars <- object$numvariablelist[loading$covars]
          DBI::dbWriteTable(object$db.con, "group.loading", loading, append = T, overwrite = F)

          scores <- data.frame(temp_object$ALASCA$score$group[temp_object$ALASCA$score$group$PC %in% limPC_group, ], model = ii)
          # scores$time <- object$numtimelist[scores$time]
          # scores$group <- object$numgrouplist[scores$group]
          DBI::dbWriteTable(object$db.con, "group.score", scores, append = T, overwrite = F)
        }
        if (object$validate_regression) {
          DBI::dbWriteTable(object$db.con, "model_prediction", data.frame(temp_object$model_prediction, model = ii), append = T, overwrite = F)
        }

        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(object$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((object$n_validation_runs - ii) * time_all, 2), " seconds"))
        fname
      })
    } else {
      temp_object <- lapply(seq_len(object$n_validation_runs), FUN = function(ii) {
        
        log4r::info(object$log, paste0("- Run ", ii, " of ", object$n_validation_runs))
        start.time.this <- Sys.time()

        # Make resampled model
        temp_object <- prepare_validation_run(object)
        if (nrow(object$df) > 100000) gc()

        # Rotate new loadings/scores to the original model
        if (object$optimize_score) {
          temp_object <- rotate_matrix_optimize_score(object = temp_object, target = object)
        } else {
          temp_object <- rotate_matrix(object = temp_object, target = object)
        }
        temp_object <- clean_alasca(temp_object)

        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(object$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((object$n_validation_runs - ii) * time_all, 2), " seconds"))
        temp_object
      })
    }
  } else {
    # Reuse previous samples
    log4r::info(object$log, "Using predefined samples")

    if (object$save_to_disk) {
      limPC_time <- get_relevant_pcs(object = object, effect = "time")
      if (object$separate_time_and_group) {
        limPC_group <- get_relevant_pcs(object = object, effect = "group")
      }
      temp_object <- lapply(seq_len(object$n_validation_runs), FUN = function(ii) {
        log4r::info(object$log, paste0("- Run ", ii, " of ", object$n_validation_runs))
        start.time.this <- Sys.time()

        # Make resampled model
        temp_object <- prepare_validation_run(object, runN = ii)
        gc()

        # Rotate new loadings/scores to the original model
        if (object$optimize_score) {
          temp_object <- rotate_matrix_optimize_score(object = temp_object, target = object)
        } else {
          temp_object <- rotate_matrix(object = temp_object, target = object)
        }

        temp_object <- clean_alasca(temp_object)

        # Save to disk
        fname <- paste0("val_", ii)
        loading <- data.frame(temp_object$ALASCA$loading$time[temp_object$ALASCA$loading$time$PC %in% limPC_time, ], model = ii)
        # loading$covars <- object$numvariablelist[loading$covars]
        DBI::dbWriteTable(object$db.con, "time.loading", loading, append = T, overwrite = F)
        scores <- data.frame(temp_object$ALASCA$score$time[temp_object$ALASCA$score$time$PC %in% limPC_time, ], model = ii)
        # scores$time <- object$numtimelist[scores$time]
        # if(!object$separate_time_and_group){
        #  scores$group <- object$numgrouplist[scores$group]
        # }
        DBI::dbWriteTable(object$db.con, "time.score", scores, append = T, overwrite = F)
        if (nrow(get_covars(object)) > 0) DBI::dbWriteTable(object$db.con, "covars", data.frame(get_covars(temp_object), model = ii), append = T, overwrite = F)

        if (object$separate_time_and_group) {
          loading <- data.frame(temp_object$ALASCA$loading$group[temp_object$ALASCA$loading$group$PC %in% limPC_group, ], model = ii)
          # loading$covars <- object$numvariablelist[loading$covars]
          DBI::dbWriteTable(object$db.con, "group.loading", loading, append = T, overwrite = F)

          scores <- data.frame(temp_object$ALASCA$score$group[temp_object$ALASCA$score$group$PC %in% limPC_group, ], model = ii)
          # scores$time <- object$numtimelist[scores$time]
          # scores$group <- object$numgrouplist[scores$group]
          DBI::dbWriteTable(object$db.con, "group.score", scores, append = T, overwrite = F)
        }
        if (object$validate_regression) {
          DBI::dbWriteTable(object$db.con, "model_prediction", data.frame(temp_object$model_prediction, model = ii), append = T, overwrite = F)
        }

        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(object$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((object$n_validation_runs - ii) * time_all, 2), " seconds"))
        fname
      })
    } else {
      temp_object <- lapply(seq_len(object$n_validation_runs), FUN = function(ii) {
        log4r::info(object$log, paste0("- Run ", ii, " of ", object$n_validation_runs))
        start.time.this <- Sys.time()

        # Make resampled model
        temp_object <- prepare_validation_run(object, runN = ii)
        gc()

        # Rotate new loadings/scores to the original model
        if (object$optimize_score) {
          temp_object <- rotate_matrix_optimize_score(object = temp_object, target = object)
        } else {
          temp_object <- rotate_matrix(object = temp_object, target = object)
        }
        temp_object <- clean_alasca(temp_object)

        time_all <- difftime(Sys.time(), start.time.all, units = c("secs")) / ii
        log4r::info(object$log, paste0("--- Used ", round(difftime(Sys.time(), start.time.this, units = c("secs")), 2), " seconds. Est. time remaining: ", round((object$n_validation_runs - ii) * time_all, 2), " seconds"))
        temp_object
      })
    }
  }

  log4r::info(object$log, "Calculating percentiles for score and loading")
  object <- get_validation_percentiles(object, objectlist = temp_object)

  if (object$keep_validation_objects) {
    object$validation$temp_objects <- temp_object
  }
  return(object)
}

.procrustes <- function(loadings, target) {
  s <- t(loadings) %*% target
  w1 <- s %*% t(s)
  v1 <- t(s) %*% s
  w <- eigen(w1)$vectors
  ew <- diag(eigen(w1)$values)
  v <- eigen(v1)$vectors
  ev <- diag(eigen(v1)$values)
  o <- t(w) %*% s %*% v
  k <- diag(((diag(o)) / abs(diag(o))), nrow = nrow(o), ncol = nrow(o))
  ww <- w %*% k
  out <- list()
  out$t1 <- ww %*% t(v) # Rotation matrix
  out$procrust <- loadings %*% out$t1 # Rotated loadings
  return(out)
}

#' Rotate PCA
#'
#' This function rotates loadings and scores during validation
#'
#' Optimizes the rotation for lowest possible difference in score
#'
#' @param object ALASCA object to be rotated (and returned)
#' @param target ALASCA object acting as target
#' @return An ALASCA object
rotate_matrix_optimize_score <- function(object, target) {
  log4r::debug(object$log, "Starting rotation of time components")
  # We are only looking at components explaining more than a predefined value
  PCloading <- get_relevant_pcs(target, effect = "time")
  PCloading_t <- paste0("PC", PCloading)

  # PCA can give loadings with either sign, so we have to check whether switching signs improves the rotation
  N <- length(PCloading) # Number of components to look at
  
  # Create matrix with all possible combinations of signs
  vec <- c(-1, 1)
  lst <- lapply(numeric(N), function(x) vec)
  signMatrix <- as.matrix(expand.grid(lst))
  
  # Test all combinations and calculate residuals
  signVar <- vapply(seq_len(nrow(signMatrix) / 2), function(i) {
    c <- .procrustes(
      loadings = as.matrix(t(t(object$pca$loading$time[target$pca$loading$time, ..PCloading_t]) * signMatrix[i, ])),
      target = as.matrix(target$pca$loading$time[, ..PCloading_t])
    )
    sum((target$pca$score$time[, ..PCloading_t] - 
           as.matrix(t(t(object$pca$score$time[target$pca$score$time, ..PCloading_t]) * signMatrix[i, ])) %*% solve(c$t1))^2)
  }, FUN.VALUE = numeric(1))
  
  # Find the combination that minimizes the sum of squares
  minSignVar <- which(signVar == min(signVar))[1]
  
  # Switch signs
  for (i in PCloading){
    set(object$pca$loading$time, j = PCloading_t[i], value = object$pca$loading$time[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
    set(object$pca$score$time, j = PCloading_t[i], value = object$pca$score$time[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
  }

  # Final rotation
  c <- .procrustes(
    loadings = as.matrix(object$pca$loading$time[target$pca$loading$time, ..PCloading_t]),
    target = as.matrix(target$pca$loading$time[, ..PCloading_t])
  )

  object$pca$loading$time[target$pca$loading$time, (PCloading_t) := as.data.table(c$procrust)]
  object$pca$score$time[target$pca$score$time, (PCloading_t) := as.data.table(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]

  log4r::debug(object$log, "Completed rotation of time components")
  if (object$separate_time_and_group) {
    log4r::debug(object$log, "Starting rotation of group components")
    # We are only looking at components explaining more than a set limit
    PCloading <- get_relevant_pcs(target, effect = "group")
    PCloading_t <- paste0("PC", PCloading)

    # PCA can give loadings with either sign, so we have to check whether swithcing signs improves the rotation
    N <- length(PCloading)
    vec <- c(-1, 1)
    lst <- lapply(numeric(N), function(x) vec)
    signMatrix <- as.matrix(expand.grid(lst))
    signVar <- vapply(seq_len(nrow(signMatrix) / 2), function(i) {
      c <- .procrustes(
        loadings = as.matrix(t(t(object$pca$loading$group[target$pca$loading$group, ..PCloading_t]) * signMatrix[i, ])),
        target = as.matrix(target$pca$loading$group[, ..PCloading_t])
      )
      sum((target$pca$score$group[, ..PCloading_t] - 
             as.matrix(t(t(object$pca$score$group[target$pca$score$group, ..PCloading_t]) * signMatrix[i, ])) %*% solve(c$t1))^2)
    }, FUN.VALUE = numeric(1))
    
    minSignVar <- which(signVar == min(signVar))[1]
    for (i in PCloading){
      set(object$pca$loading$group, j = PCloading_t[i], value = object$pca$loading$group[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
      set(object$pca$score$group, j = PCloading_t[i], value = object$pca$score$group[, get(PCloading_t[i])] * signMatrix[minSignVar, i])
    }
    
    c <- .procrustes(
      loadings = as.matrix(object$pca$loading$group[target$pca$loading$group, ..PCloading_t]),
      target = as.matrix(target$pca$loading$group[, ..PCloading_t])
    )
    object$pca$loading$group[target$pca$loading$group, (PCloading_t) := as.data.table(c$procrust)]
    object$pca$score$group[target$pca$score$group, (PCloading_t) := as.data.table(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]
    log4r::debug(object$log, "Completed rotation of group components")
  }

  return(object)
}

#' Rotate PCA
#'
#' This function rotates loadings and scores during validation
#'
#' @param object ALASCA object to be rotated (and returned)
#' @param target ALASCA object acting as target
#' @return An ALASCA object
rotate_matrix <- function(object, target) {
  PCloading <- get_relevant_pcs(target, effect = "time")
  PCloading_t <- paste0("PC", PCloading)
  c <- .procrustes(
    loadings = as.matrix(object$pca$loading$time[target$pca$loading$time, ..PCloading_t]),
    target = as.matrix(target$pca$loading$time[, ..PCloading_t])
  )

  object$pca$loading$time[target$pca$loading$time, (PCloading_t) := as.data.frame(c$procrust)]
  object$pca$score$time[target$pca$score$time, (PCloading_t) := as.data.frame(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]

  if (object$separate_time_and_group) {
    PCloading <- get_relevant_pcs(target, effect = "group")
    PCloading_t <- paste0("PC", PCloading)
    c <- .procrustes(
      loadings = as.matrix(object$pca$loading$group[target$pca$loading$group, ..PCloading_t]),
      target = as.matrix(target$pca$loading$group[, ..PCloading_t])
    )

    object$pca$loading$group[target$pca$loading$group, (PCloading_t) := as.data.frame(c$procrust)]
    object$pca$score$group[target$pca$score$group, (PCloading_t) := as.data.frame(as.matrix(.SD) %*% solve(c$t1)), .SDcols = PCloading_t]
  }

  return(object)
}

#' Extract percentiles
#'
#' This function extract percentiles during validation
#'
#' @param object ALASCA object
#' @param objectlist List of ALASCA objects
#' @return An ALASCA object
get_validation_percentiles <- function(object, objectlist) {
  if ("low" %in% colnames(object$ALASCA$loading$time)) {
    object$ALASCA$loading$time$low <- NULL
    object$ALASCA$loading$time$high <- NULL
    object$ALASCA$score$time$low <- NULL
    object$ALASCA$score$time$high <- NULL
    if (object$separate_time_and_group) {
      object$ALASCA$loading$group$low <- NULL
      object$ALASCA$loading$group$high <- NULL
      object$ALASCA$score$group$low <- NULL
      object$ALASCA$score$group$high <- NULL
    }
  }
  log4r::debug(object$log, message = "Starting get_validation_percentiles_loading")
  object <- get_validation_percentiles_loading(object, objectlist)
  log4r::debug(object$log, message = "Completed get_validation_percentiles_loading")
  log4r::debug(object$log, message = "Starting get_validation_percentiles_score")
  object <- get_validation_percentiles_score(object, objectlist)
  log4r::debug(object$log, message = "Completed get_validation_percentiles_score")
  if (object$validate_regression) {
    log4r::debug(object$log, message = "Starting get_validation_percentiles_regression")
    object <- get_validation_percentiles_regression(object, objectlist)
    log4r::debug(object$log, message = "Completed get_validation_percentiles_regression")
  }
  if (nrow(get_covars(object)) > 0) {
    log4r::debug(object$log, message = "Starting get_validation_percentiles_covars")
    object <- get_validation_percentiles_covars(object, objectlist)
    log4r::debug(object$log, message = "Completed get_validation_percentiles_covars")
  }

  return(object)
}

#' Extract percentiles for regressions
#'
#' This function extract percentiles for validation of regression
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_regression <- function(object, objectlist) {
  if (object$save_to_disk) {
    res <- DBI::dbSendQuery(object$db.con, "SELECT * FROM 'model_prediction'")
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df <- rbindlist(lapply(objectlist, function(x) x$model_prediction))
  }
  df <- df[, as.list(quantile(pred, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(group, time, variable)]
  colnames(df) <- c("group", "time", "variable", "low", "high")

  object$model_prediction <- merge(object$model_prediction, df)
  return(object)
}

#' Extract percentiles for covariates
#'
#' This function extract percentiles for validation of covariates
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_covars <- function(object, objectlist) {
  if (object$save_to_disk) {
    res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'covars'"))
    df <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df <- rbindlist(lapply(objectlist, function(x) get_covars(x)))
  }
  df <- df[, as.list(quantile(estimate, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(covar, variable)]
  colnames(df) <- c("covar", "variable", "low", "high")

  object$covar_coefficients <- merge(object$covar_coefficients, df)
  return(object)
}

#' Extract percentiles for loading
#'
#' This function extract percentiles during validation of loadings
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_loading <- function(object, objectlist) {
  PC_time <- get_relevant_pcs(object, effect = "time")
  if (object$save_to_disk) {
    res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.loading' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
    df_time <- setDT(DBI::dbFetch(res))
    DBI::dbClearResult(res)
  } else {
    df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$loading$time[x$ALASCA$loading$time$PC %in% PC_time, ]), fill = TRUE)
  }

  object$validation$time$loading <- df_time[, as.list(quantile(loading, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(PC, covars)]
  colnames(object$validation$time$loading) <- c("PC", "covars", "low", "high")
  object$ALASCA$loading$time <- merge(object$ALASCA$loading$time, object$validation$time$loading, all.x = TRUE)

  if (object$separate_time_and_group) {
    PC_group <- get_relevant_pcs(object, effect = "group")
    if (object$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'group.loading' WHERE PC IN(", paste(PC_group, collapse = ", "), ")"))
      df_group <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_group <- rbindlist(lapply(objectlist, function(x) x$ALASCA$loading$group[x$ALASCA$loading$group$PC %in% PC_group, ]), fill = TRUE)
    }

    object$validation$group$loading <- df_group[, as.list(quantile(loading, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(PC, covars)]
    colnames(object$validation$group$loading) <- c("PC", "covars", "low", "high")
    object$ALASCA$loading$group <- merge(object$ALASCA$loading$group, object$validation$group$loading, all.x = TRUE)
  }
  return(object)
}

#' Extract percentiles for score
#'
#' This function extract percentiles during validation of scores
#'
#' @inheritParams get_validation_percentiles
#' @return An ALASCA object
get_validation_percentiles_score <- function(object, objectlist) {
  if (object$separate_time_and_group) {
    # Separate time and group effects

    PC_time <- get_relevant_pcs(object, effect = "time")
    if (object$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$time[x$ALASCA$score$time$PC %in% PC_time, ]), fill = TRUE)
    }

    object$validation$time$score <- df_time[, as.list(quantile(score, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(PC, time)]
    colnames(object$validation$time$score) <- c("PC", "time", "low", "high")
    object$ALASCA$score$time <- merge(object$ALASCA$score$time, object$validation$time$score, all.x = TRUE)
    object$ALASCA$score$time[, time := factor(time, levels = object$timelist), ]
    object$ALASCA$score$time[, group := factor(group, levels = object$grouplist), ]

    PC_group <- get_relevant_pcs(object, effect = "group")
    if (object$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'group.score' WHERE PC IN(", paste(PC_group, collapse = ", "), ")"))
      df_group <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_group <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$group[x$ALASCA$score$group$PC %in% PC_group, ]), fill = TRUE)
    }

    object$validation$group$score <- df_group[, as.list(quantile(score, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(PC, time, group)]
    colnames(object$validation$group$score) <- c("PC", "time", "group", "low", "high")
    object$ALASCA$score$group <- merge(object$ALASCA$score$group, object$validation$group$score, all.x = TRUE)
    object$ALASCA$score$group[, time := factor(time, levels = object$timelist), ]
    object$ALASCA$score$group[, group := factor(group, levels = object$grouplist), ]
  } else {
    # Pooled time and groups effects
    PC_time <- get_relevant_pcs(object, effect = "time")
    if (object$save_to_disk) {
      res <- DBI::dbSendQuery(object$db.con, paste0("SELECT * FROM 'time.score' WHERE PC IN(", paste(PC_time, collapse = ", "), ")"))
      df_time <- setDT(DBI::dbFetch(res))
      DBI::dbClearResult(res)
    } else {
      df_time <- rbindlist(lapply(objectlist, function(x) x$ALASCA$score$time[x$ALASCA$score$time$PC %in% PC_time, ]), fill = TRUE)
    }

    object$validation$time$score <- df_time[, as.list(quantile(score, probs = object$limitsCI, type = object$validation_quantile_method)), by = .(PC, time, group)]
    colnames(object$validation$time$score) <- c("PC", "time", "group", "low", "high")
    object$ALASCA$score$time <- merge(object$ALASCA$score$time, object$validation$time$score, all.x = TRUE)
    object$ALASCA$score$time[, time := factor(time, levels = object$timelist), ]
    object$ALASCA$score$time[, group := factor(group, levels = object$grouplist), ]
  }
  return(object)
}

#' Get relevant PCs
#'
#' This function extract percentiles during validation
#'
#' @param x Explanatory power of PC
#' @return A vector with relevant PCs
get_relevant_pcs <- function(object, effect = "time") {
  if (effect == "time") {
    PC <- object$ALASCA$loading$explained$time >= object$explanatorylimit
  } else {
    PC <- object$ALASCA$loading$explained$group >= object$explanatorylimit
  }
  PC[1:2] <- TRUE
  return(which(PC))
}

#' Validate underlying regression models
#'
#' This function calcuates predictions from each regression model
#'
#' @param object An ALASCA object
#' @return An ALASCA object
get_regression_predictions <- function(object) {
  if (!object$minimize_object) {
    # This is not a validation run
    log4r::info(object$log, message = "Calculating predictions from regression models")
  }
  regCoeffAll <- dcast(data = object[["regression_coefficients"]], variable~covar, value.var = "estimate")
  rownames(regCoeffAll) <- regCoeffAll$variable
  regModel <- unique(model.matrix(object$formula, data = object$df))
  if (object$equal_baseline) {
    regModel <- regModel[, !grepl(paste0("time", levels(object$df$time)[1]), colnames(regModel), fixed = TRUE), drop = FALSE]
  }
  regModel <- regModel[, !grepl("|", colnames(regModel), fixed = TRUE)]
  if (object$keep_terms != "") {
    regModel <- regModel[, grepl(paste0(c("time", "group", object$keep_terms), collapse = "|"), colnames(regModel)), drop = FALSE]
  } else {
    regModel <- regModel[, grepl(paste0(c("time", "group"), collapse = "|"), colnames(regModel)), drop = FALSE]
  }
  regModel <- unique(regModel)
  
  if (object$keep_terms != "") {
    object$model_prediction <- melt(
      cbind(
        as.matrix(regModel) %*% as.matrix(regCoeffAll[colnames(regModel), -1]),
        object$df[as.numeric(rownames(regModel)), .SD, .SDcols = c("time", "group", object$keep_terms)]
      ),
      id.vars = c("time", "group", object$keep_terms), value.name = "pred"
    )
  } else {
    object$model_prediction <- melt(
      cbind(
        as.matrix(regModel) %*% as.matrix(regCoeffAll[colnames(regModel), -1]),
        object$df[as.numeric(rownames(regModel)), .SD, .SDcols = c("time", "group")]
        ),
      id.vars = c("time", "group"), value.name = "pred"
      )
  }
  
  if (object$keep_terms != "") {
    object$model_prediction[, group := apply(.SD, 1, paste, collapse = " - "), .SDcols = c("group", object$keep_terms)]
    object$model_prediction[, group := factor(group, levels = object$grouplist)]
  } 

  if (!object$minimize_object) {
    # This is not a validation run
    log4r::info(object$log, message = "Finished calculating predictions from regression models!")
  }
  return(object)
}

#' Make a single validation run
#'
#' This function ...
#'
#' @param object An ALASCA object
#' @return An ALASCA object
prepare_validation_run <- function(object, runN = NA) {
  if (object$validation_method %in% c("loo", "jack-knife", "jackknife")) {
    # Use leave-one-out validation
    selectedParts <- data.frame()

    if (object$method %in% c("LMM")) {
      if (any(is.na(object$validation_ids))) {
        # For each group, divide the participants into n_validation_folds groups, and select n_validation_folds-1 of them
        selectedParts <- lapply(unique(object$stratification_vector), function(gr) {
          selectedParts_temp_all <- unique(object$df[object$stratification_vector == gr, ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$n_validation_folds
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% unlist(selectedParts)
        )
      } else {
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% object$validation_ids[runN, ]
        )
      }
    } else if (object$method %in% c("LM")) {
      object$df$ID <- c(seq_len(nrow(object$df)))
      if (any(is.na(object$validation_ids))) {
        # For each group, divide the participants into n_validation_folds groups, and select n_validation_folds-1 of them
        selectedParts <- lapply(unique(object$stratification_vector), function(gr) {
          selectedParts_temp_all <- unique(object$df[object$stratification_vector == gr, ID])
          selectedParts_temp_ticket <- seq_along(selectedParts_temp_all) %% object$n_validation_folds
          selectedParts_temp_ticket <- selectedParts_temp_ticket[sample(seq_along(selectedParts_temp_ticket), length(selectedParts_temp_ticket))]
          selectedParts_temp_all[selectedParts_temp_ticket != 1]
        })

        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% unlist(selectedParts)
        )
      } else {
        temp_object <- ALASCA(
          validation_object = object,
          validation_participants = object$df_raw[, ID] %in% object$validation_ids[runN, ]
        )
      }
    }
  } else if (object$validation_method == "bootstrap") {
    # Use bootstrap validation
    # When using bootstrapping, we resample participants with replacement
    
    log4r::debug(object$log, message = "Completed preparation of bootstrap sample")
    participants_in_bootstrap <- get_bootstrap_ids(object, runN = runN)
    
    # Create bootstrap object without data
    bootobject <- object[!names(object) %in% c("df", "df_raw")]
    bootobject <- get_bootstrap_data(bootobject, df_raw = object$df_raw, participants_in_bootstrap)
    
    if (object$validation_assign_new_ids) {
      bootobject$df_raw[, ID := uniqueIDforBootstrap] # Replace ID
    }
    
    if (object$save_validation_ids) {
      write(paste0(participants_in_bootstrap$old_id, collapse = ";"),
            file = get_filename(object = object, prefix = "bootstrapID_", filetype = ".csv", overwrite = TRUE), append = TRUE
      )
    }
    
    log4r::debug(object$log, message = "Completed preparation of bootstrap sample")
    
    temp_object <- ALASCA(
      validation_object = bootobject,
      validation_participants = rep(TRUE, nrow(bootobject$df_raw))
    )
    temp_object$validation$original_ids <- participants_in_bootstrap$old_id
  } 

  return(temp_object)
}

#' Make bootstrap sample
#'
#' Get data frame with new and old IDs
#'
#' @param object An ALASCA object
#' @return A data frame
get_bootstrap_ids <- function(object, runN = NA) {
  if (is.na(runN)) {
    participants_in_bootstrap <- data.frame(
      new_id = seq(object$df[, uniqueN(ID)]),
      old_id = rbindlist(
        lapply(unique(object$stratification_vector), function(stratification_group){
          list(
            old_id = sample(
              unique(object$df[object$stratification_vector == stratification_group, ID]),
              replace = TRUE)
          )
        })
      )
    )
  } else {
    participants_in_bootstrap <- data.frame(
      new_id = seq(object$validation_ids[runN, ]),
      old_id = matrix(object$validation_ids[runN, ])
    )
  }
  return(participants_in_bootstrap)
}

#' Make bootstrap data set
#'
#' Get data frame
#'
#' @param object An ALASCA object
#' @return A data frame
get_bootstrap_data <- function(bootobject, df_raw, participants_in_bootstrap) {
  selected_rows <- rbindlist(
    lapply(participants_in_bootstrap$new_id, function(participant){
      list(
        new_id = participant,
        row_nr = df_raw[, .I[ID == participants_in_bootstrap$old_id[participant]]]
      )
    })
  )
  bootobject[["df_raw"]] <- df_raw[selected_rows$row_nr]
  bootobject[["df_raw"]][, uniqueIDforBootstrap := selected_rows$new_id]
  bootobject[["df_raw"]][, originalIDbeforeBootstrap := ID]
  bootobject$modmat <- bootobject$modmat[selected_rows$row_nr,]
  return(bootobject)
}
