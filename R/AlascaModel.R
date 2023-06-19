#' R6 Class Representing an ALASCA model
#'
#' @description
#' The object builds the ALASCA model and contains the data
AlascaModel <- R6::R6Class("AlascaModel",
  lock_objects = FALSE,
  public = list(
    #' @field df Data table/frame. The data to analyze
    df = NULL,
    #' @field formula An AlascaForula object
    formula = NULL,
    #' @field wide Boolean. Whether the provided data is in wide format
    wide = FALSE,
    #' @field scale_function How to scale the data. Options are `NULL`, custom function, or `"sdall"`, `"sdref"`, `"sdt1"`, `"sdreft1"`
    scale_function = "sdall",

    #' @field ignore_missing If TRUE, ignore missing predictive values
    ignore_missing = FALSE,
    #' @field ignore_missing_covars If TRUE, ignore missing covariate values
    ignore_missing_covars = FALSE,
    #' @field version Version number
    version = "1.0.11",
    #' @field update_date Date of latest update
    update_date = "2023-06-19",

    # Effect matrices
    #' @field separate_effects If TRUE, try to separate the effects
    separate_effects = FALSE,
    #' @field equal_baseline If TRUE, remove interaction between baselines
    equal_baseline = FALSE,
    #' @field effect_list List. Contains info related to the effects
    effect_list = list(
      expr = NULL,
      terms = NULL,
      model_matrix = NULL,
      effect_matrix = NULL,
      pca = NULL,
      saved_scores = NULL,
      saved_loadings = NULL
    ),

    # Validation settings
    #' @field n_validation_folds Integer. If using jack-knife validation, exclude 1/n_validation_folds of the participants at each iteration
    n_validation_folds = 7L,
    #' @field n_validation_runs Integer. Number of iterations to use for validation
    n_validation_runs = 1000L,
    #' @field validation_quantile_method Integer between 1 and 9. See [stats::quantile()] for details
    validation_quantile_method = 2L,
    #' @field save_validation_ids If TRUE, save the participants in each validation iteration to a csv file
    save_validation_ids = FALSE,
    #' @field optimize_score If TRUE, test all combinations of signs for the most important loadings, and choose the combination being the best fit
    optimize_score = TRUE,
    #' @field optimize_PCs If TRUE, test if principal components have to be re-ordered during validation
    optimize_PCs = FALSE,
    #' @field validate If TRUE, validate the model
    validate = FALSE,
    #' @field validate_regression If TRUE, validate get marginal means
    validate_regression = TRUE,
    #' @field validation Boolean. Synonym to `validate`
    validation = FALSE,
    #' @field validation_method String. Defines the validation method; `"bootstrap"` (default) or `"jack-knife"`
    validation_method = "bootstrap",
    #' @field validation_ids A data frame where each row contains the ids for one validation iteration
    validation_ids = NULL,
    validation_object = NA,
    #' @field validation_assign_new_ids Logical. Assign new IDs during validation. Must be `TRUE` for reduce_dimensions to work
    validation_assign_new_ids = TRUE,
    #' @field limitsCI Lower and upper quantile to use for validation
    limitsCI = c(0.025, 0.975),
    #' @field compress_validation Integer between 0 and 100. See [fst::write_fst()] for details
    compress_validation = 80,
    permutation_within_participants = NULL,
    permutation_across_participants = NULL,

    # Reduce dimensions
    #' @field reduce_dimensions Boolean. Use PCA to reduce data dimensions prior to analysis
    reduce_dimensions = FALSE,
    #' @field pca_function String or custom function. Which pca function to use for dimension reduction, either "prcomp" or "irlba" or "princomp" or custom function
    pca_function = "prcomp",
    reduced_df = list(
      explanatory_power = NULL,
      nComps = NULL,
      limit = 0.95,
      loading = NULL,
      score = NULL,
      df = NULL,
      variables = NULL
    ),

    # Save to disk
    #' @field save_to_disk Write model data to disk to reduce memory usage
    save_to_disk = FALSE,
    #' @field db_method String. Use a `"duckdb"` or a `"SQLite"` database for validation data
    db_method = "duckdb", # SQLite

    # Save
    #' @field filename Filename for the saved model
    filename = "ALASCA",
    #' @field filepath Where to save the model. Defaults to `ALASCA/<timestamp>`
    filepath = NULL,
    #' @field save Save model data and plots
    save = FALSE,

    #' @field method String. Can be `"LM"` or `"LMM"`
    method = NULL,
    #' @field max_PC Integer. The maximal number of principal components to keep for further analysis
    max_PC = 20L,
    #' @field use_Rfast Boolean. If `TRUE` (default), use Rfast, else use lm or lme4
    use_Rfast = TRUE,
    #' @field p_adjust_method String. See [stats::p.adjust()]
    p_adjust_method = NULL,
    #' @field participant_column String. The column used for IDs. If not provided, it will guess based on random effect or `ID`
    participant_column = NULL,
    scale_function.center = FALSE,
    #' @field stratification_column String. Name of the column to use for stratification during validation
    stratification_column = NULL,
    stratification_vector = NA,
    minimize_object = FALSE,
    #' @field explanatory_limit Only validate components explaining more than `explanatory_limit` of the variance
    explanatory_limit = 0.05,
    #' @field init_time The time when the object is initialized
    init_time = NULL,

    # Log settings
    #' @field log_to String deciding logging target: `"all"` (default), `"file"`, `"console"`, `"none"`
    log_to = "all",
    log_file = NULL,
    logger = NULL,
    #' @field log_level String. What level of log messages to print; `"DEBUG"`, `"INFO"`, `"WARN"`, `"ERROR"`
    log_level = NULL,
    uselog = TRUE,
    #' @field do_debug Boolean. Log more details
    do_debug = FALSE,
    covars = NULL,
    new_formula = NULL,
    my_df_rows = NULL,
    modmat = NULL,
    cnames_modmat = NULL,
    covar_coefficients = NULL,
    #' @field finished Boolean. Indicates whether the model has been successfully initiated
    finished = FALSE,

    #' @field ALASCA List. Contains all model outputs: `score`, `loading`, `explained` and `significant_PCs`
    ALASCA = list(
      score = NULL,
      loading = NULL,
      explained = NULL,
      significant_PCs = NULL
    ),
    initialize = function(df, formula, effects, ...) {

      # Fetch user inputs
      inputs <- list(...)
      self$effect_list$expr <- effects
      self$splot <- AlascaPlot$new(model = self)

      for (i in seq_along(inputs)) {
        if (substr(names(inputs)[i], 1, 5) == "plot.") {
          self$splot[[gsub("plot.", "", names(inputs)[i], fixed = TRUE)]] <- inputs[[i]]
        } else if (substr(names(inputs)[i], 1, 18) == "reduce_dimensions.") {
          self$reduced_df[[gsub("reduce_dimensions.", "", names(inputs)[i], fixed = TRUE)]] <- inputs[[i]]
        } else {
          self[[names(inputs)[i]]] <- inputs[[i]]
        }
      }
      if(!is.data.table(self$df)) {
        self$df <- setDT(self$df)
      }

      self$init_time <- Sys.time()
      if (is.null(self$filepath)) {
        self$filepath <- paste0("ALASCA/", strftime(self$init_time, format = "%Y-%m-%dT%H%M%S"), "/")
      } else {
        if (substr(self$filepath, nchar(self$filepath), nchar(self$filepath)) != "/") {
          self$filepath <- paste0(self$filepath, "/")
        }
      }

      dir.create(paste0(self$filepath, "plot/"), recursive = TRUE)
      self$log_file <- paste0(self$filepath, "ALASCA.log")

      # Start logging
      self$log_level <- ifelse(self$do_debug, "DEBUG", "INFO")
      if (self$log_to == "file") {
        self$logger <- log4r::logger(self$log_level, appenders = list(log4r::file_appender(self$log_file)))
      } else if (self$log_to == "console") {
        self$logger <- log4r::logger(self$log_level, appenders = list(log4r::console_appender()))
      } else if (self$log_to == "none") {
        self$uselog <- FALSE
      } else {
        self$logger <- log4r::logger(self$log_level, appenders = list(log4r::console_appender(), log4r::file_appender(self$log_file)))
      }
      self$log(paste("Initializing", self$print_version()))

      # Process the formula
      self$formula <- AlascaFormula$new(formula, model = self)
      self$set_effect_terms()

      # Keep a copy of unscaled data
      self$df_raw <- AlascaDataset$new(data_df = df, model = self)
      self$formula$get_regression_formula()
      self$my_df_rows <- self$df_raw$rows_to_serve

      # self$stratification_vector <- self$df_raw$data_df[, get(self$stratification_column)]

      # Scale data
      self$get_scaling_function()
      self$df <- self$scale_function(self$df_raw$df)

      self$get_pca_function()
      self$splot$group <- self$get_plot_group

      self$validate <- self$validate || self$validation
      self$validation <- NULL
      self$n_validation_runs <- ifelse(is.null(self$validation_ids),
        self$n_validation_runs,
        nrow(self$validation_ids)
      )
      self$filepath <- ifelse(is.na(self$filepath), NA, ifelse(substr(self$filepath, nchar(self$filepath), nchar(self$filepath)) == "/", self$filepath, paste0(self$filepath, "/")))
      self$save_to_disk <- self$validate && self$save_to_disk
      if (self$save_to_disk) {
        if (self$db_method == "SQLite") {
          self$db_driver <- RSQLite::dbDriver("SQLite")
          self$db_filename <- paste0(self$filepath, "ALASCA.db")
          self$db_con <- DBI::dbConnect(self$db_driver, dbname = self$db_filename)
        } else {
          self$db_filename <- paste0(self$filepath, "ALASCA.duckdb")
          self$db_con <- DBI::dbConnect(duckdb::duckdb(), dbdir = self$db_filename, read_only = FALSE)
        }

        DBI::dbWriteTable(
          self$db_con, "variables",
          data.table(
            id = seq_along(self$get_levels("variable", reduced = FALSE)),
            variable = self$get_levels("variable", reduced = FALSE)
          )
        )
      }
      
      if (self$validate && self$validation_method == "permutation") {
        if (is.null(self$permutation_within_participants)) {
          self$permutation_within_participants <- self$effect_terms[1]
        }
        if (is.null(self$permutation_across_participants)) {
          self$permutation_across_participants <- self$effect_terms[-1]
        }
      }

      # self$ALASCA.version <- print_version(get = "version")
      # self$ALASCA.version.date <- print_version(get = "date")

      # Clean input ----
      self$log("Has initialized the ALASCA model. Next step is to clean it and check input", level = "DEBUG")

      # Build the ALASCA model ----
      self$build_model()
    },
    finalize = function() {
      if (self$save_to_disk && !self$minimize_object) {
        if (self$db_method == "SQLite") {
          DBI::dbDisconnect(self$db_con)
        } else {
          DBI::dbDisconnect(self$db_con, shutdown = TRUE)
        }
      }
    },
    #' @description
    #' Update the current model (used for validation)
    update = function() {

      ## Avoid recursion
      self$validate <- FALSE

      ## Save space
      self$minimize_object <- TRUE

      # Build the ALASCA model ----

      self$build_model()

      # To save space, we remove unnecessary embedded data ----
      self$log("Starting to remove embedded data", level = "DEBUG")
      self$remove_embedded_data()
      self$log("Completed remove embedded data", level = "DEBUG")

      # invisible(self)
    },
    #' @description
    #' Function for logging messages using the log4r package
    #' @param message The message to log
    #' @param level Level of the message; DEBUG, INFO, WARN, ERROR, FATAL
    log = function(message, level = "INFO") {
      if (self$uselog) {
        log4r::levellog(self$logger, level = level, message)
      }
    },
    #' @description
    #' Main function for plots
    #' @param effect Integer or vector. Which(s) effect(s) to plot
    #' @param component Integer or vector. Which(s) component(s) to plot
    #' @inheritParams AlascaPlot
    plot = function(effect = 1, component = 1, ...) {
      self$splot$effect_i <- effect
      self$splot$component <- component
      self$splot$call_plot(...)
    },
    #' @description
    #' Resample participants for validation
    get_validation_ids = get_validation_ids,
    get_residuals = function(variable = NULL) {
      if (self$use_Rfast) {
        self$log("Residuals are not calculated when using Rfast. Re-build the model with `use_Rfast = TRUE`", level = "ERROR")
        stop()
      }
      if (is.null(variable)) {
        variable <- self$get_levels("variable")
      }
      rbindlist(lapply(variable, function(x) {
        list(
          variable = x,
          residuals = residuals(self$regression_model[[x]])
        )
      }))
    },
    set_effect_terms = function() {
      self$effect_list$terms <- lapply(self$effect_list$expr, function(x) {
        x <- gsub(" ", "", x)
        x <- unlist(strsplit(x, split = "\\:|\\+|\\||\\*"))
        unique(x)
      })
      if (is.null(self$splot$x_label)) {
        # If the x_label has not been provided, use the first regression term
        self$splot$x_label <- self$splot$capitalize(self$effect_list$terms[[1]][[1]])
      }
    },
    set_design_matrices = function() {
      #' TODO: There is a bug when the order of the terms in `self$effect_list$expr` differs from `t1:t2:t3` in  a three-way interaction
      if (is.null(self$effect_list$model_matrix)) {
        self$effect_list$model_matrix <- lapply(self$effect_list$expr, function(x) {
          mm <- model.matrix(as.formula(paste0("value ~ ", x)), data = self$df)
          mm <- mm[, colnames(mm) %in% self$cnames_modmat]
          if (ncol(mm) > 2) {
            mm[, -1]
          } else {
            mm
          }
        })
      }

      return(
        lapply(self$effect_list$model_matrix, function(mm) {
          # log4r::debug(self$log, paste("mm_dim: ", dim(mm)))
          # log4r::debug(self$log, paste("mm_dim: ", dim(self$df)))
          mm
        })
      )
    },
    #' @description
    #' Switch the sign of scores and loadings
    #' @param effect_i The effect to reflect, `0` or `NULL` reflects the entire model
    #' @param component The component to reflect, `0` or `NULL` reflects the entire model
    flip = function(effect_i = 0, component = 0) {
      if (effect_i[[1]] == 0 || is.null(effect_i)) {
        effect_i <- seq_along(self$effect_list$expr)
      }
      if (component[[1]] == 0 || is.null(component)) {
        component <- unique(self$ALASCA$loading[[1]]$PC)
      }

      for (effect_k in effect_i) {
        self$ALASCA$loading[[effect_k]][PC %in% component, loading := -loading]
        self$ALASCA$score[[effect_k]][PC %in% component, score := -score]

        if (self$validate) {
          self$ALASCA$loading[[effect_k]][PC %in% component, high := -high]
          self$ALASCA$loading[[effect_k]][PC %in% component, low := -low]

          self$ALASCA$score[[effect_k]][PC %in% component, high := -high]
          self$ALASCA$score[[effect_k]][PC %in% component, low := -low]

          if (self$save_to_disk) {
            DBI::dbSendQuery(
              self$db_con,
              paste0("UPDATE loading
                          SET loading = -loading,
                          high = -high,
                          low = -low
                          WHERE effect = ", effect_k, " AND PC IN(", paste(component, collapse = ", "), ")")
            )
            DBI::dbSendQuery(
              self$db_con,
              paste0("UPDATE score
                                SET score = -score,
                                high = -high,
                                low = -low
                                WHERE effect = ", effect_k, " AND PC IN(", paste(component, collapse = ", "), ")")
            )
          } else {
            tmp <- fst::read_fst(self$effect_list$saved_scores[[effect_k]], as.data.table = TRUE)
            tmp[PC %in% component, score := -score]
            fst::write_fst(tmp, path = self$effect_list$saved_scores[[effect_k]])

            tmp <- fst::read_fst(self$effect_list$saved_loadings[[effect_k]], as.data.table = TRUE)
            tmp[PC %in% component, loading := -loading]
            fst::write_fst(tmp, path = self$effect_list$saved_loadings[[effect_k]])
          }
        }
      }
    },
    #' @description
    #' Returns the reference level of a given column
    #' @param columns A column containing factors
    #' @return The reference level
    get_ref = function(columns) {
      vapply(columns, FUN = function(column) {
        if (self$reduce_dimensions && column == "variable" && !self$finished) {
          self$reduced_df[["variables"]][[1]]
        } else {
          self$df_raw$level_list[[column]][[1]]
        }
      }, FUN.VALUE = character(1))
    },
    #' @description
    #' Returns all the levels of a given column
    #' @param column A column containing factors
    #' @param reduced Boolean. If `TRUE`, returns the PCs instead of actual variables if dimensions are reduced with PCA
    #' @return A vector with factor levels
    get_levels = function(column, reduced = TRUE) {
      if (self$reduce_dimensions && reduced && column == "variable" && !self$finished) {
        self$reduced_df[["variables"]]
      } else {
        self$df_raw$level_list[[column]]
      }
    },
    #' @description
    #' Returns the most interesting principal components (i.e., components explaining more than a given limit of variance: `explanatory_limit`)
    #' @param x Index corresponding to the effect of interest
    #' @return A vector with principal components
    get_PCs = function(x) self$ALASCA$significant_PCs[[x]],
    get_predictions = function() {
      
      if (!self$get_plot_group %in% colnames(self$model_prediction)) {
        self$model_prediction[, (self$splot$group):= do.call(paste, c(.SD, sep = "-")), .SDcols = self$effect_terms[-1]]
      }
      
      if (self$equal_baseline) {
        # Must add baselines for missing interactions (baseline x group)
        baseline_to_add <- self$model_prediction[get(self$effect_terms) %in% self$get_ref(self$effect_terms)]
        
        if (length(self$effect_terms) == 2) {
          # In this case, there are no main effects except time
          baseline_to_add <- rbindlist(lapply(self$get_levels(self$get_plot_group)[-1], function(gr) {
            baseline_to_add[, (self$get_plot_group) := gr]
          }))
        } else if (length(self$effect_terms) == 3) {
          # In this case, there may be additional main effects to take into account
          terms_in_formula <- attr(terms(self$formula$formula_wo_random), "term.labels")
          
          # Find main terms
          terms_in_formula <- terms_in_formula[!grepl(pattern = ":", terms_in_formula, fixed = TRUE)]
          
          # Remove x_data
          terms_in_formula <- terms_in_formula[-1]
          
          # These are the groups missing
          terms_to_add <- self$get_levels(self$get_plot_group)
          terms_to_add <- terms_to_add[!terms_to_add %in% baseline_to_add[, get(self$get_plot_group)]]
          match_groups <- self$model_prediction[, unique(.SD), .SDcols = c(self$effect_terms[-1], self$get_plot_group)]
          
          terms_to_update <- self$effect_terms[-1]
          terms_to_update <- terms_to_update[!terms_to_update %in% terms_in_formula]
          
          baseline_to_add <- rbindlist(lapply(terms_to_add, function(gr) {
            baseline_to_add[get(terms_in_formula) == match_groups[get(self$get_plot_group) == gr, get(terms_in_formula)],
                            (self$get_plot_group) := gr]
            
          }))
          
          baseline_to_add <- unique(baseline_to_add)
          
          for (gr in terms_to_add) {
            rf <- match_groups[get(self$get_plot_group) == gr, get(terms_to_update)]
            baseline_to_add[get(self$get_plot_group) == gr, (terms_to_update) := rf]
          }
          
        } else {
          self$model$log("Sorry, this is not supported at the moment", level = "ERROR")
          stop()
        }
        

        predictions <- rbind(self$model_prediction, baseline_to_add)
      } else {
        predictions <- self$model_prediction
      }
      return(predictions)
    },
    #' @description
    #' Return scores
    #' @param effect_i The effect to reflect, `0` or `NULL` reflects the entire model
    #' @param component The component to reflect, `0` or `NULL` reflects the entire model
    get_scores = function(effect_i = 1, component = 1) {
      if (length(effect_i) == 1) {
        if (effect_i < 1) {
          effect_i <- length(self$ALASCA$score)
          lapply(effect_i, function(i) self$ALASCA$score[[i]][PC %in% component])
        } else {
          self$ALASCA$score[[effect_i]][PC %in% component]
        }
      } else {
        lapply(effect_i, function(i) self$ALASCA$score[[i]][PC %in% component])
      }
    },
    #' @description
    #' Return loadings
    #' @param effect_i The effect to reflect, `0` or `NULL` reflects the entire model
    #' @param component The component to reflect, `0` or `NULL` reflects the entire model
    #' @param n_limit Return only two times this number of loadings (the `n_limit` highest and lowest loadings). Use `0` to return all (default)
    get_loadings = function(effect_i = 1, component = 1, n_limit = 0) {
      if (n_limit > 0) {
        # Fetch top and bottom for each requested effect and component
        if (length(effect_i) == 1) {
          if (effect_i < 1) {
            effect_i <- length(self$ALASCA$loading)
            lapply(effect_i, function(i) {
              lapply(component, function(j) {
                self$ALASCA$loading[[i]][PC == j][c(
                  Rfast::nth(loading, k = n_limit, num.of.nths = n_limit, descending = FALSE, index.return = TRUE),
                  Rfast::nth(loading, k = n_limit, num.of.nths = n_limit, descending = TRUE, index.return = TRUE)
                )]
              })
            })
          } else {
            lapply(component, function(j) {
              self$ALASCA$loading[[effect_i]][PC == j][c(
                Rfast::nth(loading, k = n_limit, num.of.nths = n_limit, descending = FALSE, index.return = TRUE),
                Rfast::nth(loading, k = n_limit, num.of.nths = n_limit, descending = TRUE, index.return = TRUE)
              )]
            })
          }
        } else {
          lapply(effect_i, function(i) {
            lapply(component, function(j) {
              self$ALASCA$loading[[i]][PC == j][c(
                Rfast::nth(loading, k = n_limit, num.of.nths = n_limit, descending = FALSE, index.return = TRUE),
                Rfast::nth(loading, k = n_limit, num.of.nths = n_limit, descending = TRUE, index.return = TRUE)
              )]
            })
          })
        }
      } else {
        if (length(effect_i) == 1) {
          if (effect_i < 1) {
            effect_i <- length(self$ALASCA$loading)
            lapply(effect_i, function(i) self$ALASCA$loading[[i]][PC %in% component])
          } else {
            list(self$ALASCA$loading[[effect_i]][PC %in% component])
          }
        } else {
          lapply(effect_i, function(i) self$ALASCA$loading[[i]][PC %in% component])
        }
      }
    },
    #' @description
    #' @param n_limit Return only two times this number of coefficients (the `n_limit` highest and lowest coefficients). Use `0` to return all (default)
    get_covars = function(n_limit = 0) {
      if (n_limit > 0) {
        rbindlist(lapply(unique(self$covar_coefficients$variable), function(x) {
          self$covar_coefficients[variable == x][c(
            Rfast::nth(estimate, k = n_limit, num.of.nths = n_limit, descending = FALSE, index.return = TRUE),
            Rfast::nth(estimate, k = n_limit, num.of.nths = n_limit, descending = TRUE, index.return = TRUE)
          )]
        }))
      } else {
        return(self$covar_coefficients)
      }
    },
    #' @description
    #' Print ALASCA version
    #' @return String with version number and date of latest update
    print_version = function() {
      return(
        paste0("ALASCA (v", self$version, ", ", self$update_date, ")")
      )
    },
    do_validate = do_validate,
    get_scaling_function = get_scaling_function,
    get_default_scaling_function = get_default_scaling_function,
    get_pca_function = get_pca_function,
    build_model = build_model,
    remove_embedded_data = remove_embedded_data,
    do_reduce_dimensions = do_reduce_dimensions,
    run_regression = run_regression,
    get_regression_coefficients = get_regression_coefficients,
    remove_covars = remove_covars,
    get_effect_matrix = get_effect_matrix,
    do_pca = do_pca,
    clean_pca = clean_pca,
    clean_alasca = clean_alasca,
    get_regression_predictions = get_regression_predictions,
    get_bootstrap_data = get_bootstrap_data,
    prepare_validation_run = prepare_validation_run,
    get_validation_percentiles = get_validation_percentiles,
    get_validation_percentiles_loading = get_validation_percentiles_loading,
    get_validation_percentiles_score = get_validation_percentiles_score,
    get_validation_percentiles_regression = get_validation_percentiles_regression,
    get_validation_percentiles_covars = get_validation_percentiles_covars,
    #' @description
    #' Write scores, loadings, covars and predictions from validation run to database and remove data from memory
    #' @param ii Number of the validation run
    save_validation = function(ii) {
      for (i in seq_along(self$ALASCA$score)) {
        DBI::dbWriteTable(self$db_con, "loading",
          data.table(
            self$ALASCA$loading[[i]],
            model = ii,
            effect = i
          )[, covars := factor(covars, levels = self$get_levels("variable", reduced = FALSE))][, covars := as.integer(covars)],
          append = TRUE
        )

        scores_to_append <- data.table(
          self$ALASCA$score[[i]],
          model = ii,
          effect = i
        )

        for (new_cols in self$effect_terms[!self$effect_terms %in% colnames(scores_to_append)]) {
          set(scores_to_append, j = new_cols, value = self$get_ref(new_cols))
        }

        DBI::dbWriteTable(self$db_con, "score",
          scores_to_append,
          append = TRUE
        )

        DBI::dbWriteTable(self$db_con, "explained",
          data.table(
            self$ALASCA$explained[[i]],
            model = ii,
            effect = i
          ),
          append = TRUE
        )
      }
      self$ALASCA <- NULL

      DBI::dbWriteTable(self$db_con, "prediction",
        data.table(
          self$model_prediction,
          model = ii
        )[, variable := factor(variable, levels = self$get_levels("variable", reduced = FALSE))][, variable := as.integer(variable)],
        append = TRUE
      )
      self$model_prediction <- NULL

      if (length(self$covar_coefficients) > 0) {
        DBI::dbWriteTable(self$db_con, "covars",
          data.table(
            self$covar_coefficients,
            model = ii
          )[, covar := factor(covar, levels = self$get_levels("variable", reduced = FALSE))][, covar := as.integer(covar)],
          append = TRUE
        )
        self$covar_coefficients <- NULL
      }


      #   for (i in seq_along(self$ALASCA$score)) {
      #     fst::write_fst(self$ALASCA$score[[i]], paste0("val/iteration_", ii, "_effect_", i, "_score.fst"))
      #     fst::write_fst(self$ALASCA$loading[[i]], paste0("val/iteration_", ii, "_effect_", i, "_loading.fst"))
      #   }
      #   fst::write_fst(self$model_prediction, paste0("val/iteration_", ii, "_effect_", i, "_model_prediction.fst"))
    },
    get_validation_scores = get_validation_scores,
    get_validation_loadings = get_validation_loadings,
    #' @description
    #' Save scores, loading, covars, and predictions to csv files
    save_to_csv = function() {
      for (i in seq_along(self$ALASCA$loading)) {
        fwrite(
          self$ALASCA$loading[[i]],
          paste0(self$filepath, "loadings_effect_", i, ".csv")
        )
        fwrite(
          self$ALASCA$score[[i]],
          paste0(self$filepath, "scores_effect_", i, ".csv")
        )
        expl <- data.table(explained = self$ALASCA$explained[[i]])
        expl$PC <- seq_len(nrow(expl))
        fwrite(
          expl,
          paste0(self$filepath, "explained_effect_", i, ".csv")
        )
      }
      fwrite(
        self$model_prediction,
        paste0(self$filepath, "model_prediction.csv")
      )
      if (length(self$covar_coefficients) > 0) {
        fwrite(
          self$covar_coefficients,
          paste0(self$filepath, "covars.csv")
        )
      }
    },
    #' @description
    #' Save model to RDS and scores, loading, covars, and predictions to csv files
    save_model = function() {
      self$log("Exporting data")
      self$save_to_csv()
      saveRDS(self, paste0(self$filepath, self$filename, ".RDS"))
    },
    #' @description
    #' Rotate model loadings and scores with procrustes. This function checks all combinations of signs to minimize variation
    #' @param target Rotate model (`self`) with this as target
    rotate_matrix_optimize_score = function(target) {
      target$log("Starting rotation", level = "DEBUG")

      # PCA can give loadings with either sign, so we have to check whether switching signs improves the rotation
      for (effect_i in seq_along(self$ALASCA$loading)) {
        # Number of components to look at
        cols_to_look_at <- paste0("PC", self$get_PCs(effect_i))
        N <- length(cols_to_look_at)

        # Create matrix with all possible combinations of signs
        vec <- c(-1, 1)
        lst <- lapply(numeric(N), function(x) vec)
        signMatrix <- as.matrix(expand.grid(lst))
        
        # Test all combinations and calculate residuals
        signVar <- vapply(seq_len(nrow(signMatrix) / 2), function(i) {
          c <- .procrustes(
            loadings = as.matrix(
              t(t(self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], ..cols_to_look_at]) * signMatrix[i, ])
            ),
            target = as.matrix(target$ALASCA$loading[[effect_i]][, ..cols_to_look_at])
          )
          sum((target$ALASCA$score[[effect_i]][, ..cols_to_look_at] -
            as.matrix(t(t(self$ALASCA$score[[effect_i]][target$ALASCA$score[[effect_i]], ..cols_to_look_at]) * signMatrix[i, ])) %*% solve(c$t1))^2)
        }, FUN.VALUE = numeric(1))

        # Find the combination that minimizes the sum of squares
        minSignVar <- which(signVar == min(signVar))[1]

        # Switch signs
        for (i in seq_along(cols_to_look_at)) {
          set(self$ALASCA$loading[[effect_i]], j = cols_to_look_at[i], value = self$ALASCA$loading[[effect_i]][, get(cols_to_look_at[i])] * signMatrix[minSignVar, i])
          set(self$ALASCA$score[[effect_i]], j = cols_to_look_at[i], value = self$ALASCA$score[[effect_i]][, get(cols_to_look_at[i])] * signMatrix[minSignVar, i])
        }

        # Final rotation
        rotated_loadings <- .procrustes(
          loadings = as.matrix(self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], ..cols_to_look_at]),
          target = as.matrix(target$ALASCA$loading[[effect_i]][, ..cols_to_look_at])
        )

        self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], (cols_to_look_at) := as.data.table(rotated_loadings$procrust)]
        self$ALASCA$score[[effect_i]][target$ALASCA$score[[effect_i]], (cols_to_look_at) := as.data.table(as.matrix(.SD) %*% solve(rotated_loadings$t1)), .SDcols = cols_to_look_at]
      }

      target$log("Completed rotation", level = "DEBUG")

      # invisible(self)
    },
    #' @description
    #' Check if components have to be re-ordered to better match the target
    #' @param target Compare components of the model (`self`) with this as target
    optimize_components = function(target) {
      target$log("Comparing the order of PCs", level = "DEBUG")

      for (effect_i in seq_along(self$ALASCA$loading)) {
        # Number of components to look at
        cols_to_look_at <- paste0("PC", self$get_PCs(effect_i))
        
        # Compute the cost matrix for column reordering
        cost_matrix <- 1 - abs(
          cor(
            target$ALASCA$loading[[effect_i]][, ..cols_to_look_at],
            self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], ..cols_to_look_at]
          )
        )
        
        # Apply the Hungarian algorithm to find the optimal column permutation
        column_permutation <- clue::solve_LSAP(cost_matrix)
        
        if ( any(column_permutation != seq(1, length(cols_to_look_at))) ) {
          target$log(paste("Changing the order of PCs for effect", effect_i, ":", column_permutation), level = "WARNING")
          column_names <- colnames(self$ALASCA$loading[[effect_i]])
          setcolorder(self$ALASCA$loading[[effect_i]], neworder = paste0("PC", column_permutation))
          colnames(self$ALASCA$loading[[effect_i]]) <- column_names
        }
        
      }
    },
    #' @description
    #' Rotate model loadings and scores with procrustes
    #' @param target Rotate model (`self`) with this as target
    rotate_matrix = function(target) {
      target$log("Starting rotation (not optimized)", level = "DEBUG")
      for (effect_i in seq_along(self$ALASCA$loading)) {

        # Number of components to look at
        cols_to_look_at <- paste0("PC", self$get_PCs(effect_i))

        rotated_loadings <- .procrustes(
          loadings = as.matrix(self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], ..cols_to_look_at]),
          target = as.matrix(target$ALASCA$loading[[effect_i]][, ..cols_to_look_at])
        )

        self$ALASCA$loading[[effect_i]][target$ALASCA$loading[[effect_i]], (cols_to_look_at) := as.data.table(rotated_loadings$procrust)]
        self$ALASCA$score[[effect_i]][target$ALASCA$score[[effect_i]], (cols_to_look_at) := as.data.table(as.matrix(.SD) %*% solve(rotated_loadings$t1)), .SDcols = cols_to_look_at]
      }

      target$log("Completed rotation", level = "DEBUG")

      # invisible(self)
    }
  ),
  active = list(
    #' @field get_plot_group Name of the grouping factor (used for plotting)
    get_plot_group = function() {
      if (is.null(self$splot$group)) {
        if (length(self$effect_terms) == 1) {
          # The first term is used for time
          self$splot$group <- self$stratification_column
        } else if (length(self$effect_terms) == 2) {
          self$splot$group <- self$effect_terms[2]
        } else {
          self$splot$make_group_column <- TRUE
          self$splot$palette <- NULL
          self$splot$linetypes <- NULL
          self$splot$group <- paste(self$effect_terms[-1], collapse = "-")
          self$df_raw$level_list[[self$splot$group]] <- self$df[, unique(do.call(paste, c(.SD, sep = "-"))), .SDcols = self$effect_terms[-1]]
        }
        self$splot$group_label <- self$splot$capitalize(self$splot$group)
      }
      self$splot$group
    },
    #' @field effect_terms List of the terms in the effect matrices
    effect_terms = function() {
      unique(unlist(self$effect_list$terms))
    }
  ),
  private = list()
)
