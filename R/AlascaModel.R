# source("R/buildModel.R")
# source("R/export.R")
# source("R/pca.R")
# source("R/validation.R")
AlascaModel <- R6::R6Class("AlascaModel",
  class = FALSE,
  lock_objects = FALSE,
  public = list(
    df = NULL,
    formula = NULL,
    wide = FALSE,
    separate_time_and_group = FALSE,
    scale_function = "sdall",
    reduce_dimensions = FALSE,
    equal_baseline = FALSE,
    sum_coding = FALSE,
    method = NA,
    use_Rfast = TRUE,
    p_adjust_method = NA,
    participant_column = "ID",
    x_column = "time",
    ignore_missing = FALSE,
    ignore_missing_covars = FALSE,
    silent = FALSE,
    scale_function.center = FALSE,
    stratification_column = "group",
    stratification_vector = NA,
    minimize_object = FALSE,
    limitsCI = c(0.025, 0.975),
    pca_function = "prcomp",
    plot.x_label = "Time",
    plot.group_label = "Group",
    plot.figsize = c(180, 120, 300),
    plot.figunit = "mm",
    plot.filetype = "png",
    plot.palette = NA,
    plot.loadinggroupcolumn = NA,
    plot.loadinggroup_label = "Variable group",
    plot.palette.end = 0.8,
    explanatorylimit = 0.05,
    do_debug = FALSE,
    n_validation_folds = 7,
    n_validation_runs = 1000,
    validation_quantile_method = 2,
    plot.my_theme = ggplot2::theme_classic(),
    keep_columns = c(""),
    keep_terms = c(""),
    filename = NA,
    filepath = NA,
    reduce_dimensions.nComps = NULL,
    reduce_dimensions.limit = 0.95,
    save = FALSE,
    save_to_disk = FALSE,
    save_validation_ids = FALSE,
    optimize_score = TRUE,
    validate = FALSE,
    validate_regression = TRUE,
    validation = FALSE,
    validation_method = "bootstrap",
    validation_ids = NA,
    validation_object = NA,
    validation_assign_new_ids = FALSE,
    validation_participants = NA,
    init_time = Sys.time(),
    rawFormula = NULL,
    keep_validation_objects = FALSE,
    log_file = tempfile(),
    log_level = NULL,
    log = NULL,
    valCol = NULL,
    formula_terms = NULL,
    all_formula_terms = NULL,
    hasGroupTerm = NULL,
    hasInteractionTerm = NULL,
    covars = NULL,
    variablelist = NULL,
    timelist = NULL,
    grouplist = NULL,
    new_formula = NULL,
    my_df_rows = NULL,
    initialize = function(df, formula, ...) {
      
      inputs <- list(...)
      for (i in seq_along(inputs)) {
        self[[names(inputs)[i]]] <- inputs[[i]]
      }
      
      self$df <- setDT(df)
      self$formula <- formula
      self$validate <- self$validate || self$validation
      self$validation <- NULL
      self$n_validation_runs <- ifelse(any(is.na(self$validation_ids)), self$n_validation_runs, nrow(self$validation_ids))
      self$filepath <- ifelse(is.na(self$filepath), NA, ifelse(substr(self$filepath, nchar(self$filepath), nchar(self$filepath)) == "/", self$filepath,  paste0(self$filepath, "/")))
      self$save_to_disk <- self$validate && self$save_to_disk
      self$rawFormula <- self$formula
      #self$ALASCA.version <- print_version(get = "version")
      #self$ALASCA.version.date <- print_version(get = "date")
      self$log_level <- ifelse(self$do_debug, "DEBUG", "INFO")
      
      if (self$silent) {
        self$log <- log4r::logger(self$log_level, appenders = list(log4r::file_appender(self$log_file)))
      } else {
        self$log <- log4r::logger(self$log_level, appenders = list(log4r::console_appender(), log4r::file_appender(self$log_file)))
      }
      
      log4r::info(self$log, "Initializing ALASCA.")
      
      # Clean input ----
      log4r::info(self$log, "Has initialized the ALASCA model. Next step is to clean it and check input")
      self$sanitize_object()
      
      # Build the ALASCA model ----
      self$build_model()
      
      
    },
    update = function() {

      ## Unscaled values
      
      if(self$reduce_dimensions){
        self$Limm$main$pca <- self$Limm$pca
        self$Limm$pca <- NULL
      }
      
      ## Avoid recursion
      self$validate <- FALSE
      
      ## Save space
      self$minimize_object <- TRUE
      
      self$variablelist <- unique(self$df$variable)
      self$timelist <- levels(self$df$time)
      self$grouplist <- levels(self$df$group)
      
      # Build the ALASCA model ----
      
      self$build_model()
      
      # To save space, we remove unnecessary embedded data ----
        log4r::debug(self$log, "Starting to remove embedded data")
        self$remove_embedded_data()
        log4r::debug(self$log, "Completed remove embedded data")
      
      #invisible(self)
    },
    sanitize_object = sanitize_object,
    get_info_from_formula = get_info_from_formula,
    LM_or_LMM = LM_or_LMM,
    do_validate = do_validate,
    rename_columns_to_standard = rename_columns_to_standard,
    wide_to_long = wide_to_long,
    check_that_columns_are_valid = check_that_columns_are_valid,
    adjust_design_matrix = adjust_design_matrix,
    get_scaling_function = get_scaling_function,
    get_default_scaling_function = get_default_scaling_function,
    get_pca_function = get_pca_function,
    find_missing_predictor_variables = find_missing_predictor_variables,
    find_missing_response_variables = find_missing_response_variables,
    build_model = build_model,
    remove_embedded_data = remove_embedded_data,
    saveALASCA = saveALASCA,
    do_reduce_dimensions = do_reduce_dimensions,
    run_regression = run_regression,
    get_relevant_pcs = get_relevant_pcs,
    get_regression_coefficients = get_regression_coefficients,
    remove_covars = remove_covars,
    separate_regression_coefficients = separate_regression_coefficients,
    get_effect_matrix = get_effect_matrix,
    do_pca = do_pca,
    clean_pca = clean_pca,
    clean_alasca = clean_alasca,
    get_regression_predictions = get_regression_predictions,
    get_bootstrap_ids = get_bootstrap_ids,
    get_bootstrap_data = get_bootstrap_data,
    prepare_validation_run = prepare_validation_run,
    rotate_matrix_optimize_score = rotate_matrix_optimize_score,
    rotate_matrix = rotate_matrix,
    get_validation_percentiles = get_validation_percentiles,
    get_validation_percentiles_loading = get_validation_percentiles_loading,
    get_validation_percentiles_score = get_validation_percentiles_score,
    get_validation_percentiles_regression = get_validation_percentiles_regression,
    get_validation_percentiles_covars = get_validation_percentiles_covars,
    get_covars = get_covars
  ),
  private = list()
)
