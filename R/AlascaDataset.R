AlascaDataset <- R6::R6Class("AlascaDataset",
  class = FALSE,
  public = list(
    #' @field data_df The original dataset
    data_df = NULL,
    model = NULL,
    rows_to_serve = NULL,
    IDs = NULL,
    rows_by_ID = NULL,
    #' @field level_list List of factor levels by column
    level_list = NULL,
    #' @field column_types Type of columns in `data_df`
    column_types = NULL,
    reduced_df = list(
      loading = NULL
    ),
    initialize = function(data_df, model) {
      self$model <- model
      self$data_df <- setDT(data_df)

      # Get ID column
      if (self$model$formula$ID == "needs_to_add_ID") {
        # No mixed effects + the participant_column has not been properly defined
        if ("ID" %in% colnames(self$data_df)) {
          self$model$participant_column <- "ID"
        } else {
          if (self$model$wide) {
            self$data_df[ID := seq_len(nrow(self$data_df))]
          } else {
            self$model$log("Please provide a column with ID, either named `ID` or specify the `participant_column` argument", level = "ERROR")
            stop()
          }
        }
      }

      if (self$model$validate && typeof(self$data_df[[self$model$formula$ID]]) != "integer") {
        self$model$log("Converting IDs to integer values", level = "WARN")
        set(self$data_df, j = self$model$formula$ID, value = as.integer(as.factor(self$data_df[[self$model$formula$ID]])))
      }

      # Convert from wide to long
      if (self$model$wide) {
        private$wide_to_long()
      }

      self$data_df[, variable := factor(variable)]

      ## We need to keep original IDs to have a unique identifier later on
      if (self$model$validation_method == "bootstrap") {
        self$model$formula$add("originalIDbeforeBootstrap")
        self$data_df[, originalIDbeforeBootstrap := -1]

        self$model$formula$add("uniqueIDforBootstrap")
        self$data_df[, uniqueIDforBootstrap := -1]
      }

      # Find stratification column from formula
      private$guess_stratification_column()
      self$model$stratification_vector <- self$data_df[, get(self$model$stratification_column)]

      # Keep variable labels
      if (!is.null(self$model$splot$loading_group_column)) {
        self$model$formula$add(self$model$splot$loading_group_column)
        self$model$variable_labels <- unique(self$data_df[, .SD, .SDcols = c("variable", self$model$splot$loading_group_column)])
        colnames(self$model$variable_labels) <- c("covars", "covargroup")
      }

      # Remove surplus columns
      self$data_df <- self$data_df[, .SD, .SDcols = c(self$model$formula$all_terms, "variable", self$model$formula$lhs)]
      self$rows_to_serve <- seq_len(nrow(self$data_df))

      self$IDs <- self$data_df[[self$model$formula$ID]]
      self$rows_by_ID <- lapply(unique(self$IDs), function(x) which(x == self$IDs))
      names(self$rows_by_ID) <- unique(self$IDs)

      self$column_types <- vapply(self$data_df, class, FUN.VALUE = character(1))
      if (any(self$column_types == "character")) {
        self$model$log("Converting `character` columns to factors", level = "WARN")
        for (char_to_factor in names(self$column_types)[self$column_types == "character"]) {
          set(self$data_df, j = char_to_factor, value = as.factor(self$data_df[[char_to_factor]]))
        }
        self$column_types[self$column_types == "character"] <- "factor"
      }
      self$level_list <- lapply(self$data_df[, .SD, .SDcols = self$column_types == "factor"], function(x) levels(x))

      self$model$log("Checking for missing information", level = "DEBUG")
      private$check_against_formula()
      # private$find_missing_response_variables()
      private$find_missing_predictor_variables()
      if (self$model$save_to_disk) {
        fst::write_fst(self$data_df, paste0(self$model$filepath, self$model$filename, ".fst"))
        self$data_df <- NULL
      }
    },
    count_participants = function() {
      wide_df <- dcast(
        data = self$df,
        as.formula(paste(
          self$model$formula$ID, "+", self$model$stratification_column, "+variable",
          "~",
          paste(self$model$effect_terms[self$model$effect_terms != self$model$stratification_column], collapse = "+")
        )),
        fun.aggregate = length
      )
      samples_by_participant <- unique(
        data.frame(
          ID = wide_df[[self$model$formula$ID]],
          count = rowSums(wide_df[, 4:ncol(wide_df)])
        )
      )
      cases_by_group <- as.data.frame(
        table(
          merge(samples_by_participant,
            wide_df[!duplicated(get(self$model$formula$ID)), .SD, .SDcols = c(self$model$formula$ID, self$model$stratification_column)],
            all.x = TRUE, all.y = FALSE, by = self$model$formula$ID
          )[, c(self$model$stratification_column, "count")]
        )
      )

      list(
        participants_by_group = self$df[, .(count = uniqueN(get(self$model$formula$ID))), by = get(self$model$stratification_column)],
        cases_by_group = cases_by_group,
        participants_by_effect_1 = self$df[, .(count = uniqueN(get(self$model$formula$ID))), by = get(self$model$effect_terms[[1]])],
        participants_by_effect_1_and_group = ifelse(length(self$model$effect_terms) < 2, NULL,
          dcast(
            data = self$df[, .(count = uniqueN(get(self$model$formula$ID))), by = c(self$model$effect_terms)],
            formula(paste(paste(self$model$effect_terms[[-1]], collapse = "+"), "~", self$model$effect_terms[[1]])),
            value.var = "count"
          )
        ),
        participants_by_variable_and_effect_1_and_group = dcast(
          data = self$df[, .(count = uniqueN(get(self$model$formula$ID))), by = c("variable", self$model$effect_terms)],
          formula(paste(paste(self$model$effect_terms, collapse = "+"), "~ variable")),
          value.var = "count"
        ),
        samples_by_variable_and_effect_1_and_group = dcast(
          data = self$df[, .(count = .N), by = c("variable", self$model$effect_terms)],
          formula(paste(paste(self$model$effect_terms, collapse = "+"), "~ variable")),
          value.var = "count"
        ),
        samples_by_participant = samples_by_participant,
        samples_by_participant_and_effect_1 = dcast(
          data = self$df,
          formula(paste(self$model$formula$ID, "+", self$model$effect_terms[[1]], "~ .")),
          value.var = "variable",
          fun.aggregate = length
        )
      )
    }
  ),
  active = list(
    df = function() {
      if (self$model$save_to_disk && is.null(self$data_df)) {
        setDT(fst::read_fst(paste0(self$model$filepath, self$model$filename, ".fst")))
      } else {
        self$data_df
      }
    }
  ),
  private = list(
    wide_to_long = function() {
      self$model$log("Converting from wide to long!")

      self$data_df <- melt(self$data_df,
        id.vars = self$model$formula$all_terms, value.name = self$model$formula$lhs
      )

      self$model$log(paste0("Found ", length(unique(self$data_df$variable)), " variables"))
      self$model$wide <- FALSE

      # invisible(self)
    },
    check_against_formula = function() {
      if (any(!self$model$formula$all_terms %in% colnames(self$df))) {
        self$model$log("Some formula terms does match with the provided data", level = "ERROR")
        stop()
      }
    },
    guess_stratification_column = function() {
      # Check or guess stratification column
      if (!is.null(self$model$stratification_column) && !self$model$stratification_column %in% self$model$formula$all_fixed_terms) {
        self$formula$add(self$model$stratification_column)
        self$model$log(paste0("The `", self$model$stratification_column, "` column is used for stratification"))
      } else if (is.null(self$model$stratification_column) && "group" %in% colnames(self$data_df)) {
        self$model$formula$add("group")
        self$model$stratification_column <- "group"
        self$model$log("The `group` column is used for stratification", level = "WARN")
      } else if (is.null(self$model$stratification_column)) {
        self$model$stratification_column <- self$model$formula$all_fixed_terms[[1]]
        self$model$log(paste0("The `", self$model$formula$all_fixed_terms[[1]], "` column is used for stratification"), level = "WARN")
      }
    },
    find_missing_response_variables = find_missing_response_variables,
    find_missing_predictor_variables = find_missing_predictor_variables
  )
)
