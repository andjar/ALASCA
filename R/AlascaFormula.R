AlascaFormula <- R6::R6Class("AlascaFormula",
  class = FALSE,
  public = list(
    df = NULL,
    model = NULL,
    raw_formula = NULL,
    formula = NULL,
    regression_formula = NULL,
    formula_wo_random = NULL,
    fixed_terms = NULL,
    random_terms = NULL,
    additional_terms = NULL,
    compatible_with_Rfast = NULL,
    #' @field lhs String. Left hand side of regression formula
    lhs = NULL,
    #' @field rhs String. Right hand side of regression formula
    rhs = NULL,
    covars = NULL,
    initialize = function(formula, model) {
      self$raw_formula <- formula
      self$model <- model
      self$formula <- formula
      self$update()
      if (is.null(self$model$effect_list$expr)) {
        private$guess_effects()
      }
      
      private$set_method()
      if (self$model$use_Rfast) {
        self$model$log("Will use Rfast!")
      }
      
    },
    add = function(new_term) {
      self$additional_terms <- unique(c(self$additional_terms, new_term))
    },
    replace = function(old_term, new_term) {
      #   self$random_terms[self$random_terms == old_term] <- new_term
      # }
      # if (!is.null(self$fixed_terms)) {
      #   self$fixed_terms[self$fixed_terms == old_term] <- new_term
      # }
      # if (!is.null(self$lhs)) {
      #   if (self$lhs == old_term) self$lhs <- new_term
      # }
      self$lhs <- gsub(old_term, new_term, self$lhs)
      self$rhs <- gsub(old_term, new_term, self$rhs)
      private$update_formula()
      private$find_terms()
    },
    update = function() {
      private$find_sides()
      private$find_terms()
      private$check_Rfast()
    },
    has_random = function() length(self$random_terms) > 0,
    has_covars = function() {
      self$covars <- self$fixed_terms[!self$fixed_terms %in% unlist(self$model$effect_list$terms)]
      length(self$covars) > 0
    },
    get_regression_formula = get_regression_formula
  ),
  active = list(
    all_terms = function() {
      unique(
        c(self$fixed_terms,
          unlist(strsplit(gsub("1\\||\\(|\\)", "", self$random_terms), split = "\\|")),
          self$additional_terms)
      )
    },
    ID = function() {
      if (is.null(self$model$participant_column)) {
        if (length(self$random_terms) == 0) {
          self$model$df_raw$df[, ID := seq_len(nrow(self$model$df_raw$df))]
          self$model$participant_column <- "ID"
          self$model$participant_column
        } else if (length(self$random_terms) > 1) {
          self$model$log("Unable to find ID", level = "ERROR")
          stop()
        } else if (!grepl("1|", self$random_terms, fixed = TRUE)) {
          self$model$log("Unable to find ID", level = "ERROR")
          stop()
        } else {
          gsub("\\(|\\)|1\\|","",self$random_terms[[1]])
        }
      } else {
        self$model$participant_column
      }
    },
    all_fixed_terms = function() {
      unique(
        c(self$fixed_terms,
          self$additional_terms)
      )
    }
  ),
  private = list(
    find_sides= function() {
      self$lhs <- as.character(self$formula)[2]
      self$rhs <- gsub(" ", "", as.character(self$formula)[3])
    },
    find_terms = function() {
      str_terms <- unlist(strsplit(self$rhs, split = "\\:|\\+|\\*"))
      self$fixed_terms <- unique(str_terms[!grepl("|", str_terms, fixed = TRUE)])
      self$random_terms <- unique(str_terms[grepl("|", str_terms, fixed = TRUE)])
      
      str_terms <- unlist(strsplit(self$rhs, split = "+", fixed = TRUE))
      self$formula_wo_random <- formula(
        paste(self$lhs, "~", paste(str_terms[!grepl("|", str_terms, fixed = TRUE)], collapse = "+"))
      )
    },
    update_formula = function() {
      self$formula <- as.formula(paste(self$lhs, "~", self$rhs))
      str_terms <- unlist(strsplit(self$rhs, split = "+", fixed = TRUE))
      self$formula_wo_random <- formula(
        paste(self$lhs, "~", paste(str_terms[!grepl("|", str_terms, fixed = TRUE)], collapse = "+"))
      )
    },
    set_method = set_method,
    check_Rfast = function() {
      if (length(self$random_terms) > 1) {
        self$compatible_with_Rfast <- FALSE
      } else if (length(self$random_terms) == 1 && !grepl("1|", self$random_terms, fixed = TRUE)) {
        self$compatible_with_Rfast <- FALSE
      } else {
        self$compatible_with_Rfast <- TRUE
      }
    },
    guess_effects = function() {
      if (length(self$fixed_terms) == 1) {
        # Use the only fixed effect provided
        self$model$effect_list$expr <- self$fixed_terms[[1]]
      } else {
        # Use the first fixed term as basis
        str_all_terms <- colnames(attr(terms(self$formula_wo_random),"factors"))
        str_all_terms <- str_all_terms[!grepl("|", str_all_terms, fixed = TRUE)]
        
        # This will typically be a main effect and an interaction
        str_terms <- str_all_terms[grepl(self$fixed_terms[[1]], str_all_terms)]
        
        # Look for the main term of the interaction term and add other terms involving it
        if (!self$model$equal_baseline) {
          additional_terms <- unique(unlist(strsplit(str_terms, split= ":", fixed = TRUE)))
          additional_terms <- additional_terms[additional_terms != self$fixed_terms[[1]]]
          for(i in str_all_terms) {
            if (i %in% additional_terms) {
              str_terms <- c(str_terms, additional_terms)
            }
          }
        } 
        
        if (self$model$separate_effects) {
          self$model$effect_list$expr <- c(str_terms[1], paste(str_terms[-1], collapse = "+"))
        } else {
          self$model$effect_list$expr <- paste(str_terms, collapse = "+")
        }
      }
      self$model$log(paste0("Guessing effects: `", paste0(self$model$effect_list$expr, collapse = "` and `"), "`"), level = "WARN")
    }
  )
)
