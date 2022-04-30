AlascaPlot <- R6::R6Class("AlascaPlot",
    class = FALSE,
    public = list(
      model = NULL,
      my_theme = ggplot2::theme_classic() + ggplot2::theme(legend.position = "bottom"),
      variable_label = "Variable",
      variable = NULL,
      x_label = NULL,
      group_label = "Group",
      n_bins = NULL,
      group = NULL,
      ribbon = TRUE,
      dodgewidth = 0.5,
      height = NULL,
      width = NULL,
      dheight = NULL,
      dwidth = NULL,
      dpi = 300,
      units = "mm",
      filetype = "png",
      make_group_column = FALSE,
      palette = NULL,
      file_counter = 1,
      linetypes = NULL,
      loading_group_column = NULL,
      loadinggroup_label = "Variable group",
      sort_by_loading_group = TRUE,
      palette.end = 0.8,
      effect_i = 1,
      component = 1,
      n_col = 1,
      save = FALSE,
      flip_axis = TRUE,
      x_angle = 45,
      x_v_just = 1,
      x_h_just = 1,
      n_limit = 12,
      labels = "AUTO",
      type = "effect",
      initialize = function(model) {
        self$model <- model
      },
      call_plot = function(...) {
        inputs <- list(...)
        for (i in seq_along(inputs)) {
          self[[names(inputs)[i]]] <- inputs[[i]]
        }
        
        private$check_call()
        
        self$model$log(
          paste0(self$capitalize(self$type), " plot. Selected effect (nr ", paste0(self$effect_i, collapse = " and "), "): `", paste0(unlist(lapply(self$effect_i, function(i) self$model$effect_list$expr[[i]])), collapse = "` and `"), "`. Component: ", paste(self$component, collapse = " and "), ".")
        )
        
        if (self$type == "effect") {
          self$dheight <- 120*length(self$effect_i)*length(self$component)
          self$dwidth <- 180
          g <- self$plot_effect()
        } else if (self$type == "score") {
          self$dheight <- 90
          self$dwidth <- 90
          g <- self$plot_effect_score(effect_i = self$effect_i, component = self$component)
        } else if (self$type == "loading") {
          self$dheight <- 120*length(self$effect_i)*length(self$component)
          self$dwidth <- 90
          g <- self$plot_effect_loading(effect_i = self$effect_i, component = self$component)
        } else if (self$type == "validation") {
          self$dheight <- 120*length(self$effect_i)*length(self$component)
          self$dwidth <- 180
          g <- self$plot_effect_validation()
        } else if (self$type == "2D") {
          self$dheight <- 240
          self$dwidth <- 240
          g <- self$plot_2D()
        } else if (self$type == "histogram") {
          self$dheight <- 180
          self$dwidth <- 180
          if (is.null(self$n_bins)) {
            self$n_bins <- self$model$n_validation_runs/5
          }
          g <- self$plot_histogram()
        } else if (self$type == "scree") {
          self$dheight <- 60
          self$dwidth <- 60
          g <- self$plot_scree()
        } else if (self$type == "covars") {
          self$dheight <- 90
          self$dwidth <- 90
          g <- self$plot_covars()
        } else if (self$type == "residuals") {
          self$dheight <- 180
          self$dwidth <- 180
          g <- self$plot_residuals()
        } else if (self$type == "prediction") {
          self$dheight <- 120
          self$dwidth <- 180
          g <- self$plot_prediction()
        } else {
          self$model$log(paste("Unkown plot type:", self$type), level = "ERROR")
          stop()
        }
        
        private$post_process(g)
        return(g)
      },
      plot_effect = plot_effect,
      plot_effect_score = plot_effect_score,
      plot_effect_loading = plot_effect_loading,
      plot_effect_validation = plot_effect_validation,
      plot_effect_validation_score = plot_effect_validation_score,
      plot_effect_validation_loading = plot_effect_validation_loading,
      plot_histogram = plot_histogram,
      plot_histogram_score = plot_histogram_score,
      plot_histogram_loading = plot_histogram_loading,
      plot_2D = plot_2D,
      plot_2D_score = plot_2D_score,
      plot_2D_loading_1 = plot_2D_loading_1,
      plot_2D_loading_2 = plot_2D_loading_2,
      plot_scree = function() {
        if (length(self$effects) > 1) {
          g_list <- lapply(self$effects, function(i) {
            data_to_plot <- data.table(explained = self$model$ALASCA$explained[[i]])
            data_to_plot$PC <- seq_len(nrow(data_to_plot))
            data_to_plot <- data_to_plot[PC <= max(self$component)]
            g <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = factor(PC), y = explained, group = NA)) +
              ggplot2::geom_point() +
              ggplot2::geom_line() +
              ggplot2::labs(x = "Principal component", y = "Variance explained") +
              self$my_theme
            g
          })
          do.call(ggpubr::ggarrange, list(plotlist = g_list, labels = self$labels))
        } else {
          data_to_plot <- data.table(explained = self$model$ALASCA$explained[[self$effects]])
          data_to_plot$PC <- seq_len(nrow(data_to_plot))
          data_to_plot <- data_to_plot[PC <= max(self$component)]
          ggplot2::ggplot(data_to_plot, ggplot2::aes(x = factor(PC), y = explained, group = NA)) +
            ggplot2::geom_point() +
            ggplot2::geom_line() +
            ggplot2::labs(x = "Principal component", y = "Variance explained") +
            self$my_theme
        }
        
      },
      plot_residuals = function() {
        data_to_plot <- self$model$get_residuals(variable = self$variable)
        if (length(unique(data_to_plot$variable)) > 20) {
          self$model$log("Trying to plot too many variables. Please provide a list with < 20 elements as `variable`", level = "ERROR")
          stop()
        }
        ggplot2::ggplot(
          data_to_plot,
          ggplot2::aes(sample = residuals)
        ) + 
          ggplot2::stat_qq() +
          ggplot2::stat_qq_line() +
          ggplot2::facet_wrap(~variable) +
          ggplot2::labs(x = "Theoretical", y = "Sample") +
          self$my_theme
      },
      plot_prediction = plot_prediction,
      plot_covars = function() {
        
        # Note: Here `variable` refers to the regression term and `covar` to the variable/marker/...
        
        if (self$n_limit > 0) {
          self$model$log(paste("Only showing", self$n_limit*2, "variables. Adjust the number with `n_limit`"), level = "WARN")
        }
        data_to_plot <- self$model$get_covars(n_limit = self$n_limit)
        
        # Prettify terms
        all_variables <- unique(data_to_plot[, variable])
        for (i in self$model$formula$covars) {
          term_to_look_at <- all_variables[substr(all_variables, 1, nchar(i)) == i]
          for(j in term_to_look_at) {
            if (nchar(j) > nchar(i)) {
              # This term needs some processing (probably a factor)
              data_to_plot[variable == j, variable := self$prettify_covar(effect = i, txt = j)]
            }
          }
        }
        
        # Sort terms
        data_to_plot[, covar := factor(covar, levels = covar[order(estimate, decreasing = TRUE)])]
        
        if (is.null(self$loading_group_column)) {
          if(self$model$validate) {
            g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covar", y = "estimate", ymin = "low", ymax = "high")) + ggplot2::geom_pointrange()
          } else {
            g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covar", y = "estimate")) + ggplot2::geom_point()
          }
        } else {
          data_to_plot <- merge(data_to_plot, self$model$variable_labels, by.x = "covar", by.y = "covars")
          if (self$sort_by_loading_group) {
            data_to_plot[, covar := factor(covar, levels = covar[order(covargroup, estimate, decreasing = TRUE)])]
          } 
          if(self$model$validate) {
            g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covar", y = "estimate", ymin = "low", ymax = "high", shape = "covargroup", color = "covargroup")) + ggplot2::geom_pointrange()
          } else {
            g <- ggplot2::ggplot(data_to_plot, ggplot2::aes_string(x = "covar", y = "estimate", shape = "covargroup", color = "covargroup")) + ggplot2::geom_point()
          }
        }
        
        if (!is.null(self$loading_group_column)) {
          g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
            ggplot2::labs(color = self$loadinggroup_label, shape = self$loadinggroup_label)
        }
        
        g <- g + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
          ggplot2::labs(x = "Covariate", y = "Coefficient") + 
          ggplot2::facet_wrap(~variable, scales = "free_y") + 
          self$my_theme + self$xflip
        
        return(g)
      },
      capitalize = function(txt) {
        paste0(toupper(substr(txt, 1, 1)), substr(txt, 2, nchar(txt)))
      },
      prettify_covar = function(effect, txt) {
        partA <- substr(txt, 1, nchar(effect))
        partB <- substr(txt, 1+ nchar(effect), nchar(txt))
        paste0(self$capitalize(partA), ": ", self$capitalize(partB))
      },
      get_plot_linetypes = function() {
        if (is.null(self$linetypes)) {
          self$linetypes <- scales::linetype_pal()(length(self$get_levels(self$get_plot_group)))
          names(self$linetypes) <- self$get_levels(self$get_plot_group)
        }
        return(self$linetypes)
      },
      get_plot_palette = function() {
        if (is.null(self$palette)) {
          self$palette <- scales::viridis_pal(end = self$palette.end)(length(self$get_levels(self$get_plot_group)))
          names(self$palette) <- self$get_levels(self$get_plot_group)
        }
        return(self$palette)
      },
      get_explained_label = function(effect_i = 1, component = 1, type= "Score") {
        paste0(type, " PC", component, " (", round(100 * self$model$ALASCA$explained[[effect_i]][[component]], 2), "%)")
      },
      get_levels = function(x) self$model$get_levels(x),
      get_ref = function(x) self$model$get_ref(x)
    ),
    active = list(
      validate = function() self$model$validate,
      h = function() ifelse(is.null(self$height), self$dheight, self$height),
      w = function() ifelse(is.null(self$width), self$dwidth, self$width),
      xflip = function() {
        if (self$flip_axis) {
          ggplot2::coord_flip()
        } else {
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = self$x_angle, vjust = self$x_v_just, hjust = self$x_h_just))
        }
      },
      effects = function() {
        if (length(self$effect_i) == 1 && self$effect_i == 0) {
          self$effect_i <- seq_along(self$model$ALASCA$loading)
        }
        self$effect_i
      },
      get_plot_group = function() self$model$get_plot_group
    ),
    private = list(
      check_call = function() {
        if (max(self$effect_i) > length(self$model$ALASCA$loading)) {
          self$model$log("The effect you wanted to plot does not exist!", level = "ERROR")
          stop()
        }
        self$n_limit <- min(self$n_limit, length(self$model$get_levels("variable"))/2)
      },
      post_process = function(g) {
        if (self$save || self$model$save) {
          fname <- paste0(ifelse(is.null(self$model$filepath), "", self$model$filepath), "plot/", sprintf("%02d", self$file_counter), "-", self$model$filename, "_", self$type, ".", self$filetype)
          ggplot2::ggsave(
            plot = g,
            filename = fname,
            height = self$h,
            width = self$w,
            dpi = self$dpi,
            units = self$units
          )
          self$model$log(paste("Plot saved:", fname))
          self$file_counter <- self$file_counter + 1
        }
      }
    )
)

