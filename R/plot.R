#' Get an ALASCA object
#'
#' This function plots your ALASCA model
#'
#' @inheritParams plotDevelopment
#' @inheritParams plotComponents
#' @return An ggplot2 object (or a list og ggplot objects)
#'
#' @examples
#' load("PE.Rdata")
#' plot(model.val)
#' plot(model, component = "PC2")
#' plot(model, only = "score", effect = "time")
#' plot(model, tooDense = 5)
#' plot(model, highlight = c("PlGF", "IL-1b", "IL-6"))
#' @export
plot.ALASCA <- function(object,
                        component = 1,
                        ...) {
  if (length(component) == 1) {
    plotDevelopment(object = object, component = component, ...)
  } else if (length(component) == 2) {
    plotComponents(object = object, comps = component, ...)
  } else {
    stop("Please provide exactly 1 or 2 components to plot")
  }
}



#' Get an ALASCA object
#'
#' This function plots your ALASCA model
#'
#' @param object An [ALASCA()] object
#' @param component Integer stating which component to return (1 is default)
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param decreasingLoadings Sort the loadings in decreasing (`TRUE`, default) or increasing order (`FALSE`)
#' @param only String stating which plot to return; `both` (default), `score` or `loading`
#' @param enlist Logical. If `TRUE`, the plots are returned as a list and not as a composed figure (default)
#' @param tooDense Integer, If > 0, only name this number of covariables
#' @param xlabel Defaults to "Time" if not specified here or during model setup
#' @param grouplabel Defaults to "Group" if not specified here or  during model setup
#' @param flipaxes When `TRUE` (default), list the variable loadings vertical instead of horizontal
#' @param plotzeroline When `TRUE` (default), plot a zero line in the loading plot
#' @param limit_loading Only list robust loadings
#' @param filetype Which file type you want to save the figure to (default: `png`)
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit Unit for figure size (default: `mm`)
#' @param highlight Vector of strings with variables to highlight in the loadings plot
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return An ggplot2 object (or a list og ggplot objects)
#'
#' @examples
#' load("PE.Rdata")
#' plot(model.val)
#' plot(model, component = "PC2")
#' plot(model, only = "score", effect = "time")
#' plot(model, tooDense = 5)
#' plot(model, highlight = c("PlGF", "IL-1b", "IL-6"))
#' @export
plotDevelopment <- function(object,
                            component = 1,
                            effect = "both",
                            decreasingLoadings = TRUE,
                            only = "both",
                            enlist = FALSE,
                            tooDense = NA,
                            n.limit = 0,
                            highlight = NA,
                            xlabel = NA,
                            grouplabel = NA,
                            flipaxes = TRUE,
                            plotzeroline = TRUE,
                            filename = NA,
                            filetype = NA,
                            figsize = NA,
                            variables = NA,
                            dodgewidth = 0.35,
                            figunit = NA,
                            plotribbon = TRUE,
                            loadinggroup = NA,
                            limit_loading = FALSE,
                            sortbyloadinggroup = TRUE,
                            myTheme = object$plot.myTheme) {
  
  if (!(effect %in% c("both", "time", "group"))) stop("`effect` has to be `both`, `time` or `group`")
  if (!object$separateTimeAndGroup) effect <- "time"
  if (!is.na(xlabel)) object$plot.xlabel <- xlabel
  if (!is.na(filename)) object$filename <- filename
  if (!is.na(grouplabel)) object$plot.grouplabel <- grouplabel
  if (!is.na(filetype)) object$plot.filetype <- filetype
  if (!is.na(figsize)) object$plot.figsize <- figsize
  if (!is.na(figunit)) object$plot.figunit <- figunit
  
  if (flipaxes) {
    plotwidths <- c(2, 3, 2, 3)
    plotalign <- "h"
    decreasingLoadings <- !decreasingLoadings
  } else {
    plotwidths <- c(1, 2, 1, 2)
    plotalign <- "hv"
  }
  if (only == "score") {
    if (effect == "both") {
      g_score_time <- getScorePlot(object, component = component, effect = "time", plotribbon = plotribbon, dodgewidth = dodgewidth, myTheme = myTheme)
      g_score_group <- getScorePlot(object, component = component, effect = "group", plotribbon = plotribbon, dodgewidth = dodgewidth, myTheme = myTheme)
      g <- list(g_score_time, g_score_group)
    } else {
      g <- getScorePlot(object, component = component, effect = effect, plotribbon = plotribbon, dodgewidth = dodgewidth, myTheme = myTheme)
    }
  } else if (only == "loading") {
    if (effect == "both") {
      g_loading_time <- getLoadingPlot(object,
        component = component,
        effect = "time",
        decreasingLoadings = decreasingLoadings,
        flipaxes = flipaxes,
        plotzeroline = plotzeroline,
        tooDense = tooDense,
        n.limit = n.limit,
        highlight = highlight,
        variables = variables,
        loadinggroup = loadinggroup,
        limit_loading = limit_loading,
        sortbyloadinggroup = sortbyloadinggroup,
        myTheme = myTheme
      )
      g_loading_group <- getLoadingPlot(object,
        component = component,
        effect = "group",
        decreasingLoadings = decreasingLoadings,
        flipaxes = flipaxes,
        plotzeroline = plotzeroline,
        tooDense = tooDense,
        n.limit = n.limit,
        variables = variables,
        highlight = highlight,
        limit_loading = limit_loading,
        loadinggroup = loadinggroup,
        sortbyloadinggroup = sortbyloadinggroup,
        myTheme = myTheme
      )
      g <- list(g_loading_time, g_loading_group)
    } else {
      g <- getLoadingPlot(object,
        component = component,
        effect = effect,
        decreasingLoadings = decreasingLoadings,
        flipaxes = flipaxes,
        plotzeroline = plotzeroline,
        tooDense = tooDense,
        n.limit = n.limit,
        variables = variables,
        limit_loading = limit_loading,
        loadinggroup = loadinggroup,
        sortbyloadinggroup = sortbyloadinggroup,
        myTheme = myTheme
      )
    }
  } else {
    if (effect == "both") {
      g_loading_time <- getLoadingPlot(object,
        component = component,
        effect = "time",
        decreasingLoadings = decreasingLoadings,
        flipaxes = flipaxes,
        plotzeroline = plotzeroline,
        tooDense = tooDense,
        n.limit = n.limit,
        limit_loading = limit_loading,
        highlight = highlight,
        loadinggroup = loadinggroup,
        sortbyloadinggroup = sortbyloadinggroup,
        myTheme = myTheme
      )
      g_loading_group <- getLoadingPlot(object,
        component = component,
        effect = "group",
        decreasingLoadings = decreasingLoadings,
        flipaxes = flipaxes,
        plotzeroline = plotzeroline,
        tooDense = tooDense,
        n.limit = n.limit,
        variables = variables,
        highlight = highlight,
        limit_loading = limit_loading,
        loadinggroup = loadinggroup,
        sortbyloadinggroup = sortbyloadinggroup,
        myTheme = myTheme
      )
      g_score_time <- getScorePlot(object, component = component, effect = "time", plotribbon = plotribbon, dodgewidth = dodgewidth, myTheme = myTheme)
      g_score_group <- getScorePlot(object, component = component, effect = "group", plotribbon = plotribbon, dodgewidth = dodgewidth, myTheme = myTheme)
      if (enlist) {
        g <- list(g_score_time, g_loading_time, g_score_group, g_loading_group)
      } else {
        if (is.na(loadinggroup) & is.na(object$plot.loadinggroupcolumn)) {
          g <- ggpubr::ggarrange(g_score_time,
            g_loading_time,
            g_score_group,
            g_loading_group,
            nrow = 2, ncol = 2,
            widths = plotwidths,
            align = plotalign,
            common.legend = TRUE,
            legend.grob = ggpubr::get_legend(g_score_group),
            legend = "bottom"
          )
        } else {
          g <- ggpubr::ggarrange(g_score_time,
            g_loading_time,
            g_score_group,
            g_loading_group,
            nrow = 2, ncol = 2,
            widths = plotwidths,
            align = plotalign
          )
        }
      }
    } else {
      g_loading <- getLoadingPlot(object,
        component = component,
        effect = effect,
        decreasingLoadings = decreasingLoadings,
        flipaxes = flipaxes,
        plotzeroline = plotzeroline,
        tooDense = tooDense,
        n.limit = n.limit,
        variables = variables,
        limit_loading = limit_loading,
        highlight = highlight,
        loadinggroup = loadinggroup,
        sortbyloadinggroup = sortbyloadinggroup,
        myTheme = myTheme
      )
      g_score <- getScorePlot(object, component = component, effect = effect, plotribbon = plotribbon, dodgewidth = dodgewidth, myTheme = myTheme)
      if (enlist) {
        g <- list(g_score, g_loading)
      } else {
        if (is.na(loadinggroup) & is.na(object$plot.loadinggroupcolumn)) {
          g <- ggpubr::ggarrange(g_score,
            g_loading,
            nrow = 1,
            ncol = 2,
            widths = plotwidths,
            align = plotalign,
            common.legend = TRUE,
            legend = "bottom"
          )
        } else {
          g <- ggpubr::ggarrange(g_score,
            g_loading,
            nrow = 1,
            ncol = 2,
            widths = plotwidths,
            align = plotalign
          )
        }
      }
    }
  }
  if (object$save) {
    saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  }
  return(g)
}

#' Get screeplot
#'
#' This function returns a screeplot for an ALASCA model showing what proportion of the variance each component of the model explains
#'
#' @param object An ALASCA object
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(12, 8, 300)`)
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return An ggplot2 object (or a list og ggplot objects)
#'
#' @examples
#' load("PE.Rdata")
#' screeplot(model)
#' @export
screeplot.ALASCA <- function(object,
                             effect = "both",
                             nComps = NA,
                             filename = "scree_plot",
                             filetype = object$plot.filetype,
                             figsize = object$plot.figsize,
                             figunit = object$plot.figunit,
                             myTheme = object$plot.myTheme) {
  explained <- as.data.frame(get_scores(object)$explained)
  explained$component <- seq_len(nrow(explained))
  if (!is.na(filename)) object$filename <- filename
  
  if (any(!is.na(nComps))) {
    if (length(nComps) == 1) {
      explained <- subset(explained, component <= nComps)
    } else {
      explained <- subset(explained, component %in% nComps)
      explained$component <- factor(explained$component)
    }
  }
  g <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = time, group = NA)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    myTheme +
    ggplot2::labs(x = "Principal Component", y = paste0("Relative Expl. of ", object$plot.xlabel, " Var."))
  if (object$separateTimeAndGroup) {
    gg <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = group, group = NA)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      myTheme +
      ggplot2::labs(x = "Principal Component", y = "Relative Expl. of Group Var.")
    if (effect == "both") {
      g <- ggpubr::ggarrange(g, gg)
    }
  }
  if (effect == "group") {
    if (object$save) {
      saveALASCAPlot(object = object, g = gg, filetype = filetype, figsize = figsize, figunit = figunit)
    }
    return(gg)
  } else {
    if (object$save) {
      saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
    }
    return(g)
  }
}

#' Get loadings
#'
#' This function  returns the loadings for an ALASCA model
#'
#' @param object An ALASCA object
#' @param n.limit Returns the n highest and lowest loadings by PC (i.e., 2*n.limit loadings per PC)
#' @return A list with loadings for time (and group), and the exploratory power for each component
#' @export
get_loadings <- function(object, limit_loading = FALSE, n.limit = 0L, component = c(0)) {
  dfl <- list()
  if (component[[1]] > 0 || length(component) > 1) {
    dfl$time <- object$ALASCA$loading$time[PC %in% component]
    if (object$separateTimeAndGroup) {
      dfl$group <- object$ALASCA$loading$group[PC %in% component]
    }
  } else {
    dfl$time <- object$ALASCA$loading$time
    if (object$separateTimeAndGroup) {
      dfl$group <- object$ALASCA$loading$group
    }
  }
  if (limit_loading && object$validate) {
    dfl$time <- dfl$time[!is.na(low) & sign(low) == sign(high)]
    if (object$separateTimeAndGroup) {
      dfl$group <- dfl$group[!is.na(low) & sign(low) == sign(high)]
    }
  } else {
    if (n.limit > 0L) {
      index_head_and_tail <- c(seq(n.limit), length(object$variablelist)+1-seq(n.limit))
      dfl$time <- dfl$time[dfl$time[order(loading, decreasing = TRUE), .I[index_head_and_tail], by = PC]$V1]
      if (object$separateTimeAndGroup) {
        dfl$group <- dfl$group[dfl$group[order(loading, decreasing = TRUE), .I[index_head_and_tail], by = PC]$V1]
      }
    }
  }

  return(dfl)
}

#' Get scores
#'
#' This function returns the scores for an ALASCA model
#'
#' @inheritParams get_loadings
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
get_scores <- function(object, component = 0) {
  if(any(component > 0)){
    object$ALASCA$score$time <- object$ALASCA$score$time[ PC %in% component]
    if (object$separateTimeAndGroup) {
      object$ALASCA$score$group <- object$ALASCA$score$group[ PC %in% component]
    }
  }
  return(object$ALASCA$score)
}

#' Get covariables
#'
#' This function returns the other covariables in an ALASCA model
#'
#' @inheritParams get_loadings
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
get_covars <- function(object, n.limit = 0) {
  if (n.limit > 0) {
    return(
      rbind(
        object$covar_coefficients[order(estimate, decreasing = TRUE), head(.SD, n.limit), by = variable],
        object$covar_coefficients[order(estimate, decreasing = FALSE), head(.SD, n.limit), by = variable]
      )
    )
  } else {
    return(object$covar_coefficients)
  }
}

#' Get loading plot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An ALASCA object
#' @param component Which component to plot?
#' @param effect Plot time or group
#' @param decreasingLoadings Logical. Should loadings be sorted in decreasing order?
#' @param flipaxes When `TRUE` (default), list the variable loadings vertical instead of horizontal
#' @param plotzeroline When `TRUE` (default), plot a zero line in the loading plot
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(12, 8, 300)`)
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot object
getLoadingPlot <- function(object,
                           component = 1,
                           effect = "time",
                           decreasingLoadings = TRUE,
                           tooDense = NA,
                           n.limit = 0,
                           highlight = NA,
                           filetype = NA,
                           flipaxes = TRUE,
                           plotzeroline = TRUE,
                           figsize = NA,
                           figunit = NA,
                           variables = NA,
                           loadinggroup = NA,
                           limit_loading = FALSE,
                           sortbyloadinggroup = TRUE,
                           point_size = 0.4,
                           myTheme = NA) {
  
  if (any(is.na(myTheme))) myTheme <- object$plot.myTheme
  if (!is.na(filetype)) object$plot.filetype <- filetype
  if (any(is.na(variables))) variables <- object$variablelist
  if (!is.na(figsize)) object$plot.figsize <- figsize
  if (!is.na(figunit)) object$plot.figunit <- figunit
  
  if (effect == "time") {
    loadings <- subset(get_loadings(object, limit_loading = limit_loading, n.limit = n.limit, component = component)$time, covars %in% variables)
  } else {
    loadings <- subset(get_loadings(object, limit_loading = limit_loading, n.limit = n.limit, component = component)$group, covars %in% variables)
  }
  if (!is.na(object$plot.loadinggroupcolumn)) {
    loadings <- merge(loadings, object$variable_labels)
  } else {
    loadings$covargroup <- NA
  }
  if (sortbyloadinggroup & !is.na(object$plot.loadinggroupcolumn)) {
    loadings$covars <- factor(loadings$covars, levels = unique(loadings$covars[order(loadings$covargroup, loadings$loading, decreasing = decreasingLoadings)]))
  } else {
    loadings$covars <- factor(loadings$covars, levels = unique(loadings$covars[order(loadings$loading, decreasing = decreasingLoadings)]))
  }
  if (object$validate) {
    if (any(colnames(loadings) == "model")) {
      g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high, shape = model)) +
        ggplot2::geom_pointrange(size = point_size)
    } else {
      if (is.na(object$plot.loadinggroupcolumn)) {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high)) +
          ggplot2::geom_pointrange(size = point_size)
      } else {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high, color = covargroup, shape = covargroup)) +
          ggplot2::geom_pointrange(size = point_size)
      }
    }
  } else {
    if (any(colnames(loadings) == "model")) {
      g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, shape = model)) +
        ggplot2::geom_point()
    } else {
      if (is.na(object$plot.loadinggroupcolumn)) {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading)) +
          ggplot2::geom_point()
      } else {
        g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, color = covargroup, shape = covargroup)) +
          ggplot2::geom_point()
      }
    }
  }
  g <- g + myTheme +
    ggplot2::labs(
      x = "Variable",
      y = .getExpLabel(object, component = component, effect = effect, type = "Loading")
    )
  if (plotzeroline) {
    g <- g + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  }
  if (n.limit > 0 & !sortbyloadinggroup) {
    g <- g + ggplot2::geom_vline(xintercept = n.limit + 0.5, linetype = "dotted")
  }
  if (flipaxes) {
    g <- g + ggplot2::coord_flip()
  } else {
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = c(0.8, 0.8))
  }
  if (!is.na(object$plot.loadinggroupcolumn)) {
    g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
      ggplot2::labs(color = object$plot.loadinggrouplabel, shape = object$plot.loadinggrouplabel) +
      ggplot2::theme(legend.position = "bottom") # ggplot2::scale_color_brewer(palette = "Dark2")
  }
  if (!any(is.na(highlight))) {
    g <- g + ggplot2::geom_point(color = ifelse(loadings$covars %in% highlight, "red", "grey")) +
      ggrepel::geom_text_repel(data = subset(loadings, covars %in% highlight), ggplot2::aes(label = covars), max.iter = 5000) +
      myTheme +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )
  } else if (!is.na(tooDense) & tooDense > 0) {
    limUpper <- unique(loadings$loading[order(loadings$loading, decreasing = TRUE)[tooDense]])
    limLower <- unique(loadings$loading[order(loadings$loading, decreasing = FALSE)[tooDense]])
    g <- g + ggplot2::geom_point(color = ifelse(loadings$loading <= limLower | loadings$loading >= limUpper, "red", "grey")) +
      ggrepel::geom_text_repel(
        data = subset(loadings, loading <= limLower | loading >= limUpper),
        ggplot2::aes(label = covars),
        max.iter = 5000
      ) +
      myTheme
    if (flipaxes) {
      g <- g + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      )
    } else {
      g <- g + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )
    }
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit, suffix = "_loading")

  return(g)
}

#' Get score plot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An ALASCA object
#' @param component Which component to plot?
#' @param effect Plot time or group
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::scale_color_viridis_d(end = 0.9) + ggplot2::theme_bw()`
#' @return A ggplot object
getScorePlot <- function(object,
                         component = 1,
                         effect = "time",
                         filetype = NA,
                         figsize = NA,
                         figunit = NA,
                         plotribbon = TRUE,
                         dodgewidth = 0.35,
                         point_size = 2,
                         myTheme = NA) {
  if (any(is.na(myTheme))) myTheme <- object$plot.myTheme
  if (!is.na(filetype)) object$plot.filetype <- filetype
  if (any(!is.na(figsize))) object$plot.figsize <- figsize
  if (!is.na(figunit)) object$plot.figunit <- figunit

  if (effect == "time") {
    if (object$separateTimeAndGroup) {
      score <- get_scores(object, component = component)$time
      if (object$validate) {
        # Show error bars
        if (grepl("permutation", object$validationMethod)) {
          pvals <- object$pvals
          score <- merge(score, pvals, by.x = "time", by.y = "effect", all.x = TRUE, all.y = FALSE)
          score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, label = p.value.str)) +
            ggplot2::geom_point(position = ggplot2::position_dodge(width = dodgewidth)) +
            ggplot2::geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = dodgewidth), show.legend = FALSE) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
        } else {
          if (any(colnames(score) == "model")) {
            g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = model, linetype = model, shape = model, ymin = low, ymax = high)) +
              ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
              ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
          } else {
            g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, ymin = low, ymax = high)) +
              ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth))
            if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
            if (plotribbon && object$method %in% c("LMM")) {
              g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                alpha = .1,
                position = ggplot2::position_dodge(width = dodgewidth), color = NA
              ) +
                ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
                ggplot2::labs(fill = object$plot.grouplabel)
            }
          }
        }
      } else {
        # No validation - no error bars
        if (any(colnames(score) == "model")) {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = model, linetype = model)) +
            ggplot2::geom_point() +
            ggplot2::geom_line()
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group)) +
            ggplot2::geom_point()
          if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line()
        }
      }
    } else {
      score <- get_scores(object, component = component)$time
      if (object$validate) {
        if (grepl("permutation", object$validationMethod)) {
          pvals <- object$pvals
          score$effect <- paste(score$time, score$group)
          score <- merge(score, pvals, by.x = "effect", by.y = "effect", all.x = TRUE, all.y = FALSE)
          score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, label = p.value.str)) +
            ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
            geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth))
            if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
          if (plotribbon && object$method %in% c("LMM")) {
            g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
              alpha = .1,
              position = ggplot2::position_dodge(width = dodgewidth), color = NA
            ) +
              ggplot2::scale_fill_manual(values = getPlotPalette(object)) + ggplot2::labs(fill = object$plot.grouplabel)
          }
        }
      } else {
        if (any(colnames(score) == "model")) {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = model)) +
            ggplot2::geom_point() +
            ggplot2::geom_line()
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group)) +
            ggplot2::geom_point()
          if (object$method %in% c("LMM")) g <- g + ggplot2::geom_line()
        }
      }
    }
    g <- g + myTheme +
      ggplot2::scale_color_manual(values = getPlotPalette(object))
    if (object$method %in% c("LMM")) g <- g + ggplot2::scale_linetype_manual(values = getPlotLinetypes(object))
    g <- g + ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        x = object$plot.xlabel,
        group = object$plot.grouplabel, color = object$plot.grouplabel, linetype = object$plot.grouplabel,
        y = .getExpLabel(object, component = component, effect = "time")
      )
  } else {
    # Group effect
    score <- get_scores(object, component = component)$group
    if (object$validate) {
      if (grepl("permutation", object$validationMethod)) {
        pvals <- object$pvals
        score$effect <- paste(score$time, score$group)
        score <- merge(score, pvals, by.x = "effect", by.y = "effect", all.x = TRUE, all.y = FALSE)
        score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, label = p.value.str)) +
          ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
          geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
      } else {
        if (any(colnames(score) == "model")) {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, linetype = model, shape = model, color = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
        } else {
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth))
          if (plotribbon) {
            g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
              alpha = .1,
              position = ggplot2::position_dodge(width = dodgewidth), color = NA
            ) +
              ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
              ggplot2::labs(fill = object$plot.grouplabel)
          }
        }
      }
    } else {
      if (any(colnames(score) == "model")) {
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = model)) +
          ggplot2::geom_point() +
          ggplot2::geom_line()
      } else {
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, linetype = group)) +
          ggplot2::geom_point() +
          ggplot2::geom_line()
      }
    }
    g <- g + myTheme +
      ggplot2::scale_color_manual(values = getPlotPalette(object)) +
      ggplot2::scale_linetype_manual(values = getPlotLinetypes(object)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        x = object$plot.xlabel,
        group = object$plot.grouplabel, color = object$plot.grouplabel, linetype = object$plot.grouplabel,
        y = .getExpLabel(object, component = component, effect = "group")
      )
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit, suffix = "_score")
  
  return(g)
}

#' Plot participants
#'
#' This function returns the scores for an ALASCA model
#'
#' @param object An ALASCA object or a data frame. If a data frame, you need to specify the column names for participant and value. This also applies if you have not specified the participant column in the ALASCA model before.
#' @param variable List of variable names to print. If `NA`, return all (default).
#' @param participantColumn Specify the column with participant identifier. Not necessary if you have already provided it to the ALASCA object
#' @param valueColumn Specify column with values (y axis). Not necessary to provide if you are plotting an ALASCA object.
#' @param timeColumn Specify column with times (x axis). Defaults to `time`.
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(width,height,dpi)` (default: `c(120, 80, 300)`)
#' @param addSmooth. Specify which geom_smooth model you want to apply, eg. `lm`, `glm`, `gam`, `loess` (default). Set to `NA` to remove.
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#'
#' @examples
#' load("PE.Rdata")
#' plotParts(df, variable = "IL-6", participantColumn = "ID", valueColumn = "value")
#' do.call(
#'   ggpubr::ggarrange,
#'   c(plotParts(df,
#'     participantColumn = "ID", timeColumn = "GA", valueColumn = "value", addSmooth = NA,
#'     variable = c("PlGF", "IL-6", "IL-1b", "IFN-g", "Eotaxin-2", "Eotaxin")
#'   ),
#'   common.legend = TRUE, legend = "bottom"
#'   )
#' )
#' plotParts(model, variable = "IL-6", participantColumn = "ID", timeColumn = "GA")[[1]] + ggplot2::labs(x = "GA")
#' @export
plotParts <- function(object,
                      variables = NA,
                      participantColumn = FALSE,
                      xlabel = NA,
                      grouplabel = NA,
                      valueColumn = FALSE,
                      timeColumn = "time",
                      addSmooth = "loess",
                      filename = NA,
                      filetype = object$plot.filetype,
                      figunit = object$plot.figunit,
                      plot.ylabel = "value",
                      as.list = FALSE,
                      figsize = object$plot.figsize,
                      myTheme = object$plot.myTheme) {
  
  if (!is.na(filename)) object$filename <- filename
  
  if(as.list){
    if (is.data.frame(object)) {
      df <- object
      if (any(participantColumn == FALSE) | any(valueColumn == FALSE)) {
        stop("You need to specify participant and value columns")
      } else {
        participantColumn <- participantColumn
        valueColumn <- valueColumn
      }
      plotFunction <- function(df, timeColumn, valueColumn, participantColumn, xi, addSmooth, xlabel, grouplabel, myTheme) {
        g <- ggplot2::ggplot(subset(df, variable == xi), ggplot2::aes_string(x = timeColumn, y = valueColumn, color = "group", group = participantColumn)) +
          ggplot2::geom_point(alpha = 0.7) +
          ggplot2::geom_line(alpha = 0.3) +
          ggplot2::scale_color_manual(values = getPlotPalette(list(df = df))) +
          ggplot2::scale_fill_manual(values = getPlotPalette(list(df = df))) +
          myTheme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = xlabel, y = xi, color = grouplabel, fill = grouplabel)
        if (!any(is.na(addSmooth))) {
          g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group, fill = group), se = TRUE)
        }
        return(g)
      }
    } else if (is(object, "ALASCA")) {
      df <- object$df_raw
      valueColumn <- as.character(object$formula)[2]
      if (any(participantColumn == FALSE)) {
        if (any(object$participantColumn == FALSE)) {
          stop("You need to specify the participant column")
        } else {
          participantColumn <- object$participantColumn
        }
      }
      
      if (is.na(xlabel)) xlabel <- object$plot.xlabel
      if (is.na(grouplabel)) grouplabel <- object$plot.grouplabel
      
      plotFunction <- function(df, timeColumn, valueColumn, participantColumn, xi, addSmooth, xlabel, grouplabel, myTheme) {
        g <- ggplot2::ggplot(subset(df, variable == xi), ggplot2::aes_string(x = timeColumn, y = valueColumn, color = "group", group = participantColumn)) +
          ggplot2::geom_point(alpha = 0.7) +
          ggplot2::geom_line(alpha = 0.3) +
          ggplot2::scale_color_manual(values = getPlotPalette(object)) +
          ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
          myTheme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = xlabel, y = xi, color = grouplabel, fill = grouplabel)
        if (!any(is.na(addSmooth))) {
          g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group, fill = group), se = TRUE)
        }
        return(g)
      }
    } else {
      stop("Wrong input object: must be a ALASCA model or a data frame")
    }
  
    if (any(is.na(variables))) variables <- unique(df$variable)

    g <- lapply(variables, function(xi) {
      plotFunction(df, timeColumn, valueColumn, participantColumn, xi, addSmooth, xlabel = xlabel, grouplabel = grouplabel, myTheme = myTheme)
    })
    names(g) <- variables
    if (is(object, "ALASCA")) {
      if (object$save) {
        for (i in seq_along(g)) {
          saveALASCAPlot(object = object, g = g[[i]], filetype = filetype, figsize = figsize, figunit = figunit, suffix = names(g)[i])
        }
      }
    }
    return(g)
  } else {
      if (any(is.na(variables))) {
        variables <- object$variablelist
      }
      df <- object$df_raw
      if (is.na(xlabel)) {
        xlabel <- object$plot.xlabel
      }
      if (is.na(grouplabel)) {
        grouplabel <- object$plot.grouplabel
      }
      g <- ggplot2::ggplot(subset(df, variable %in% variables), ggplot2::aes(x = time, y = value, color = group, group = ID)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::geom_line(alpha = 0.3) +
        ggplot2::scale_color_manual(values = getPlotPalette(object)) +
        ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
        myTheme +
        ggplot2::facet_wrap(~variable, scales = "free_y") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = xlabel, y = plot.ylabel, color = grouplabel, fill = grouplabel)
      if (!any(is.na(addSmooth))) {
        g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group, fill = group), se = TRUE)
      }

    if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
    
    return(g)
  }
}

#' Plot model predictions
#'
#' This function returns the scores for an ALASCA model
#'
#' @param object An ALASCA object or a data frame. If a data frame, you need to specify the column names for participant and value. This also applies if you have not specified the participant column in the ALASCA model before.
#' @param variable List of variable names to print. If `NA`, return all (default).
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm",
#' @return A list with ggplot2 objects.
#'
#' @examples
#' load("PE.Rdata")
#' plotPred(model, variable = "IL-6")[[1]] + ggplot2::theme_bw()
#' do.call(
#'   ggpubr::ggarrange,
#'   c(plotPred(model, variable = c("PlGF", "IL-6", "IL-1b", "IFN-g", "Eotaxin-2", "Eotaxin")),
#'     common.legend = TRUE, legend = "bottom"
#'   )
#' )
#' @export
plotPred <- function(object,
                     variables = object$variablelist,
                     filetype = NA,
                     figsize = NA,
                     figunit = NA,
                     filename = NA,
                     dodgewidth = 0.35,
                     plotribbon = TRUE,
                     as.list = FALSE,
                     plot.ylabel = "value",
                     myTheme = object$plot.myTheme) {
  
  if (!is.na(filename)) object$filename <- filename
  
  if ( as.list ){
    if (object$validate) {
      gg <- lapply(variables, function(x) {
        g <- ggplot2::ggplot(subset(object$model_prediction, variable == x), ggplot2::aes(
          x = time,
          y = pred,
          color = group,
          group = group,
          linetype = group,
          ymin = low,
          ymax = high
        )) +
          ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth)) +
          ggplot2::scale_color_manual(values = getPlotPalette(object)) +
          myTheme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = object$plot.xlabel, y = plot.ylabel, color = object$plot.grouplabel, linetype = object$plot.grouplabel)
        if (plotribbon) {
          g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
            alpha = .1,
            position = ggplot2::position_dodge(width = dodgewidth), color = NA
          ) +
            ggplot2::scale_fill_manual(values = getPlotPalette(object)) + ggplot2::labs(fill = object$plot.grouplabel)
        }
        g
      })
    } else {
      gg <- lapply(variables, function(x) {
        g <- ggplot2::ggplot(subset(object$model_prediction, variable == x), ggplot2::aes(x = time, y = pred, color = group, linetype = group, group = group)) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_color_manual(values = getPlotPalette(object)) +
          myTheme +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(x = object$plot.xlabel, y = plot.ylabel, color = object$plot.grouplabel, linetype = object$plot.grouplabel)
        g
      })
    }
    names(gg) <- variable
    if (object$save) {
      for (i in seq_along(gg)) {
        saveALASCAPlot(object = object, g = gg[[i]], filetype = filetype, figsize = figsize, figunit = figunit, suffix = names(gg)[i])
      }
    }
  } else {
    if (object$validate) {
    g <- ggplot2::ggplot(object$model_prediction[object$model_prediction$variable %in% variables,], ggplot2::aes(
      x = time,
      y = pred,
      color = group,
      group = group,
      linetype = group,
      ymin = low,
      ymax = high
    )) +
      ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = dodgewidth)) +
      ggplot2::geom_line(position = ggplot2::position_dodge(width = dodgewidth)) +
      ggplot2::scale_color_manual(values = getPlotPalette(object)) +
      myTheme +
      ggplot2::facet_wrap(~variable, scales = "free_y") +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = object$plot.xlabel, y = plot.ylabel, color = object$plot.grouplabel, linetype = object$plot.grouplabel)
    if (plotribbon) {
      g <- g + ggplot2::geom_ribbon(ggplot2::aes(fill = group),
                                    alpha = .1,
                                    position = ggplot2::position_dodge(width = dodgewidth), color = NA
      ) +
        ggplot2::scale_fill_manual(values = getPlotPalette(object)) + ggplot2::labs(fill = object$plot.grouplabel)
    }
    } else {
      g <- ggplot2::ggplot(object$model_prediction[object$model_prediction$variable %in% variables,], ggplot2::aes(x = time, y = pred, color = group, linetype = group, group = group)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(values = getPlotPalette(object)) +
        myTheme +
        ggplot2::facet_wrap(~variable, scales = "free_y") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = object$plot.xlabel, y = plot.ylabel, color = object$plot.grouplabel, linetype = object$plot.grouplabel)
    }
    if (object$save) {
        saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
    }
    return(g)
  }
  
}

#' Plot validations models
#'
#' This function returns a plot of the validation models
#'
#' @param object A validated ALASCA object
#' @param component Which component to plot (default: 1)
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#'
#' @export
plotVal <- function(object,
                    component = 1,
                    filetype = NA,
                    filename = NA,
                    figsize = NA,
                    figunit = NA,
                    n.limit = 0,
                    decreasingLoadings = FALSE,
                    plotzeroline = TRUE,
                    myTheme = NA,
                    flipaxes = TRUE,
                    plot.alpha = 0.3) {
  if (!object$validate) stop("You must validate the model first")
  if (any(is.na(myTheme)))  myTheme <- object$plot.myTheme
  if (!is.na(filename)) object$filename <- filename
  
  if (object$separateTimeAndGroup) {
    # Score plots
    ## Time
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    }))
    dfm <- get_scores(object, component = component)$time
    gst <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = model, color = group, linetype = group)) +
      ggplot2::geom_point(alpha = plot.alpha)
    if (object$method %in% c("LMM")) gst <- gst + ggplot2::geom_line(alpha = plot.alpha)
    gst <- gst + ggplot2::geom_point(data = dfm, group = NA, alpha = 1, color = "black")
    if (object$method %in% c("LMM")) gst <- gst + ggplot2::geom_line(data = dfm, group = dfm$group, alpha = 1, color = "black")
    gst <- gst + ggplot2::labs(x = object$plot.xlabel,
                               color = object$plot.grouplabel, linetype = object$plot.grouplabel,
                               y = .getExpLabel(object, component = component, effect = "time", type = "Score")) +
      ggplot2::scale_color_manual(values = getPlotPalette(object)) +
      ggplot2::scale_linetype_manual(values = getPlotLinetypes(object)) +
      myTheme + ggplot2::theme(legend.position = "none")
    
    ## Group
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$group, PC == component),
        model = x
      )
    }))
    dfm <- get_scores(object, component = component)$group
    dff$plotGroup <- paste0(dff$model, "-", dff$group)
    gsg <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = plotGroup, color = group, linetype = group)) +
      ggplot2::geom_point(alpha = plot.alpha)
    if (object$method %in% c("LMM")) gsg <- gsg + ggplot2::geom_line(alpha = plot.alpha)
    gsg <- gsg + ggplot2::geom_point(data = dfm, group = NA, alpha = 1, color = "black")
    if (object$method %in% c("LMM")) gsg <- gsg + ggplot2::geom_line(data = dfm, group = dfm$group, alpha = 1, color = "black")
    gsg <- gsg +   ggplot2::scale_color_manual(values = getPlotPalette(object)) +
      ggplot2::scale_linetype_manual(values = getPlotLinetypes(object)) +
      ggplot2::labs(x = object$plot.xlabel,
                    y = .getExpLabel(object, component = component, effect = "group", type = "Score")) +
      myTheme + ggplot2::theme(legend.position = "bottom")

    # Loading plot
    ## Time
    dfm <- get_loadings(object, component = component, n.limit = n.limit)$time
    dfm$covars <- factor(dfm$covars, levels = unique(dfm$covars[order(dfm$loading, decreasing = decreasingLoadings)]))
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_loadings(object$validation$temp_objects[[x]], component = component)$time, covars %in% dfm$covars)
      )
    }))
    dff$covars <- factor(dff$covars, levels = unique(dff$covars[order(dff$loading, decreasing = decreasingLoadings)]))
    
    glt <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) +
      ggplot2::geom_point(data = dfm, alpha = 1, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(ymax = high, ymin = low), alpha = 0.3, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(xmax = as.numeric(covars) + 0.5, xmin = as.numeric(covars) - 0.5), alpha = 1, color = "black") +
      ggplot2::geom_point(alpha = 0.2, color = "black")
    if (plotzeroline) glt <- glt + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    if (flipaxes) glt <- glt + ggplot2::coord_flip()
      glt <- glt + ggplot2::labs(x = "Variable",
                    y = .getExpLabel(object, component = component, effect = "time", type = "Loading")) +
      ggplot2::theme(legend.position = "none") +
      myTheme

    ## Group
    dfm <- get_loadings(object, component = component, n.limit = n.limit)$group
    dfm$covars <- factor(dfm$covars, levels = unique(dfm$covars[order(dfm$loading, decreasing = decreasingLoadings)]))
    dff <- rbindlist(lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_loadings(object$validation$temp_objects[[x]], component = component)$group, covars %in% dfm$covars)
      )
    }))
    dff$covars <- factor(dff$covars, levels = unique(dff$covars[order(dff$loading, decreasing = decreasingLoadings)]))
    glg <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) +
      ggplot2::geom_point(data = dfm, alpha = 1, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(ymax = high, ymin = low), alpha = 0.3, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(xmax = as.numeric(covars) + 0.5, xmin = as.numeric(covars) - 0.5), alpha = 1, color = "black") +
      ggplot2::geom_point(alpha = 0.2, color = "black")
    if (plotzeroline) glg <- glg + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    if (flipaxes) glg <- glg + ggplot2::coord_flip()
    glg <- glg + ggplot2::labs(x = "Variable",
                    y = .getExpLabel(object, component = component, effect = "group", type = "Loading")) +
      ggplot2::theme(legend.position = "bottom") +
      myTheme

    g <- ggpubr::ggarrange(gst, glt, gsg, glg, nrow = 2, ncol = 2, widths = c(2, 3, 2, 3), common.legend = TRUE, legend = "bottom")
  } else {
    # Score plot
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    }))
    dff$plotGroup <- paste0(dff$model, "-", dff$group)
    dfm <- get_scores(object, component = component)$time
    gs <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = plotGroup, color = group, linetype = group)) +
      ggplot2::geom_point(alpha = plot.alpha)
      if (object$method %in% c("LMM")) gs <- gs + ggplot2::geom_line(alpha = plot.alpha)
      gs <- gs + ggplot2::geom_point(data = dfm, group = NA, alpha = 1, color = "black")
      if (object$method %in% c("LMM")) gs <- gs + ggplot2::geom_line(data = dfm, group = dfm$group, alpha = 1, color = "black")
      gs <- gs + ggplot2::scale_color_manual(values = getPlotPalette(object)) +
      ggplot2::scale_linetype_manual(values = getPlotLinetypes(object)) +
      ggplot2::labs(x = object$plot.xlabel,
                    color = object$plot.grouplabel, linetype = object$plot.grouplabel,
                    y = .getExpLabel(object, component = component, effect = "time", type = "Score")) +
      myTheme +
      ggplot2::theme(legend.position = "bottom")

    # Loading plot
    dfm <- get_loadings(object, component = component, n.limit = n.limit)$time
    dfm$covars <- factor(dfm$covars, levels = unique(dfm$covars[order(dfm$loading, decreasing = decreasingLoadings)]))
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_loadings(object$validation$temp_objects[[x]], component = component)$time, covars %in% dfm$covars)
      )
    }))
    dff$covars <- factor(dff$covars, levels = unique(dff$covars[order(dff$loading, decreasing = decreasingLoadings)]))
    
    gl <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) +
      ggplot2::geom_point(data = dfm, alpha = 1, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(ymax = high, ymin = low), alpha = 0.3, color = "black") +
      ggplot2::geom_linerange(data = dfm, ggplot2::aes(xmax = as.numeric(covars) + 0.5, xmin = as.numeric(covars) - 0.5), alpha = 1, color = "black") +
      ggplot2::geom_point(alpha = 0.2, color = "black")
    if (plotzeroline) gl <- gl + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    if (flipaxes) gl <- gl + ggplot2::coord_flip()
    gl <- gl + ggplot2::labs(x = "Variable",
                    color = object$plot.grouplabel,
                    y = .getExpLabel(object, component = component, effect = "time", type = "Loading")) +
      myTheme

    g <- ggpubr::ggarrange(gs, gl, nrow = 1, ncol = 2, widths = c(2, 3), common.legend = TRUE, legend = "bottom")
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  
  return(g)
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param xlabel Alternative names for the covariables
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm",
#' @param return_data Set to `TRUE` to return data instead of plot
#'
#' @param myTheme A ggplot2 theme
#' @return A ggplot2 objects\.
#'
#' @export
plotCovar <- function(object,
                      covar = NA,
                      xlabel = NA,
                      variables = NA,
                      n.limit = 0,
                      filename = "covars",
                      return_data = FALSE,
                      filetype = NA,
                      figsize = NA,
                      figunit = NA,
                      myTheme = NA,
                      pvalue = "star") {
  if (any(is.na(myTheme))) myTheme <- object$plot.myTheme
  if (!is.na(filename)) object$filename <- filename
  
  df <- get_covars(object, n.limit = n.limit)
  if ( nrow(df) == 0 ) stop("No covariates to plot")
  df$covar <- factor(df$covar, levels = unique(df$covar[order(df$estimate)]))
  if (any(is.na(covar))) {
    covar <- unique(df$variable)
  }
  if (any(is.na(variables))) {
    variables <- object$variablelist
  }
  df <- subset(df, covar %in% variables)
  if (any(is.na(xlabel))) {
    xlabel <- covar
  }
  for (i in seq_len(length(covar))) {
    df$xlabel[df$variable == covar[i]] <- xlabel[i]
  }
  if (!is.na(object$plot.loadinggroupcolumn)) {
    loadings <- merge(loadings, object$variable_labels, by.x = "covar", by.y = "covars")
  } else {
    df$covargroup <- NA
  }
  if (!object$useRfast) {
    # lmer and lm provide p values
    df$pvalue_label <- ifelse(df$pvalue >= 0.05, "Not significant", ifelse(df$pvalue < 0.001, "< 0.001", ifelse(df$pvalue < 0.01, "< 0.01", "< 0.05")))
    df$pvalue_sign <- ifelse(df$pvalue >= 0.05, "", ifelse(df$pvalue < 0.001, "***", ifelse(df$pvalue < 0.01, "**", "*")))
    if (return_data) {
      return(df)
    } else {
      if (pvalue == "shape") {
        object$plot.loadinggroupcolumn <- NA # Cannot use shape for both loadinggroup and significance
        g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, shape = pvalue_label)) +
          ggplot2::scale_shape_manual(values = c("Not significant" = 3, "< 0.05" = 15, "< 0.01" = 16, "< 0.001" = 17, "Baseline" = 5))
      } else if (pvalue == "star" | pvalue == "asterisk" | pvalue == "stars") {
        if (is.na(object$plot.loadinggroupcolumn)) {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, label = pvalue_sign)) +
            ggplot2::geom_text(hjust = 0.5, vjust = 0, show.legend = FALSE)
        } else {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, label = pvalue_sign, color = covargroup, shape = covargroup)) +
            ggplot2::geom_text(hjust = 0.5, vjust = 0, show.legend = FALSE)
        }
      } else {
        g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar))
      }
      g <- g +
        ggplot2::geom_point() +
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::facet_wrap(~xlabel, scales = "free_y") + ggplot2::labs(x = "Coefficient", y = "", shape = "P value") +
        myTheme + ggplot2::theme(legend.position = "bottom", legend.box = "vertical", legend.margin = ggplot2::margin())

      if (!is.na(object$plot.loadinggroupcolumn)) {
        g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
          ggplot2::labs(color = object$plot.loadinggrouplabel, shape = object$plot.loadinggrouplabel) +
          ggplot2::theme(legend.position = "bottom") # ggplot2::scale_color_brewer(palette = "Dark2")
      }
      
      if (object$save) {
        saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
      }
      return(g)
    }
  } else {
    # Rfast does not provide p values, use bootstrap intervals
    df$covar <- factor(df$covar, levels = unique(df$covar[order(df$estimate)]))
    if (return_data) {
      return(df)
    } else {
      if (object$validate) {
        if (is.na(object$plot.loadinggroupcolumn)) {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, xmin = low, xmax = high)) +
            ggplot2::geom_pointrange()
        } else {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, xmin = low, xmax = high, color = covargroup, shape = covargroup)) +
            ggplot2::geom_pointrange()
        }
      } else {
        if (is.na(object$plot.loadinggroupcolumn)) {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar)) +
            ggplot2::geom_point()
        } else {
          g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, color = covargroup, shape = covargroup)) +
            ggplot2::geom_point()
        }
      }
      
      g <- g + 
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::facet_wrap(~xlabel, scales = "free_y") +
        ggplot2::labs(x = "Coefficient", y = "") +
        myTheme
      
      if (!is.na(object$plot.loadinggroupcolumn)) {
        g <- g + ggplot2::scale_color_viridis_d(option = "A", end = 0.85) +
          ggplot2::labs(color = object$plot.loadinggrouplabel, shape = object$plot.loadinggrouplabel) +
          ggplot2::theme(legend.position = "bottom") # ggplot2::scale_color_brewer(palette = "Dark2")
      }

      if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
      
      return(g)
    }
  }
}

#' Plot PCs
#'
#' This function returns a plot of ...
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param xlab Alternative names for the covariables
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm"
#' @param plottext If `TRUE`, plot time point as text
#' @param validationshape Either `NA`, "ellipse", or "cross"
#' @param validationlevel Defaults to 0.95
#' @param return_data Set to `TRUE` to return data instead of plot
#'
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot2 object.
#'
#' @export
plotComponents <- function(object,
                           comps = c(1, 2),
                           filename = NA,
                           filetype = NA,
                           figsize = NA,
                           figunit = NA,
                           ...) {
  g_score <- plotComponentsScore(object = object, comps = comps, ...)
  g_loading <- plotComponentsLoadings(object = object, comps = comps, ...)
  g <- ggpubr::ggarrange(
    g_score, g_loading,
    ncol = 1, nrow = 2
  )
  if (!is.na(filename)) object$filename <- filename
  if (object$save) {
    saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  }
  return(g)
}

plotComponentsLoadings <- function(object,
                                   comps = c(1, 2),
                                   filetype = NA,
                                   figsize = NA,
                                   figunit = NA,
                                   validationshape = NA,
                                   myTheme = ggplot2::theme_classic(),
                                   ...) {
  if (any(!comps %in% get_relevant_pcs(object = object, effect = "time"))) {
    warning("Please note: Some components have low explanatory power and HAVE NOT BEEN rotated during rotation. Proceed with care.")
  }

  dff <- subset(get_loadings(object)$time, PC %in% comps)
  dff$PC <- paste0("PC", dff$PC)
  dff <- reshape2::melt(data = dff, id.vars = c("PC", "covars"))
  dff <- reshape2::dcast(data = dff, covars ~ PC + variable, value.var = "value")
  dff$color <- "black"
  dff$lab <- dff$covars
  if (object$validate == TRUE) {
    dff$color <- ifelse(sign(dff[paste0("PC", comps[1], "_low")]) == sign(dff[paste0("PC", comps[1], "_high")]) & sign(dff[paste0("PC", comps[2], "_low")]) == sign(dff[paste0("PC", comps[2], "_high")]), "black", "gray")
    dff$lab <- ifelse(dff$color == "black", dff$covars, NA)
  }


  g <- ggplot2::ggplot(dff, ggplot2::aes_string(
    x = paste0("PC", comps[1], "_loading"),
    y = paste0("PC", comps[2], "_loading"),
    label = "lab"
  ))
  if (object$validate == TRUE) {
    g <- g +
      ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high")), color = dff$color) +
      ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")), color = dff$color)
  } else {
    g <- g + ggplot2::geom_point()
  }
  g <- g +
    ggrepel::geom_text_repel() +
    ggplot2::labs(
      x = .getExpLabel(object, component = comps[1], effect = "time", type = "Loading"),
      y = .getExpLabel(object, component = comps[2], effect = "time", type = "Loading")
    ) +
    myTheme

  if (object$separateTimeAndGroup) {
    dff <- subset(get_loadings(object)$group, PC %in% comps)
    dff$PC <- paste0("PC", dff$PC)
    dff <- reshape2::melt(data = dff, id.vars = c("PC", "covars"))
    dff <- reshape2::dcast(data = dff, covars ~ PC + variable, value.var = "value")
    dff$color <- "black"
    dff$lab <- dff$covars
    if (object$validate == TRUE) {
      dff$color <- ifelse(sign(dff[paste0("PC", comps[1], "_low")]) == sign(dff[paste0("PC", comps[1], "_high")]) & sign(dff[paste0("PC", comps[2], "_low")]) == sign(dff[paste0("PC", comps[2], "_high")]), "black", "gray")
      dff$lab <- ifelse(dff$color == "black", dff$covars, NA)
    }

    glg <- ggplot2::ggplot(dff, ggplot2::aes_string(
      x = paste0("PC", comps[1], "_loading"),
      y = paste0("PC", comps[2], "_loading"),
      label = "lab"
    ))
    if (object$validate == TRUE) {
      glg <- glg +
        ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high")), color = dff$color) +
        ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")), color = dff$color)
    } else {
      glg <- glg + ggplot2::geom_point()
    }
    glg <- glg +
      ggrepel::geom_text_repel() +
      ggplot2::labs(
        x = .getExpLabel(object, component = comps[1], effect = "group", type = "Loading"),
        y = .getExpLabel(object, component = comps[2], effect = "group", type = "Loading")
      ) +
      myTheme

    g <- ggpubr::ggarrange(g, glg)
  }
  
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  
  return(g)
}

#' Plot scores along PCs
#'
#' This function returns a plot of ...
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param xlab Alternative names for the covariables
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param figunit = "mm"
#' @param plottext If `TRUE`, plot time point as text
#' @param validationshape Either `NA`, "ellipse", or "cross"
#' @param validationlevel Defaults to 0.95
#' @param return_data Set to `TRUE` to return data instead of plot
#'
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot2 object.
#'
#' @export
plotComponentsScore <- function(object,
                                comps = c(1, 2),
                                xlabel = NA,
                                return_data = FALSE,
                                filename = NA,
                                filetype = NA,
                                figsize = NA,
                                figunit = NA,
                                validationshape = "ellipse",
                                validationlevel = 0.95,
                                plottext = TRUE,
                                texthjust = -0.2,
                                alphavalidate = 0.4,
                                myTheme = ggplot2::theme_classic(),
                                ...) {
  
  if (!is.na(filename)) object$filename <- filename
  if (!object$validate) validationshape <- NA
  if (any(!comps %in% get_relevant_pcs(object = object, effect = "time"))) warning("Please note: Some components have low explanatory power and HAVE NOT BEEN rotated during rotation. Proceed with care.")
  
  if (validationshape == "cross" & !is.na(validationshape)) {
    dff <- subset(get_scores(object)$time, PC %in% comps)
    dff$PC <- paste0("PC", dff$PC)
    dff <- reshape2::melt(data = dff, id.vars = c("PC", "time", "group"))
    dff <- reshape2::dcast(data = dff, time + group ~ PC + variable, value.var = "value")
    dff$time <- factor(dff$time, levels = object$timelist)
    dff$group <- factor(dff$group, levels = object$grouplist)
    g <- ggplot2::ggplot(dff, ggplot2::aes_string(
      x = paste0("PC", comps[1], "_score"),
      y = paste0("PC", comps[2], "_score"),
      shape = "time",
      color = "group",
      group = "group",
      linetype = "group"
    )) +
      ggplot2::geom_line() +
      ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high"))) +
      ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")))
    g <- .getPlotHandle(
      g = g,
      object = object,
      myTheme = myTheme,
      comp = comps,
      effect = "time"
    )

    if (object$separateTimeAndGroup) {
      dff <- subset(get_scores(object)$group, PC %in% comps)
      dff$PC <- paste0("PC", dff$PC)
      dff <- reshape2::melt(data = dff, id.vars = c("PC", "time", "group"))
      dff <- reshape2::dcast(data = dff, time + group ~ PC + variable, value.var = "value")
      dff$time <- factor(dff$time, levels = object$timelist)
      dff$group <- factor(dff$group, levels = object$grouplist)
      gsg <- ggplot2::ggplot(dff, ggplot2::aes_string(
        x = paste0("PC", comps[1], "_score"),
        y = paste0("PC", comps[2], "_score"),
        shape = "time",
        color = "group",
        group = "group",
        linetype = "group"
      )) +
        ggplot2::geom_line() +
        ggplot2::geom_pointrange(ggplot2::aes_string(xmin = paste0("PC", comps[1], "_low"), xmax = paste0("PC", comps[1], "_high"))) +
        ggplot2::geom_pointrange(ggplot2::aes_string(ymin = paste0("PC", comps[2], "_low"), ymax = paste0("PC", comps[2], "_high")))
      gsg <- .getPlotHandle(
        g = gsg,
        object = object,
        myTheme = myTheme,
        comp = comps,
        effect = "group",
        legend = "bottom"
      )

      g <- ggpubr::ggarrange(g, gsg, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(gsg), align = "hv")
    }
  } else if (validationshape == "ellipse" & !is.na(validationshape)) {
    # Score plots
    ## Time
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC %in% comps),
        model = x
      )
    }))
    df_temp <- subset(get_scores(object)$time, PC %in% comps)
    df_temp$model <- 0
    df_temp$low <- NULL
    df_temp$high <- NULL
    dff <- rbind(dff, df_temp)

    dff$PC <- paste0("PC", dff$PC)
    dfff <- reshape2::dcast(data = dff, time + group + model ~ PC, value.var = "score")
    dfff$time <- factor(dfff$time, levels = object$timelist)
    dfff$group <- factor(dfff$group, levels = object$grouplist)
    dfff$alpha <- ifelse(dfff$model == 0, 1, alphavalidate)
    dfff <- dfff[order(dff$time), ]
    g <- ggplot2::ggplot(dfff, ggplot2::aes_string(paste0("PC", comps[1]), paste0("PC", comps[2]), shape = "time", color = "group")) +
      ggplot2::geom_line(data = subset(dfff, model == 0), ggplot2::aes(group = paste(model, group), linetype = group)) +
      ggplot2::geom_point(alpha = dfff$alpha) +
      ggplot2::stat_ellipse(level = validationlevel)
    if (plottext) {
      g <- g + ggplot2::geom_text(
        data = subset(dfff, model == 0),
        color = "black",
        ggplot2::aes(label = time), hjust = texthjust
      )
    }
    g <- .getPlotHandle(g = g, object = object, myTheme = myTheme, comp = comps, effect = "time")
    if (object$separateTimeAndGroup) {
      ## Group
      dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
        data.frame(
          subset(get_scores(object$validation$temp_objects[[x]])$group, PC %in% comps),
          model = x
        )
      }))
      df_temp <- subset(get_scores(object)$group, PC %in% comps)
      df_temp$model <- 0
      df_temp$low <- NULL
      df_temp$high <- NULL
      dff <- rbind(dff, df_temp)

      dff$PC <- paste0("PC", dff$PC)
      dfff <- reshape2::dcast(data = dff, time + group + model ~ PC, value.var = "score")
      dfff$time <- factor(dfff$time, levels = object$timelist)
      dfff$group <- factor(dfff$group, levels = object$grouplist)
      dfff$alpha <- ifelse(dfff$model == 0, 1, alphavalidate)
      dfff <- dfff[order(dff$time), ]
      gsg <- ggplot2::ggplot(dfff, ggplot2::aes_string(paste0("PC", comps[1]), paste0("PC", comps[2]), shape = "time", color = "group")) +
        ggplot2::geom_line(data = subset(dfff, model == 0), ggplot2::aes(group = paste(model, group), linetype = group)) +
        ggplot2::geom_point(alpha = dfff$alpha) +
        ggplot2::stat_ellipse(level = validationlevel)
      if (plottext) {
        gsg <- gsg + ggplot2::geom_text(data = subset(dfff, model == 0), color = "black", ggplot2::aes(label = time), hjust = texthjust)
      }

      gsg <- .getPlotHandle(g = gsg, object = object, myTheme = myTheme, comp = comps, effect = "group", legend = "bottom")

      g <- ggpubr::ggarrange(g, gsg, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(gsg), align = "hv")
    }
  } else {
    dff <- subset(get_scores(object)$time, PC %in% comps)
    dff$PC <- paste0("PC", dff$PC)
    dff <- reshape2::dcast(data = dff, time + group ~ PC, value.var = "score")
    dff$time <- factor(dff$time, levels = object$timelist)
    dff$group <- factor(dff$group, levels = object$grouplist)
    g <- ggplot2::ggplot(dff, ggplot2::aes_string(x = paste0("PC", comps[1]), y = paste0("PC", comps[2]), shape = "time", color = "group", group = "group", linetype = "group")) +
      ggplot2::geom_line() +
      ggplot2::geom_point()
    g <- .getPlotHandle(g = g, object = object, myTheme = myTheme, comp = comps, effect = "time")
    if (object$separateTimeAndGroup) {
      dff <- subset(get_scores(object)$group, PC %in% comps)
      dff$PC <- paste0("PC", dff$PC)
      dff <- reshape2::dcast(data = dff, time + group ~ PC, value.var = "score")
      dff$time <- factor(dff$time, levels = object$timelist)
      dff$group <- factor(dff$group, levels = object$grouplist)
      gsg <- ggplot2::ggplot(dff, ggplot2::aes_string(x = paste0("PC", comps[1]), y = paste0("PC", comps[2]), shape = "time", color = "group", "group" = "group", linetype = "group")) +
        ggplot2::geom_line() +
        ggplot2::geom_point()
      gsg <- .getPlotHandle(
        g = gsg,
        object = object,
        myTheme = myTheme,
        comp = comps,
        effect = "group",
        legend = "bottom"
      )

      g <- ggpubr::ggarrange(g, gsg, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(gsg))
    }
  }
  if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
  
  return(g)
}

#' Get color palette
#'
#' This function returns a list with colors for plotting
#'
#' @param object An ALASCA object
#' @return A list with colors
#'
#' @export
getPlotPalette <- function(object) {
  if (any(is.na(object$plot.palette)) || any(is.null(object$plot.palette))) {
    plotColors <- scales::viridis_pal(end = object$plot.palette.end)(length(object$grouplist))
    names(plotColors) <- object$grouplist
  } else {
    plotColors <- object$plot.palette
  }
  return(plotColors)
}

#' Get linetypes
#'
#' This function returns a list with linetypes for plotting
#'
#' @param object An ALASCA object
#' @return A list with linetypes
#'
#' @export
getPlotLinetypes <- function(object) {
  plotLinestypes <- scales::linetype_pal()(length(object$grouplist))
  names(plotLinestypes) <- object$grouplist
  return(plotLinestypes)
}

#' Get plotting features
#'
#' This function returns a list with colors for plotting
#'
#' @param object An ALASCA object
#' @return A list with colors
#'
#' @export
.getPlotHandle <- function(g, object, myTheme, effect = "time", comps = c(1, 2), legend = NA) {
  g <- g + ggplot2::scale_color_manual(values = getPlotPalette(object)) +
    ggplot2::labs(
      x = .getExpLabel(object, component = comps[1], effect = effect),
      y = .getExpLabel(object, component = comps[2], effect = effect)
    ) +
    myTheme
  
  if (!is.na(legend)) g <- g + ggplot2::theme(legend.position = legend)
  
  return(g)
}

#' Get exploratory power for plot label
#'
#' This function returns ...
#'
#' @param object An ALASCA object
#' @param comp Which two components to plot (default: `c(1, 2`)
#' @return A ggplot2 objects.
.getExpLabel <- function(object, component = 1, effect = "time", type = "Score") {
  if (effect == "time") {
    paste0(type, " PC", component, " (", round(100 * object$ALASCA$loading$explained$time[component], 2), "%)")
  } else {
    paste0(type, " PC", component, " (", round(100 * object$ALASCA$loading$explained$group[component], 2), "%)")
  }
}

#' Plot projection of participants
#'
#' This function returns a plot of...
#'
#' @param object An ALASCA object
#' @param comp Which two components to plot (default: `c(1, 2`)
#' @param return_data Set to `TRUE` to return data instead of plot
#' @param filetype Which filetype you want to save the figure to
#' @param figsize A vector containing `c(widht,height,dpi)` (default: `c(120, 80, 300)`)
#' @param myTheme A ggplot2 theme to use
#' @return A ggplot2 objects.
#'
#' @export
plotProjection <- function(object,
                           comp = c(1, 2),
                           return_data = FALSE,
                           filename = NA,
                           filetype = NA,
                           figsize = NA,
                           figunit = NA,
                           myTheme = NA) {
  if (any(is.na(myTheme))) myTheme <- object$plot.myTheme
  if (!is.na(filename)) object$filename <- filename

  df <- object$df
  df$ID <- df[, ID]
  loadings_Time <- subset(get_loadings(object)$time, PC %in% comp)
  loadings_Time <- reshape2::dcast(data = loadings_Time, covars ~ paste0("PC", PC), value.var = "loading")
  df_time <- merge(df, loadings_Time, by.x = "variable", by.y = "covars")
  if (object$separateTimeAndGroup) {
    df_time <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x) {
      data.frame(
        part = subset(df_time, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_time, paste0(ID, time) == x)$PC1 * subset(df_time, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_time, paste0(ID, time) == x)$PC2 * subset(df_time, paste0(ID, time) == x)$value),
        time = subset(df_time, paste0(ID, time) == x)$time
      )
    }))
    df_time <- df_time[!duplicated(df_time), ]
    colnames(df_time) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time")
    if (!return_data) {
      g_t <- ggplot2::ggplot(df_time, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part")) +
        ggplot2::geom_point() +
        ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.20, "cm"))) +
        myTheme
    }
    loadings_group <- subset(get_loadings(object)$group, PC %in% comp)
    loadings_group <- reshape2::dcast(data = loadings_group, covars ~ paste0("PC", PC), value.var = "loading")
    df_group <- merge(df, loadings_group, by.x = "variable", by.y = "covars")
    df_group <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x) {
      data.frame(
        part = subset(df_group, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_group, paste0(ID, time) == x)$PC1 * subset(df_group, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_group, paste0(ID, time) == x)$PC2 * subset(df_group, paste0(ID, time) == x)$value),
        time = subset(df_group, paste0(ID, time) == x)$time,
        group = subset(df_group, paste0(ID, time) == x)$group
      )
    }))
    df_group <- df_group[!duplicated(df_group), ]
    colnames(df_group) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time", "group")
    if (!return_data) {
      g_g <- ggplot2::ggplot(df_group, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part", color = "group")) +
        ggplot2::geom_point() +
        ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.20, "cm"))) +
        myTheme
      g <- ggpubr::ggarrange(g_t, g_g, common.legend = TRUE, legend = "bottom")
    }
    g <- list(df_time, df_group)
    names(g) <- c("time", "group")
    if (object$save) {
      for (i in seq_along(g)) {
        saveALASCAPlot(object = object, g = g[[i]], filetype = filetype, figsize = figsize, figunit = figunit, suffix = names(g)[i])
      }
    }
    return(g)
  } else {
    df_time <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x) {
      data.frame(
        part = subset(df_time, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_time, paste0(ID, time) == x)$PC1 * subset(df_time, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_time, paste0(ID, time) == x)$PC2 * subset(df_time, paste0(ID, time) == x)$value),
        time = subset(df_time, paste0(ID, time) == x)$time,
        group = subset(df_time, paste0(ID, time) == x)$group
      )
    }))
    df_time <- df_time[!duplicated(df_time), ]
    colnames(df_time) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time", "group")
    g <- ggplot2::ggplot(df_time, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part", color = "group")) +
      ggplot2::geom_point() +
      ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.20, "cm"))) +
      myTheme
    if (return_data) {
      return(df_time)
    } else {
      
      if (object$save) saveALASCAPlot(object = object, g = g, filetype = filetype, figsize = figsize, figunit = figunit)
      
      return(g)
    }
  }
}
