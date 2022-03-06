#' Get residuals
#'
#' This function returns the residual of the underlying linear mixed models
#'
#' @param object An ALASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @return A list of residuals per variable
#'
#' @examples
#' load("PE.Rdata")
#' residuals(model, variable = c("IL-6", "PlGF"))
#' @export
residuals.ALASCA <- function(object, variable = NA) {
  if (any(is.na(variable))) {
    return(lapply(object$regression_model, residuals))
  } else {
    varList <- names(object$regression_model)
    resList <- lapply(seq_along(object$regression_model), function(x) {
      if (varList[x] %in% variable) {
        residuals(object$regression_model[[x]])
      }
    })
    names(resList) <- names(object$regression_model)
    resList[vapply(resList, is.null)] <- NULL
    return(resList)
  }
}

#' Plot residuals
#'
#' This function plots the residuals of the underlying linear mixed models
#'
#' @param object An ALASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @param plottitle If `TRUE` (default), include variable name as title
#' @return A list of ggplot2 objects per variable
#'
#' @export
plotResiduals <- function(object, variable = NA, plottitle = TRUE, myTheme = ggplot2::theme_classic()) {
  resList <- residuals(object, variable = variable)
  lapply(seq_along(resList), function(x) {
    g <- ggplot2::ggplot(
      data.frame(Residuals = resList[[x]]),
      ggplot2::aes(sample = Residuals)
    ) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line()
    if (plottitle) g <- g + ggplot2::labs(title = names(resList)[x], x = "Theoretical", y = "Sample")
    g <- g + myTheme
    if (object$save) saveALASCAPlot(object, g, prefix = "plot/", suffix = paste0("_qq_plot_", names(resList)[x]))
    g
  })
}

#' Plot histogram of validation
#'
#' This function plots the validation runs as histograms
#'
#' @param object An ALASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @param plottitle If `TRUE` (default), include variable name as title
#' @return A list of ggplot2 objects per variable
#'
#' @export
plotHistogram <- function(object,
                          component = 1,
                          bins = object$nValRuns / 10,
                          n.limit = 0,
                          variable = NA,
                          filename = NA,
                          effect = "time",
                          orderbyname = FALSE) {
  if (!is.na(filename)) object$filename <- filename

  if (effect == "time" | effect == "group") {
    g_s <- plothistogram_score(object = object, component = component, bins = bins, effect = effect)
    g_l <- plothistogram_loading(object = object, component = component, bins = bins, variable = variable, effect = effect, orderbyname = orderbyname, n.limit = n.limit)
    g <- ggpubr::ggarrange(
      g_s, g_l,
      nrow = 2, heights = c(1, 3), labels = "AUTO"
    )
    if (object$save) saveALASCAPlot(object, g, prefix = "plot/", suffix = paste0("_histo"), figsize = c(200, 240, 300))
  } else {
    g <- list(
      time = plothistogram(object = object, component = component, bins = bins, variable = variable, effect = "time", orderbyname = orderbyname),
      group = plothistogram(object = object, component = component, bins = bins, variable = variable, effect = "group", orderbyname = orderbyname)
    )
  }

  return(g)
}

#' Plot histogram of validation
#'
#' This function plots the validation runs as histograms
#'
#' @param object An ALASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @param plottitle If `TRUE` (default), include variable name as title
#' @return A list of ggplot2 objects per variable
#'
#' @export
plothistogram_score <- function(object, component = 1, bins = object$nValRuns / 10, effect = "time") {
  if (!object$validate) stop("Model not validated")
  if (effect == "time" | effect == "both") {
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    }))
    df_temp <- subset(get_scores(object)$time, PC == component)
    df_temp$model <- 0
    df_temp$low <- NULL
    df_temp$high <- NULL
    if (object$separateTimeAndGroup) {
      dff$group <- levels(object$df$group)[1]
      df_temp$group <- levels(object$df$group)[1]
    }

    g <- ggplot2::ggplot(dff, ggplot2::aes(score, fill = group)) +
      ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = bins) +
      ggplot2::geom_vline(data = df_temp, ggplot2::aes(xintercept = score, color = group)) +
      ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
      ggplot2::scale_color_manual(values = getPlotPalette(object)) +
      ggplot2::facet_wrap(~time) +
      object$plot.myTheme +
      ggplot2::theme(legend.position = "bottom")
  }
  if (effect == "group" | effect == "both") {
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        subset(get_scores(object$validation$temp_objects[[x]])$group, PC == component),
        model = x
      )
    }))
    df_temp <- subset(get_scores(object)$group, PC == component)
    df_temp$model <- 0
    df_temp$low <- NULL
    df_temp$high <- NULL

    if (effect == "group") {
      g <- ggplot2::ggplot(dff, ggplot2::aes(score, fill = group)) +
        ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = bins) +
        ggplot2::geom_vline(data = df_temp, ggplot2::aes(xintercept = score, color = group)) +
        ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
        ggplot2::scale_color_manual(values = getPlotPalette(object)) +
        ggplot2::facet_wrap(~time) +
        ggplot2::labs(x = "Score", fill = object$plot.grouplabel, color = object$plot.grouplabel) +
        object$plot.myTheme +
        ggplot2::theme(legend.position = "bottom")
    } else {
      gg <- ggplot2::ggplot(dff, ggplot2::aes(score, fill = group)) +
        ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = bins) +
        ggplot2::geom_vline(data = df_temp, ggplot2::aes(xintercept = score, color = group)) +
        ggplot2::scale_fill_manual(values = getPlotPalette(object)) +
        ggplot2::scale_color_manual(values = getPlotPalette(object)) +
        ggplot2::facet_wrap(~time) +
        ggplot2::labs(x = "Score", fill = object$plot.grouplabel, color = object$plot.grouplabel) +
        object$plot.myTheme +
        ggplot2::theme(legend.position = "bottom")
      g <- ggpubr::ggarrange(g, gg, nrow = 2, labels = "AUTO", common.legend = TRUE)
    }
  }
  if (object$save) saveALASCAPlot(object, g, prefix = "plot/", suffix = paste0("_histo_score"))
  return(g)
}

#' Plot histogram of validation
#'
#' This function plots the validation runs as histograms
#'
#' @param object An ALASCA object
#' @param variable The name of the variable(s) you want. `NA` returns all (default).
#' @param plottitle If `TRUE` (default), include variable name as title
#' @return A list of ggplot2 objects per variable
#'
#' @export
plothistogram_loading <- function(object, component = 1,
                                  bins = object$nValRuns / 10,
                                  variable = NA,
                                  effect = "time",
                                  orderbyname = FALSE,
                                  n.limit = 0) {
  if (!object$validate) stop("Model not validated")
  if (any(is.na(variable))) {
    variable <- object$variablelist
  }
  if (effect == "time") {
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        get_loadings(object$validation$temp_objects[[x]], component = component, n.limit = n.limit)$time,
        model = x
      )
    }))
    df_temp <- subset(get_loadings(object)$time, PC == component & covars %in% variable)
  } else {
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x) {
      data.frame(
        get_loadings(object$validation$temp_objects[[x]], component = component, n.limit = n.limit)$group,
        model = x
      )
    }))
    df_temp <- subset(get_loadings(object)$group, PC == component & covars %in% variable)
  }

  df_temp$model <- 0
  df_temp$low <- NULL
  df_temp$high <- NULL
  if (!orderbyname) dff$covars <- factor(dff$covars, levels = df_temp$covars[order(df_temp$loading, decreasing = TRUE)])
  if (!orderbyname) df_temp$covars <- factor(df_temp$covars, levels = df_temp$covars[order(df_temp$loading, decreasing = TRUE)])
  g <- ggplot2::ggplot(dff, ggplot2::aes(loading)) +
    ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = bins) +
    ggplot2::geom_vline(data = df_temp, ggplot2::aes(xintercept = loading)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~covars) +
    ggplot2::labs(x = "Loading") +
    object$plot.myTheme +
    ggplot2::theme(legend.position = "bottom")
  if (object$save) saveALASCAPlot(object, g, prefix = "plot/", suffix = paste0("_histo_loading"))
  return(g)
}

#' Count participants and samples
#'
#' Count participants and samples
#'
#' @param object An ALASCA object
#' @return A list of data frames
#'
#' @export
countParts <- function(object) {
  wideDF <- dcast(data = object$df, ID + group ~ time + variable, fun.aggregate = length)
  samples.by.participant <- data.frame(
    ID = wideDF$ID,
    count = rowSums(wideDF[, 3:ncol(wideDF)])
  )
  complete.cases.by.group <- as.data.frame(table(wideDF[rowSums(wideDF[, 3:ncol(wideDF)]) == length(object$variablelist) * length(object$timelist), "group"]))
  colnames(complete.cases.by.group) <- c("group", "count")

  list(
    participants.by.group = object$df[, .(count = uniqueN(ID)), by = .(group)],
    complete.cases.by.group = complete.cases.by.group,
    participants.by.time = object$df[, .(count = uniqueN(ID)), by = .(time)],
    participants.by.time.and.group = dcast(data = object$df[, .(count = uniqueN(ID)), by = .(time, group)], group ~ time, value.var = "count"),
    participants.by.variable.and.time.and.group = dcast(data = object$df[, .(count = uniqueN(ID)), by = .(variable, time, group)], time + group ~ variable, value.var = "count"),
    samples.by.variable.and.time.and.group = dcast(data = object$df[, .(count = .N), by = .(variable, time, group)], time + group ~ variable, value.var = "count"),
    samples.by.participant = samples.by.participant,
    samples.by.participant.and.and.time = dcast(data = object$df, ID + time ~ ., value.var = "variable", fun.aggregate = length)
  )
}
