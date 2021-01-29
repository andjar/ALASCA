#' Get an RMASCA object
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An RMASCA object
#' @param component String stating which component to return (PC1 is default)
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param decreasingLoadings Sort the loadings in dedreasing (default) or increasing order
#' @param only String stating which plot to return; `both` (default), `score` or `loading`
#' @param enlist Logical. If TRUE, the plots are returned as a list and not as a composed figure (default)
#' @return An ggplot2 object (or a list og ggplot objects)
#' @export
plot.RMASCA <- function(object, component = "PC1", effect = "both", decreasingLoadings = TRUE, only = "both", enlist = FALSE){
  if(!(effect %in% c("both","time","group"))){
    stop("`effect` has to be `both`, `time` or `group`")
  }
  if(!object$separateTimeAndGroup){
    effect = "time"
  }
  if(only == "score"){
    if(effect == "both"){
      g_score_time <- getScorePlot(object, component = component, effect = "time")
      g_score_group <- getScorePlot(object, component = component, effect = "group")
      g <- list(g_score_time, g_score_group)
    }else{
      g <- getScorePlot(object, component = component, effect = effect)
    }
  }else if(only == "loading"){
    if(effect == "both"){
      g_loading_time <- getLoadingPlot(object, component = component, effect = "time", decreasingLoadings = decreasingLoadings)
      g_loading_group <- getLoadingPlot(object, component = component, effect = "group", decreasingLoadings = decreasingLoadings)
      g <- list(g_loading_time, g_loading_group)
    }else{
      g <- getLoadingPlot(object, component = component, effect = effect, decreasingLoadings = decreasingLoadings)
    }
  }else{
    if(effect == "both"){
      g_loading_time <- getLoadingPlot(object, component = component, effect = "time", decreasingLoadings = decreasingLoadings)
      g_loading_group <- getLoadingPlot(object, component = component, effect = "group", decreasingLoadings = decreasingLoadings)
      g_score_time <- getScorePlot(object, component = component, effect = "time")
      g_score_group <- getScorePlot(object, component = component, effect = "group")
      if(enlist){
        g <- list(g_score_time, g_loading_time, g_score_group, g_loading_group)
      }else{
        g <- ggpubr::ggarrange(g_score_time, g_loading_time, g_score_group, g_loading_group, nrow = 2, ncol = 2, align = "hv", common.legend = TRUE, legend = "bottom")
      }
    }else{
      g_loading <- getLoadingPlot(object, component = component, effect = effect, decreasingLoadings = decreasingLoadings)
      g_score <- getScorePlot(object, component = component, effect = effect)
      if(enlist){
        g <- list(g_score, g_loading)
      }else{
        g <- ggpubr::ggarrange(g_score, g_loading, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE, legend = "bottom")
      }
    }
  }
  return(g)
}

#' Get screeplot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An RMASCA object
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @return An ggplot2 object (or a list og ggplot objects)
#' @export
screeplot.RMASCA <- function(object, effect = "both"){
  explained <- as.data.frame(getScores(object)$explained)
  explained$component <- 1:nrow(explained)
  g <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = time, group = NA)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Principal Component", y = "Relative Expl. of Time Var.")
  if(length(object$RMASCA$loading) == 3){
    gg <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = group, group = NA)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Principal Component", y = "Relative Expl. of Group Var.")
    if(effect == "both"){
      g <- ggpubr::ggarrange(g, gg)
    }
  }
  if(effect == "group"){
    return(gg)
  }else{
    return(g)
  }
}

#' Get loadings
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An RMASCA object
#' @return A list with loadings for time (and group), and the exploratory power for each component
#' @export
getLoadings <- function(object){
  object$RMASCA$loading$time <- object$RMASCA$loading$time[!duplicated(object$RMASCA$loading$time),]
  if(length(object$RMASCA$loading) == 3){
    object$RMASCA$loading$group <- object$RMASCA$loading$group[!duplicated(object$RMASCA$loading$group),]
  }
  return(object$RMASCA$loading)
}

#' Get scores
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An RMASCA object
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
getScores <- function(object){
  object$RMASCA$score$time <- object$RMASCA$score$time[!duplicated(object$RMASCA$score$time),]
  if(length(object$RMASCA$score) == 3){
    object$RMASCA$score$group <- object$RMASCA$score$group[!duplicated(object$RMASCA$score$group),]
  }
  return(object$RMASCA$score)
}

#' Get loading plot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An RMASCA object
#' @param component Which component to plot?
#' @param effect Plot time or group
#' @param decreasingLoadings Logical. Should loading sbe sorted in decreasing order?
#' @return A ggplot object
getLoadingPlot <- function(object, component = "PC1", effect = "time", decreasingLoadings = TRUE){
  pointSize <- 0.4
  if(effect == "time"){
    loadings <- getLoadings(object)$time
    PC <- which(colnames(loadings) == component)
    loadings <- loadings[,colnames(loadings) %in% c(component,"covars")]
    colnames(loadings) <- c("loading", "covars")
    if(object$validate){
      loadings_unc <- object$validation$time$loading[object$validation$time$loading$PC == PC,]
      loadings <- merge(loadings, loadings_unc, by = "covars")
    }
  }else{
    loadings <- getLoadings(object)$group
    PC <- which(colnames(loadings) == component)
    loadings <- loadings[,colnames(loadings) %in% c(component,"covars")]
    colnames(loadings) <- c("loading", "covars")
    if(object$validate){
      loadings_unc <- object$validation$group$loading[object$validation$group$loading$PC == PC,]
      loadings <- merge(loadings, loadings_unc, by = "covars")
    }
  }
  loadings$covars = factor(loadings$covars, levels = unique(loadings$covars[order(loadings$loading, decreasing = decreasingLoadings)]))
  if(object$validate){
    g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high)) + ggplot2::geom_pointrange(size = pointSize)
  }else{
    g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading)) + ggplot2::geom_point()
  }

  g <- g +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
    ggplot2::labs(x = "Variable",
                  y = paste0(component, " (", round(100*ifelse(effect == "time", object$RMASCA$loading$explained$time[PC],object$RMASCA$loading$explained$group[PC]),2),"%)"))
  return(g)
}

#' Get score plot
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object An RMASCA object
#' @param component Which component to plot?
#' @param effect Plot time or group
#' @return A ggplot object
getScorePlot <- function(object, component = "PC1", effect = "time"){
  pointSize <- 0.4
  PC <- which(colnames(object$RMASCA$score$time) == component)
  if(effect == "time"){
    if(object$separateTimeAndGroup){
      score <- data.frame(
                          score = object$RMASCA$score$time[,PC],
                          time = object$parts$time
                          )
      score <- score[!duplicated(score),]
      if(object$validate){
        score_unc <- object$validation$time$score[object$validation$time$score$PC == PC,]
        score <- merge(score, score_unc, by = "time")
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = NA, ymin = low, ymax = high)) +
          ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35), size = pointSize) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35))
      }else{
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = NA)) + ggplot2::geom_point() + ggplot2::geom_line()
      }
    }else{
      score <- data.frame(
        score = object$RMASCA$score$time[,PC],
        time = object$parts$time,
        group = object$parts$group
      )
      score <- score[!duplicated(score),]
      if(object$validate){
        score_unc <- object$validation$time$score[object$validation$time$score$PC == PC,]
        score <- merge(score, score_unc, by = c("time", "group"))
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, ymin = low, ymax = high)) +
          ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35), size = pointSize) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35))
      }else{
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group)) +
          ggplot2::geom_point() + ggplot2::geom_line()
      }
    }
    g <- g +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = "Time", y = paste0(component, " (",round(100*object$RMASCA$score$explained$time[PC],2),"%)"))
  }else{
    score <- data.frame(
      score = object$RMASCA$score$group[,PC],
      time = object$parts$time,
      group = object$parts$group
    )
    if(object$validate){
      score_unc <- object$validation$group$score[object$validation$group$score$PC == PC,]
      score <- merge(score, score_unc, by = c("time", "group"))
      score <- score[!duplicated(score),]
      g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, ymin = low, ymax = high)) +
        ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35), size = pointSize) +
        ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35))
    }else{
      g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group)) +
        ggplot2::geom_point() + ggplot2::geom_line()
    }
    g <- g +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = "Time", y = paste0(component, " (",round(100*object$RMASCA$score$explained$group[PC],2),"%)"))
  }
  return(g)
}
