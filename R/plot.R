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
        g <- ggpubr::ggarrange(g_loading, g_score, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE, legend = "bottom")
      }
    }
  }
  return(g)
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

getLoadingPlot <- function(object, component = "PC1", effect = "time", decreasingLoadings = TRUE){
  if(effect == "time"){
    loadings <- getLoadings(object)$time
  }else{
    loadings <- getLoadings(object)$group
  }
  PC <- which(colnames(loadings) == component)
  loadings <- loadings[,colnames(loadings) %in% c(component,"covars")]
  colnames(loadings) <- c("loading", "covars")
  loadings$covars = factor(loadings$covars, levels = unique(loadings$covars[order(loadings$loading, decreasing = decreasingLoadings)]))
  g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Variable", y = paste0(component, " (",round(100*object$RMASCA$loading$explained$time[PC],2),"%)"))
  return(g)
}

getScorePlot <- function(object, component = "PC1", effect = "time"){
  PC <- which(colnames(object$RMASCA$score$time) == component)
  if(effect == "time"){
    if(object$separateTimeAndGroup){
      score <- data.frame(
                          score = object$RMASCA$score$time[,PC],
                          time = object$parts$time
                          )
      score <- score[!duplicated(score),]
      g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = NA)) +
              ggplot2::geom_point() +
              ggplot2::geom_line() +
              ggplot2::labs(x = "Time", y = paste0(component, " (",round(100*object$RMASCA$score$explained$time[PC],2),"%)"))
    }else{
      score <- data.frame(
        score = object$RMASCA$score$time[,PC],
        time = object$parts$time,
        group = object$parts$group
      )
      score <- score[!duplicated(score),]
      g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x = "Time", y = paste0(component, " (",round(100*object$RMASCA$score$explained$time[PC],2),"%)"))
    }
  }else{
    score <- data.frame(
      score = object$RMASCA$score$group[,PC],
      time = object$parts$time,
      group = object$parts$group
    )
    score <- score[!duplicated(score),]
    g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Time", y = paste0(component, " (",round(100*object$RMASCA$score$explained$time[PC],2),"%)"))
  }
  return(g)
}
