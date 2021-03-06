#' Get an ALASCA object
#'
#' This function plots your ALASCA model
#'
#' @param object An [ALASCA()] object
#' @param component Integer stating which component to return (1 is default)
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param decreasingLoadings Sort the loadings in decreasing (default) or increasing order
#' @param only String stating which plot to return; `both` (default), `score` or `loading`
#' @param enlist Logical. If TRUE, the plots are returned as a list and not as a composed figure (default)
#' @param tooDense Integer, If > 0, only name this number of covariables
#' @param xlabel Defaults to "Time" if not specified during model setup
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
#' 
#' @export
plot.ALASCA <- function(object, 
                        component = 1, 
                        effect = "both", 
                        decreasingLoadings = TRUE, 
                        only = "both", 
                        enlist = FALSE, 
                        tooDense = NA, 
                        highlight = NA, 
                        xlabel = NA,
                        myTheme = ggplot2::theme_bw()){
  if(!(effect %in% c("both","time","group"))){
    stop("`effect` has to be `both`, `time` or `group`")
  }
  if(!object$separateTimeAndGroup){
    effect = "time"
  }
  if(!is.na(xlabel)){
    object$plot.xlabel <- xlabel
  }
  if(only == "score"){
    if(effect == "both"){
      g_score_time <- getScorePlot(object, component = component, effect = "time", myTheme = myTheme)
      g_score_group <- getScorePlot(object, component = component, effect = "group", myTheme = myTheme)
      g <- list(g_score_time, g_score_group)
    }else{
      g <- getScorePlot(object, component = component, effect = effect)
    }
  }else if(only == "loading"){
    if(effect == "both"){
      g_loading_time <- getLoadingPlot(object, component = component, effect = "time", decreasingLoadings = decreasingLoadings, tooDense = tooDense, highlight = highlight, myTheme = myTheme)
      g_loading_group <- getLoadingPlot(object, component = component, effect = "group", decreasingLoadings = decreasingLoadings, tooDense = tooDense, highlight = highlight, myTheme = myTheme)
      g <- list(g_loading_time, g_loading_group)
    }else{
      g <- getLoadingPlot(object, component = component, effect = effect, decreasingLoadings = decreasingLoadings, tooDense = tooDense, myTheme = myTheme)
    }
  }else{
    if(effect == "both"){
      g_loading_time <- getLoadingPlot(object, component = component, effect = "time", decreasingLoadings = decreasingLoadings, tooDense = tooDense, highlight = highlight, myTheme = myTheme)
      g_loading_group <- getLoadingPlot(object, component = component, effect = "group", decreasingLoadings = decreasingLoadings, tooDense = tooDense, highlight = highlight, myTheme = myTheme)
      g_score_time <- getScorePlot(object, component = component, effect = "time", myTheme = myTheme)
      g_score_group <- getScorePlot(object, component = component, effect = "group", myTheme = myTheme)
      if(enlist){
        g <- list(g_score_time, g_loading_time, g_score_group, g_loading_group)
      }else{
        g <- ggpubr::ggarrange(g_score_time, g_loading_time, g_score_group, g_loading_group, nrow = 2, ncol = 2, widths = c(1,2,1,2), align = "hv", common.legend = TRUE, legend = "bottom")
      }
    }else{
      g_loading <- getLoadingPlot(object, component = component, effect = effect, decreasingLoadings = decreasingLoadings, tooDense = tooDense, highlight = highlight, myTheme = myTheme)
      g_score <- getScorePlot(object, component = component, effect = effect, myTheme = myTheme)
      if(enlist){
        g <- list(g_score, g_loading)
      }else{
        g <- ggpubr::ggarrange(g_score, g_loading, nrow = 1, ncol = 2, widths = c(1,2,1,2), align = "hv", common.legend = TRUE, legend = "bottom")
      }
    }
  }
  return(g)
}

#' Get screeplot
#'
#' This function returns a screeplot for an ALASCA model showing what proportion of the variance each component of the model explains
#'
#' @param object An ALASCA object
#' @param effect String stating which effect to return; `time`, `group`, `both` (default)
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return An ggplot2 object (or a list og ggplot objects)
#' 
#' @examples
#' load("PE.Rdata")
#' screeplot(model)
#' 
#' @export
screeplot.ALASCA <- function(object, effect = "both", myTheme = ggplot2::theme_bw()){
  explained <- as.data.frame(getScores(object)$explained)
  explained$component <- 1:nrow(explained)
  g <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = time, group = NA)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    myTheme + 
    ggplot2::labs(x = "Principal Component", y = paste0("Relative Expl. of ",object$plot.xlabel," Var."))
  if(length(object$ALASCA$loading) == 3){
    gg <- ggplot2::ggplot(explained, ggplot2::aes(x = component, y = group, group = NA)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      myTheme + 
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
#' This function  returns the loadings for an ALASCA model
#'
#' @param object An ALASCA object
#' @return A list with loadings for time (and group), and the exploratory power for each component
#' @export
getLoadings <- function(object){
  return(object$ALASCA$loading)
}

#' Get scores
#'
#' This function returns the scores for an ALASCA model
#'
#' @inheritParams getLoadings
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
getScores <- function(object){
  return(object$ALASCA$score)
}

#' Get covariables
#'
#' This function returns the other covariables in an ALASCA model
#'
#' @inheritParams getLoadings
#' @return A list with scores for time (and group), and the exploratory power for each component
#' @export
getCovars <- function(object){
  return(object$CovarCoefficients)
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
#' @param decreasingLoadings Logical. Should loading sbe sorted in decreasing order?
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot object
getLoadingPlot <- function(object, component = 1, effect = "time", decreasingLoadings = TRUE, tooDense = NA, highlight = NA, myTheme = ggplot2::theme_bw()){
  pointSize <- 0.4
  if(effect == "time"){
    loadings <- subset(getLoadings(object)$time, PC == component)
  }else{
    loadings <- subset(getLoadings(object)$group, PC == component)
  }
  loadings$covars = factor(loadings$covars, levels = unique(loadings$covars[order(loadings$loading, decreasing = decreasingLoadings)]))
  if(object$validate){
    g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, ymin = low, ymax = high)) + ggplot2::geom_pointrange(size = pointSize)
  }else{
    if(any(colnames(loadings) == "model")){
      g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading, shape = model)) + ggplot2::geom_point()
    }else{
      g <- ggplot2::ggplot(loadings, ggplot2::aes(x = covars, y = loading)) + ggplot2::geom_point()
    }
  }
  
  g <- g + myTheme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1), legend.position = c(0.8, 0.8)) +
    ggplot2::labs(x = "Variable",
                  y = paste0("PC",component, " (", round(100*ifelse(effect == "time", object$ALASCA$loading$explained$time[component],object$ALASCA$loading$explained$group[component]),2),"%)"))
  if(!any(is.na(highlight))){
    g <- g + ggplot2::geom_point(color = ifelse(loadings$covars %in% highlight, "red", "grey")) +
      ggrepel::geom_text_repel(data = subset(loadings, covars %in% highlight), ggplot2::aes(label=covars), max.iter	= 5000) +
      myTheme +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(),
                     legend.position = "none")
  }else if(!is.na(tooDense) & tooDense > 0){
    limUpper <- unique(loadings$loading[order(loadings$loading, decreasing = TRUE)[tooDense]])
    limLower <- unique(loadings$loading[order(loadings$loading, decreasing = FALSE)[tooDense]])
    g <- g + ggplot2::geom_point(color = ifelse(loadings$loading <= limLower | loadings$loading >= limUpper, "red", "grey")) +
      ggrepel::geom_text_repel(data = subset(loadings, loading <= limLower | loading >= limUpper), ggplot2::aes(label=covars), max.iter	= 5000) +
      myTheme +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(),
                     legend.position = "none")
  }
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
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot object
getScorePlot <- function(object, component = 1, effect = "time", myTheme = ggplot2::theme_bw()){
  pointSize <- 2
  if(effect == "time"){
    if(object$separateTimeAndGroup){
      score <- subset(getScores(object)$time, PC == component)
      if(object$validate){
        if(grepl("permutation", object$validationMethod)){
          pvals <- object$pvals
          score <- merge(score, pvals, by.x = "time", by.y = "effect", all.x=TRUE, all.y=FALSE)
          score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = NA, label=p.value.str)) +
            ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
            ggplot2::geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE)+
            ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
        }else{
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = NA, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35)) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35))
        }
      }else{
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = NA)) + ggplot2::geom_point() + ggplot2::geom_line()
      }
    }else{
      score <- subset(getScores(object)$time, PC == component)
      if(object$validate){
        if(grepl("permutation", object$validationMethod)){
          pvals <- object$pvals
          score$effect <- paste(score$time, score$group)
          score <- merge(score, pvals, by.x = "effect", by.y = "effect", all.x=TRUE, all.y=FALSE)
          score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, label = p.value.str)) +
            ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
            geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
        }else{
          g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, ymin = low, ymax = high)) +
            ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35)) +
            ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35))
        }
      }else{
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group)) +
          ggplot2::geom_point() + ggplot2::geom_line()
      }
    }
    g <- g + myTheme +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = object$plot.xlabel, y = paste0("PC",component, " (",round(100*object$ALASCA$score$explained$time[component],2),"%)"))
  }else{
    score <- subset(getScores(object)$group, PC == component)
    if(object$validate){
      if(grepl("permutation", object$validationMethod)){
        pvals <- object$pvals
        score$effect <- paste(score$time, score$group)
        score <- merge(score, pvals, by.x = "effect", by.y = "effect", all.x=TRUE, all.y=FALSE)
        score$p.value.str <- ifelse(score$p.value > .05, "", ifelse(score$p.value < .001, "***", ifelse(score$p.value < .01, "**", "*")))
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, label = p.value.str)) +
          ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
          geom_text(vjust = 0, hjust = 0.5, position = ggplot2::position_dodge(width = 0.5), show.legend = FALSE) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.5))
      }else{
        g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group, ymin = low, ymax = high)) +
          ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35)) +
          ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35))
      }
    }else{
      g <- ggplot2::ggplot(score, ggplot2::aes(x = time, y = score, group = group, color = group)) +
        ggplot2::geom_point() + ggplot2::geom_line()
    }
    g <- g + myTheme +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = object$plot.xlabel, y = paste0("PC",component, " (",round(100*object$ALASCA$score$explained$group[component],2),"%)"))
  }
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
#' @param addSmooth. Specify which geom_smooth model you want to apply, eg. `lm`, `glm`, `gam`, `loess` (default). Set to `NA` to remove.
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#' 
#' @examples
#' load("PE.Rdata")
#' plotParts(df, variable = "IL-6", participantColumn = "ID", valueColumn = "value")
#' do.call(
#'     ggpubr::ggarrange,
#'     c(plotParts(df, participantColumn = "ID", timeColumn = "GA", valueColumn = "value", addSmooth = NA, 
#'     variable = c("PlGF", "IL-6", "IL-1b", "IFN-g", "Eotaxin-2", "Eotaxin")), 
#'       common.legend = TRUE, legend = "bottom")
#'   )
#' plotParts(model, variable = "IL-6", participantColumn = "ID", timeColumn = "GA")[[1]] + ggplot2::labs(x = "GA")
#' 
#' @export
plotParts <- function(object, 
                      variable = NA, 
                      participantColumn = FALSE, 
                      valueColumn = FALSE, 
                      timeColumn = "time", 
                      addSmooth = "loess", 
                      myTheme = ggplot2::theme_bw()){
  if(is.data.frame(object)){
    df <- object
    if(any(participantColumn == FALSE) | any(valueColumn == FALSE)){
      stop("You need to specify participant and value columns")
    }else{
      participantColumn <- participantColumn
      valueColumn <- valueColumn
    }
  }else if(is(object, "ALASCA")){
    df <- object$dfRaw
    valueColumn <- as.character(object$formula)[2]
    if(any(participantColumn == FALSE)){
      if(any(object$participantColumn == FALSE)){
        stop("You need to specify the participant column")
      }else{
        participantColumn <- object$participantColumn
      }
    }
  }else{
    stop("Wrong input object: must be a ALASCA model or a data frame")
  }
  plotFunction <- function(df, timeColumn, valueColumn, participantColumn, xi, addSmooth, myTheme){
    g <- ggplot2::ggplot(subset(df, variable == xi), ggplot2::aes_string(x = timeColumn, y = valueColumn, color = "group", group = participantColumn)) + 
      ggplot2::geom_point(alpha = 0.7) + ggplot2::geom_line(alpha = 0.3)  + myTheme +
      ggplot2::theme(legend.position = "bottom") + ggplot2::labs(x = "Time", y = xi)
    if(!any(is.na(addSmooth))){
      g <- g + ggplot2::geom_smooth(method = addSmooth, ggplot2::aes(group = group), se = TRUE)
    }
    return(g)
  }
  if(any(is.na(variable))){
    variable <- unique(df$variable)
  }
  g <- lapply(variable, function(xi){
    plotFunction(df, timeColumn, valueColumn, participantColumn, xi, addSmooth, myTheme = myTheme)
  })
  names(g) <- variable
  return(g)
}

#' Plot model predictions
#'
#' This function returns the scores for an ALASCA model
#'
#' @param object An ALASCA object or a data frame. If a data frame, you need to specify the column names for participant and value. This also applies if you have not specified the participant column in the ALASCA model before.
#' @param variable List of variable names to print. If `NA`, return all (default).
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#' 
#' @examples
#' load("PE.Rdata")
#' plotPred(model, variable = "IL-6")[[1]] + ggplot2::theme_bw()
#' do.call(
#'   ggpubr::ggarrange,
#'   c(plotPred(model, variable = c("PlGF", "IL-6", "IL-1b", "IFN-g", "Eotaxin-2", "Eotaxin")), 
#'     common.legend = TRUE, legend = "bottom")
#' )
#' @export
plotPred <- function(object, variable = NA, myTheme = ggplot2::theme_bw()){
  if(any(is.na(variable))){
    variable <- unique(object$df$variable)
  }
  if(object$validateRegression){
    gg <- lapply(unique(variable), function(x){
      g <- ggplot2::ggplot(subset(object$mod.pred, variable == x), ggplot2::aes(x = time, y = pred, color = group, group = group, ymin = low, ymax = high)) +
        ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.35)) + 
        ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35)) + 
        myTheme + 
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = object$plot.xlabel, y = x)
      g
    })
  }else{
    gg <- lapply(unique(variable), function(x){
      g <- ggplot2::ggplot(subset(object$mod.pred, variable == x), ggplot2::aes(x = time, y = pred, color = group, group = group)) +
        ggplot2::geom_point() + 
        ggplot2::geom_line() + 
        myTheme + 
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = object$plot.xlabel, y = x)
      g
    })
  }
  names(gg) <- variable
  return(gg)
}

#' Plot validations models
#'
#' This function returns a plot of the validation models
#'
#' @param object A validated ALASCA object
#' @param component Which component to plot (default: 1)
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A list with ggplot2 objects.
#' 
#' @export
plotVal <- function(object, component = 1, myTheme = ggplot2::theme_bw()){
  if(!object$validate){
    stop("You must validate the model first.")
  }
  if(object$separateTimeAndGroup){
    # Score plots
    ## Time
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x){
      data.frame(
        subset(getScores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    })
    )
    gst <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = model)) + 
      ggplot2::geom_point(alpha = 0.2) + ggplot2::geom_line(alpha = 0.2) +
      ggplot2::labs(x = object$plot.xlabel) + myTheme
    
    ## Group
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x){
      data.frame(
        subset(getScores(object$validation$temp_objects[[x]])$group, PC == component),
        model = x
      )
    })
    )
    dff$plotGroup <- paste0(dff$model,"-",dff$group)
    gsg <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = plotGroup, color = group)) + 
      ggplot2::geom_point(alpha = 0.2) + ggplot2::geom_line(alpha = 0.2) +
      ggplot2::geom_point(data = subset(getScores(object)$group, PC == component), group = NA, alpha = 1, color = "black") +
      ggplot2::geom_line(data = subset(getScores(object)$group, PC == component), group = subset(getScores(object)$group, PC == component)$group, alpha = 1, color = "black") +
      ggplot2::labs(x = object$plot.xlabel) + myTheme
    
    # Loading plot
    ## Time
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x){
      data.frame(
        subset(getLoadings(object$validation$temp_objects[[x]])$time, PC == component)
      )
    })
    )
    glt <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) + 
      ggplot2::geom_point(alpha = 0.2, color = "black") + 
      ggplot2::geom_point(data = subset(getLoadings(object)$time, PC == component), alpha = 1, color = "red") +
      ggplot2::labs(x = "Variable") + myTheme
    
    ## Group
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x){
      data.frame(
        subset(getLoadings(object$validation$temp_objects[[x]])$group, PC == component)
      )
    })
    )
    glg <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) + 
      ggplot2::geom_point(alpha = 0.2, color = "black") + 
      ggplot2::geom_point(data = subset(getLoadings(object)$group, PC == component), alpha = 1, color = "red") +
      ggplot2::labs(x = "Variable") + myTheme
    
    g <- ggpubr::ggarrange(gst, glt, gsg, glg, nrow = 2, ncol = 2, widths = c(1,2,1,2))
    
  }else{
    # Score plot
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x){
      data.frame(
        subset(getScores(object$validation$temp_objects[[x]])$time, PC == component),
        model = x
      )
    })
    )
    dff$plotGroup <- paste0(dff$model,"-",dff$group)
    gs <- ggplot2::ggplot(dff, ggplot2::aes(x = time, y = score, group = plotGroup, color = group)) + 
      ggplot2::geom_point(alpha = 0.2) + ggplot2::geom_line(alpha = 0.2) +
      ggplot2::geom_point(data = subset(getScores(object)$time, PC == component), group = NA, alpha = 1, color = "black") +
      ggplot2::geom_line(data = subset(getScores(object)$time, PC == component), group = subset(getScores(object)$time, PC == component)$group, alpha = 1, color = "black") +
      ggplot2::labs(x = object$plot.xlabel) + myTheme + ggplot2::theme(legend.position = "bottom")
    
    # Loading plot
    dff <- Reduce(rbind, lapply(seq_along(object$validation$temp_objects), function(x){
      data.frame(
        subset(getLoadings(object$validation$temp_objects[[x]])$time, PC == component)
      )
    })
    )
    gl <- ggplot2::ggplot(dff, ggplot2::aes(x = covars, y = loading)) + 
      ggplot2::geom_point(alpha = 0.2, color = "black") + 
      ggplot2::geom_point(data = subset(getLoadings(object)$time, PC == component), alpha = 1, color = "red") +
      ggplot2::labs(x = object$plot.xlabel) + myTheme
    
    g <- ggpubr::ggarrange(gs, gl, nrow = 1, ncol = 2, widths = c(1,2))
  }
  return(g)
  
}

#' Plot covariate coefficients
#'
#' This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
#'
#' @param object An ALASCA object
#' @param covar Which covariable(s) to plot (default: `NA` which prints all)
#' @param tlab Alternative names for the covariables
#' @param return_data Set to `TRUE` to return data instead of plot
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot2 objects\.
#' 
#' @export
plotCovar <- function(object, covar = NA, tlab = NA, return_data = FALSE, myTheme = ggplot2::theme_bw(), pvalue = "shape"){
  df <- getCovars(object)
  if(any(is.na(covar))){
    covar <- unique(df$variable)
  }
  if(any(is.na(tlab))){
    tlab <- covar
  }
  for(i in 1:length(covar)){
    df$tlab[df$variable == covar[i]] <- tlab[i]
  }
  df$pvalue_label <- ifelse(df$pvalue >= 0.05, "Not significant", ifelse(df$pvalue < 0.001, "< 0.001", ifelse(df$pvalue < 0.01, "< 0.01", "< 0.05")))
  df$pvalue_sign <- ifelse(df$pvalue >= 0.05, "", ifelse(df$pvalue < 0.001, "***", ifelse(df$pvalue < 0.01, "**", "*")))
  df$covar <- factor(df$covar, levels = unique(df$covar[order(df$estimate)]))
  if(return_data){
    return(df)
  }else{
    if(pvalue == "shape"){
      g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, shape = pvalue_label)) + 
        ggplot2::geom_point() + ggplot2::geom_vline(xintercept = 0) +
        ggplot2::scale_shape_manual(values = c("Not significant"=3, "< 0.05"=15, "< 0.01"=16, "< 0.001"=17, "Baseline"=5)) +
        ggplot2::facet_wrap(~tlab) + ggplot2::labs(x = "Coefficient", y = "", shape = "P value") +
        myTheme + ggplot2::theme(legend.position = "bottom", legend.box="vertical", legend.margin=ggplot2::margin())
    }else if(pvalue == "star" | pvalue == "asterisk" | pvalue == "stars"){
      g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar, label = pvalue_sign)) + 
        ggplot2::geom_point() + ggplot2::geom_vline(xintercept = 0) +
        ggplot2::geom_text(hjust = 0.5, vjust = 0) +
        ggplot2::facet_wrap(~tlab) + ggplot2::labs(x = "Coefficient", y = "", shape = "P value") +
        myTheme
    }else{
      g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = covar)) + 
        ggplot2::geom_point() + ggplot2::geom_vline(xintercept = 0) +
        ggplot2::facet_wrap(~tlab) + ggplot2::labs(x = "Coefficient", y = "", shape = "P value") +
        myTheme + ggplot2::theme(legend.position = "bottom", legend.box="vertical", legend.margin=ggplot2::margin())
    }
    
    
    return(g)
  }
}


#' Plot projection of participants
#'
#' This function returns a plot of...
#'
#' @param object An ALASCA object
#' @param comp Which two components to plot (default: `c(1, 2`)
#' @param return_data Set to `TRUE` to return data instead of plot
#' @param myTheme A ggplot2 theme to use, defaults to `ggplot2::theme_bw()`
#' @return A ggplot2 objects.
#' 
#' @export
plotProjection <- function(object, comp = c(1,2), return_data = FALSE, myTheme = ggplot2::theme_bw()){
  df <- object$df
  df$ID <- df[, ID]
  loadings_Time <- subset(getLoadings(object)$time, PC %in% comp)
  loadings_Time <- reshape2::dcast(data = loadings_Time, covars ~ paste0("PC",PC), value.var = "loading")
  df_time <- merge(df, loadings_Time, by.x = "variable", by.y = "covars")
  if(object$separateTimeAndGroup){
    df_time <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x){
      data.frame(
        part = subset(df_time, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_time, paste0(ID, time) == x)$PC1 * subset(df_time, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_time, paste0(ID, time) == x)$PC2 * subset(df_time, paste0(ID, time) == x)$value),
        time = subset(df_time, paste0(ID, time) == x)$time
      )
    }))
    df_time <- df_time[!duplicated(df_time),]
    colnames(df_time) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time")
    if(!return_data){
      g_t <- ggplot2::ggplot(df_time, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part"))  + ggplot2::geom_point() + ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type="closed", length = ggplot2::unit(0.20,"cm"))) + myTheme
    }
    loadings_group <- subset(getLoadings(object)$group, PC %in% comp)
    loadings_group <- reshape2::dcast(data = loadings_group, covars ~ paste0("PC",PC), value.var = "loading")
    df_group <- merge(df, loadings_group, by.x = "variable", by.y = "covars")
    df_group <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x){
      data.frame(
        part = subset(df_group, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_group, paste0(ID, time) == x)$PC1 * subset(df_group, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_group, paste0(ID, time) == x)$PC2 * subset(df_group, paste0(ID, time) == x)$value),
        time = subset(df_group, paste0(ID, time) == x)$time,
        group = subset(df_group, paste0(ID, time) == x)$group
      )
    }))
    df_group <- df_group[!duplicated(df_group),]
    colnames(df_group) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time", "group")
    if(!return_data){
      g_g <- ggplot2::ggplot(df_group, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part", color = "group"))  + ggplot2::geom_point() + ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type="closed", length = ggplot2::unit(0.20,"cm"))) + myTheme
      g <- ggpubr::ggarrange(g_t, g_g, common.legend = TRUE, legend = "bottom")
    }
    g <- list(df_time, df_group)
    names(g) <- c("time", "group")
    return(g)
  }else{
    df_time <- Reduce(rbind, lapply(unique(paste0(df_time$ID, df_time$time)), function(x){
      data.frame(
        part = subset(df_time, paste0(ID, time) == x)$ID,
        pc1 = sum(subset(df_time, paste0(ID, time) == x)$PC1 * subset(df_time, paste0(ID, time) == x)$value),
        pc2 = sum(subset(df_time, paste0(ID, time) == x)$PC2 * subset(df_time, paste0(ID, time) == x)$value),
        time = subset(df_time, paste0(ID, time) == x)$time,
        group = subset(df_time, paste0(ID, time) == x)$group
      )
    }))
    df_time <- df_time[!duplicated(df_time),]
    colnames(df_time) <- c("part", paste0("PC", comp[1]), paste0("PC", comp[2]), "time", "group")
    g <- ggplot2::ggplot(df_time, ggplot2::aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), group = "part", color = "group"))  + 
      ggplot2::geom_point() + ggplot2::geom_line(alpha = 0.7, arrow = ggplot2::arrow(type="closed", length = ggplot2::unit(0.20,"cm"))) + myTheme
    if(return_data){
      return(df_time)
    }else{
      return(g)
    }
  }
  
}