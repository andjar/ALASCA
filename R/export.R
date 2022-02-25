#' Save ALASCA object
#'
#' @inheritParams saveALASCAModel
#' @inheritParams savetocsv
#' @return An ALASCA object
#' @export
saveALASCA <- function(object, filename = NA, filepath = NA, saveCSV = TRUE, saveScores = TRUE, saveLoadings = TRUE, saveCovars = TRUE, csv = "csv", ...) {
  saveALASCAModel(object = object, filename = filename, filepath = NA)
  savetocsv(object = object, filename = filename, filepath = filepath, saveCSV = saveCSV, saveScores = saveScores, saveLoadings = saveLoadings, saveCovars = saveCovars, csv = "csv", ...)
  summary.ALASCA(object = object, file = getFilename(object = object, filetype = "txt"), sessioninfo = TRUE)
}

#' Save ALASCA object
#'
#' @inheritParams saveALASCAModel
#' @inheritParams savetocsv
#' @return An ALASCA object
#' @export
saveBootstrapID <- function(object) {
  if (!(object$validate && object$validationMethod == "bootstrap")) stop("Please validate with bootstrapping")
  for (i in seq_along(object$validation$temp_object)) {
    write(paste0(object$validation$temp_object[[i]]$originalIDs, collapse = ";"),
      file = getFilename(object = object, prefix = "bootstrapID_", filetype = ".csv", overwrite = TRUE), append = TRUE
    )
  }
}

#' Save ALASCA object
#'
#' @inheritParams saveALASCAModel
#' @inheritParams savetocsv
#' @return An ALASCA object
#' @export
saveJackKnifeID <- function(object) {
  if (!(object$validate && object$validationMethod %in% c("loo", "jack-knife", "jackknife"))) stop("Please validate with jack-knife")
  for (i in seq_along(object$validation$temp_object)) {
    write(paste0(unique(object$validation$temp_object[[i]]$partID), collapse = ";"),
          file = getFilename(object = object, prefix = "jackknifeID_", filetype = ".csv", overwrite = TRUE), append = TRUE
    )
  }
}

#' Save ALASCA object to csv
#'
#' @param object An ALASCA object
#' @param filepath
#' @param filename
#' @export
savetocsv <- function(object, filename = NA, filepath = NA, saveCSV = TRUE, saveScores = TRUE, saveLoadings = TRUE, saveCovars = TRUE, csv = "csv", ...) {
  if (saveCSV) {
    if (csv == "csv") {
      if (object$separateTimeAndGroup) {
        if (saveScores) {
          write.csv(getScores(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_time_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv(getLoadings(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_time_loadings", filetype = ".csv"), ...
          )
        }
        if (saveScores) {
          write.csv(getScores(object)$group,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_group_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv(getLoadings(object)$group,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_group_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv(getCovars(object),
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      } else {
        if (saveScores) {
          write.csv(getScores(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv(getLoadings(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv(getCovars(object),
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      }
    } else if (csv == "csv2") {
      if (object$separateTimeAndGroup) {
        if (saveScores) {
          write.csv2(getScores(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_time_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv2(getLoadings(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_time_loadings", filetype = ".csv"), ...
          )
        }
        if (saveScores) {
          write.csv2(getScores(object)$group,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_scores_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv2(getLoadings(object)$group,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_scores_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv2(getCovars(object),
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      } else {
        if (saveScores) {
          write.csv2(getScores(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_scores", filetype = ".csv"), ...
          )
        }
        if (saveLoadings) {
          write.csv2(getLoadings(object)$time,
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_loadings", filetype = ".csv"), ...
          )
        }
        if (saveCovars) {
          write.csv2(getCovars(object),
            file = getFilename(object = object, filename = filename, filepath = filepath, suffix = "_covars", filetype = ".csv"), ...
          )
        }
      }
    } else {
      stop("Unkown input to `csv`")
    }
    cat(paste0("- Saved csv's to ", getFilename(object = object, filename = filename, filepath = filepath), "\n"))
  }
}

#' Save ALASCA object
#'
#' @param object An ALASCA object
#' @param filepath
#' @param filename
#' @return Full file name of the saved object (String)
#' @export
saveALASCAModel <- function(object, filename = NA, filepath = NA) {
  fname <- getFilename(object = object, filename = filename, filepath = filepath, filetype = ".rds")
  saveRDS(object, file = fname)
  cat(paste0("- Saved model to ", fname, "\n"))
  return(fname)
}

#' Save figure
#'
#' @param g ggplot-object
#' @return A ggplot2 objects.
#'
#' @export
saveALASCAPlot <- function(object, g, filetype = NA, figsize = NA, prefix = "plot/", suffix = "", figunit = NA) {
  if (any(is.na(filetype))) {
    filetype <- object$plot.filetype
  }
  if (any(is.na(figsize))) {
    figsize <- object$plot.figsize
  }
  if (is.na(figunit)) {
    figunit <- object$plot.figunit
  }

  for (i in filetype) {
    fname <- getFilename(object = object, prefix = prefix, suffix = suffix, filetype = i)
    ggplot2::ggsave(plot = g, filename = fname, width = figsize[1], height = figsize[2], dpi = figsize[3], unit = figunit)
    cat(paste0("- Saved ", fname, "\n"))
  }
}

#' Get file path
#'
#' @param object An ALASCA object
#' @param filename File name
#' @param filepath File path
#' @param prefix Prefix
#' @param suffix Suffix
#'
#' @return A ggplot2 objects.
#'
#' @export
getFilename <- function(object, filename = NA, filepath = NA, prefix = "", suffix = "", filetype = "", overwrite = FALSE) {
  # Use arguments if defined
  if (any(!is.na(filename))) {
    object$filename <- filename
  }
  if (any(!is.na(filepath))) {
    object$filepath <- filepath
  }

  # Set defaults if undefined
  if (any(is.na(object$filepath))) {
    object$filepath <- paste0("ALASCA/", strftime(object$initTime, format = "%Y%m%d_%H%M%S"), "/")
  }
  if (any(is.na(object$filename))) {
    object$filename <- "ALASCA"
  }

  # Create folder
  if (!dir.exists(object$filepath)) {
    dir.create(object$filepath, recursive = TRUE)
  }

  # Check if prefix implies subfolder
  if (grepl("/", prefix)) {
    if (!dir.exists(paste0(object$filepath, prefix))) {
      dir.create(paste0(object$filepath, prefix), recursive = TRUE)
    }
  }

  fname <- paste0(object$filepath, prefix, object$filename)

  # Check if file already exists
  if (!overwrite) {
    cnt <- 1
    while (file.exists(paste0(fname, suffix, ifelse(substr(filetype, 1, 1) == ".", filetype, paste0(".", filetype))))) {
      fname <- paste0(object$filepath, prefix, object$filename, "_", cnt)
      cnt <- cnt + 1
    }
  }
  fname <- paste0(fname, suffix, ifelse(substr(filetype, 1, 1) == "." | filetype == "", filetype, paste0(".", filetype)))
  return(fname)
}
