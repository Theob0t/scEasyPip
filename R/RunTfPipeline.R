#' Run the Seurat Pipeline using TF activity (Dorothea + Viper)
#'
#' This function uses Dorothea and Viper to cluster single-cells by TF activity
#'
#'
#' @inheritParams RunSeurat
#' @param species Mouse or human
#' @param dorothea.confidence Confidence levels of regulons
#'
#'
#' @return A processed Seurat Object along with visualization plots
#'
#' @importFrom stats quantile
#'
#' @import Seurat
#' @import clustree
#' @import dorothea
#' @import viper
#'
#' @export
RunTfPipeline <- function(data.dir = NULL, object = NULL, species = "mouse", dorothea.confidence = c("A",
    "B", "C"), ...) {
    object <- RunTf(data.dir = data.dir, object = object, species = species, dorothea.confidence = dorothea.confidence)

    DefaultAssay(object) <- "dorothea"

    message("TF Activity scores computed !")

    return(object)
}


#' Create Dorothea Assay with TF activity computed using Viper
#'
#' @inheritParams Read10xData
#' @inheritParams RunSeurat
#' @inheritParams RunTfPipeline
#'
#' @return A processed Seurat Object along with visualization plots based on TF activity
#'
#' @importFrom utils data
#'
#' @import Seurat
#' @import dorothea
#' @import viper
#' @import dplyr
#'
#' @export

RunTf <- function(data.dir = NULL, object = NULL, species = "mouse", dorothea.confidence = c("A",
    "B", "C"), ...) {
    if (!is.null(data.dir)) {
        object <- Read10xData(data.dir = data.dir, ...)
    }

    message("Running Dorothea - creating regulons")
    if (species == "human") {
        message("Using : human regulons")
        dorothea_regulon_hm <- get(data("dorothea_hs", package = "dorothea", envir = environment()))
        regulon <- dorothea_regulon_hm %>%
            dplyr::filter(confidence %in% dorothea.confidence)
    } else {
        message("Using : mouse regulons")
        dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea", envir = environment()))
        regulon <- dorothea_regulon_mouse %>%
            dplyr::filter(confidence %in% dorothea.confidence)
    }

    message("Run Viper - Computing TF activity scores")
    object <- dorothea::run_viper(object, regulon, options = list(method = "scale",
        minsize = 4, eset.filter = FALSE, cores = 1, verbose = FALSE))

    DefaultAssay(object) <- "dorothea"
    object@assays[["dorothea"]]@var.features <- rownames(object)
    object@assays[["dorothea"]]@scale.data <- object@assays[["dorothea"]]@data

    return(object)
}
