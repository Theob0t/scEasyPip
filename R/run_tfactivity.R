#' Create Dorothea Assay with TF activity computed using Viper
#'
#' @inheritParams read_10Xdata
#' @param object Path to the output directory - where plots and markers file will be saved
#' @param species Save the final Seurat object in RDS format (default set to True)
#' @param dorothea.confidence To run pipeline on the integrated assay
#'
#'
#' @return A processed Seurat Object along with visualization plots
#'
#' @importFrom utils data
#'
#' @import Seurat
#' @import dorothea
#' @import viper
#' @import dplyr
#'
#' @export

run_tfactivity <-
  function(data.dir = NULL,
           object = NULL,
           species = 'mouse',
           dorothea.confidence = c("A", "B", "C"),
           ...) {
    if (!is.null(data.dir)) {
      object <- read_10Xdata(data.dir = data.dir, ...)
    }

    message('Run Dorothea - creating regulons')
    if (species == 'mouse') {
      dorothea_regulon_mouse <-
        get(data("dorothea_mm", package = "dorothea"))
      regulon <- dorothea_regulon_mouse %>%
        dplyr::filter(confidence %in% dorothea.confidence)
    }
    if (species == 'human') {
      dorothea_regulon_hm <-
        get(data("dorothea_hs", package = "dorothea"))
      regulon <- dorothea_regulon_hm %>%
        dplyr::filter(confidence %in% dorothea.confidence)
    }
    else {
      message('Available species are human and mouse - Using default: mouse')
    }

    message('Run Viper - Computing TF activity score')
    object <- dorothea::run_viper(
      object,
      regulon,
      options = list(
        method = "scale",
        minsize = 4,
        eset.filter = FALSE,
        cores = 1,
        verbose = FALSE
      )
    )

    DefaultAssay(object) <- "dorothea"
    object@assays[["dorothea"]]@var.features <- rownames(object)
    object@assays[["dorothea"]]@scale.data <-
      object@assays[["dorothea"]]@data

    return(object)
  }



#' Run the Seurat Pipeline using TF activity (Dorothea + Viper)
#'
#' This function uses Dorothea and Viper to cluster single-cells by TF activity
#'
#'
#' @inheritParams run_tfactivity
#' @param output.dir Path to the output directory - where plots and markers file will be saved
#' @param save.rds Save the final Seurat object in RDS format (default set to True)
#' @param integrated.assay To run pipeline on the integrated assay
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
run_tfactivity_pipeline <-
  function(data.dir = NULL,
           output.dir = NULL,
           object = NULL,
           species = 'mouse',
           dorothea.confidence = c("A", "B", "C"),
           ...) {

    object <-
      run_tfactivity(
        data.dir = data.dir,
        object = object,
        species = species,
        dorothea.confidence = dorothea.confidence
      )

    DefaultAssay(object) <- 'dorothea'

    message('Run Seurat Pipeline')

    object <-
      run_seurat(
        object = object,
        output.dir = output.dir,
        tf.activity = TRUE,
        ...
      )

    return(object)
  }
