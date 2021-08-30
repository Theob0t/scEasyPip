#' Create Seurat Object from sparse data matrices provided by 10X genomics
#'
#' This function create a Seurat object from a 10x output directory by loading
#' a matrix (matrix.mtx.gz), a barcode file (barcodes.tsv.gz) and a feature file
#' (features.tsv.gz) from the current directory.
#'
#'
#' @inheritParams Seurat::Read10X
#' @inheritParams Seurat::CreateSeuratObject
#' @param output.dir Path to the output directory - where plots and markers file will be saved
#' @param mt.pattern The regex pattern of the mitochondrial genes ('^MT-' or '^mt-')
#'
#' @return A Seurat Object and its QC plots
#'
#' @importFrom patchwork wrap_plots
#' @import Seurat
#'
#' @export
read_10Xdata <-
  function(data.dir,
           output.dir = NULL,
           min.cells = 3,
           min.features = 200,
           mt.pattern = '^mt-',
           project.name = 'SeuratProject',...) {

    obj.data <- Seurat::Read10X(data.dir = data.dir)
    obj <-
      Seurat::CreateSeuratObject(
        counts = obj.data,
        project = project.name,
        min.cells = min.cells,
        min.features = min.features
      )

    message('reading 10X')


    obj[["percent.mt"]] <-
      Seurat::PercentageFeatureSet(obj, pattern = mt.pattern)

    return(obj)
  }
