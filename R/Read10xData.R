#' Create Seurat Object from sparse data matrices provided by 10X genomics
#'
#' This function create a Seurat object from a 10x output directory by loading
#' a matrix (matrix.mtx.gz), a barcode file (barcodes.tsv.gz) and a feature file
#' (features.tsv.gz) from the current directory.
#'
#'
#' @inheritParams Seurat::Read10X
#' @inheritParams RunSeurat
#' @inheritParams SeuratObject::CreateAssayObject
#'
#' @return A Seurat Object and its QC plots
#'
#' @importFrom patchwork wrap_plots
#' @import Seurat
#'
#' @export
Read10xData <- function(data.dir, output.dir = NULL, min.cells = 3, min.features = 200,
    mt.pattern = "^mt-", project.name = "SeuratProject", ...) {
    obj.data <- Seurat::Read10X(data.dir = data.dir)
    obj <- Seurat::CreateSeuratObject(counts = obj.data, project = project.name, min.cells = min.cells,
        min.features = min.features)

    message("Reading 10X")


    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt.pattern)

    return(obj)
}
