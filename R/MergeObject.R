#' Merge a list of rds file Seurat object
#'
#' @inheritParams RunSeurat
#'
#' @param rds.list A vector of rds files paths
#' @param obj.list A vector of Seurat objects
#'
#' @return A merged Seurat Object
#'
#' @import Seurat
#' @import patchwork
#' @import sctransform
#' @import ggplot2
#' @export
#'
MergeObject <- function(rds.list, obj.list = NULL, output.dir = getwd(), save.rds = TRUE,
    ...) {
    if (is.null(obj.list)) {
        rds.list <- lapply(X = rds.list, FUN = function(x) {
            readRDS(x)
        })
    } else {
        rds.list <- obj.list
    }
    merged_obj <- merge(rds.list[[1]], rds.list[-1])

    if (save.rds) {
        saveRDS(merged_obj, paste0(output.dir, "/merged.obj.Rds"))
    }
    return(merged_obj)
}
