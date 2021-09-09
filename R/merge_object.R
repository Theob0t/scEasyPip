#' Merge a list of rds file Seurat object
#'
#'
#'
#' @inheritParams run_seurat
#'
#' @param rds.list A vector of rds files paths
#'
#'
#' @return A merged Seurat Object
#'
#' @import Seurat
#' @import patchwork
#' @import sctransform
#' @import ggplot2
#' @export
merge_object <- function(rds.list, output.dir=getwd(),save.rds=TRUE, ...){

  rds.list <- lapply(X = rds.list, FUN = function(x){readRDS(x)})

  merged_obj <- merge(rds.list[[1]], rds.list[-1])

  if(save.rds){saveRDS(merged_obj,paste0(output.dir,'/merged.obj.Rds'))}
  return(merged_obj)
}


