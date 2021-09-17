#' Run Seurat Anchors Integration Pipeline
#'
#' This function run the following steps of the Seurat Anchors Integration pipeline :
#' - Create Seurat Object
#' - QC
#' - Normalization
#' - Dimensional Reduction
#' - Clustering
#' - Visualization
#' - Differential Gene Expression
#'
#' @inheritParams run_seurat
#'
#' @param rds.list A vector of rds files paths
#'
#'
#' @return An integrated Seurat Object along with its visualization plots
#'
#' @importFrom utils write.csv
#' @importFrom stringr str_to_title
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#'
#' @import Seurat
#' @import patchwork
#' @import sctransform
#' @import ggplot2
#' @export
run_integration <- function(rds.list, output.dir=getwd(), resolution=0.5, grouping.var = 'orig.ident', find.conserved.markers = FALSE, ident.1 = NULL,
                            ident.2 = NULL, save.rds=TRUE, ...){

  rds.list <- lapply(X = rds.list, FUN = function(x){readRDS(x)})

  if ('SCT' %in% Assays(rds.list[[1]])){
    message('SCT array found')
    features <- SelectIntegrationFeatures(object.list = rds.list, nfeatures = 3000)
    obj.list <- PrepSCTIntegration(object.list = rds.list, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = "rpca", normalization.method = "SCT")
    obj.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

    obj.integrated <- RunPCA(obj.integrated, verbose = FALSE)
    obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:30)

    save_plot(
      output.dir = output.dir, filename = 'integrated_sct_umap',
      Seurat::DimPlot(obj.integrated, reduction = 'umap', group.by = grouping.var)+ggtitle(paste('Integrated object \n', grouping.var))
    )
  }

  else {
    message('No SCT array found')

    obj.list <- lapply(X = rds.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

    print('features')
    features <- SelectIntegrationFeatures(object.list = obj.list)

    obj.list <- lapply(X = rds.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })

    print('anchors')
    anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = "rpca")
    print('integrate')
    obj.integrated <- IntegrateData(anchorset = anchors)

    # Run the standard workflow for visualization and clustering
    obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
    obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
    obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:30)

    save_plot(
      output.dir = output.dir, filename = 'integrated_umap',
      Seurat::DimPlot(obj.integrated, reduction = 'umap', group.by = grouping.var)+ggtitle(paste('Integrated object \n', grouping.var))
    )

  }

  if(find.conserved.markers){
  DefaultAssay(obj.integrated) <- "RNA"
  write.csv(FindConservedMarkers(obj.integrated, ident.1 = ident.1, ident.2 = ,grouping.var = grouping.var, verbose = FALSE), paste(output.dir, '/conserved_markers.csv', sep=''))
  }

  if (save.rds){
    print('saving rds')
    saveRDS(object = obj.integrated, file =paste0(output.dir,'obj.integrated.Rds'))
    }

  message(paste('DONE ! All your plots and the final integrated object are in :', output.dir))

  return(obj.integrated)

}
