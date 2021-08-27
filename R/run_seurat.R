#' Run Seurat Pipeline
#'
#' This function run the following steps of the Seurat pipeline :
#' - Create Seurat Object
#' - QC
#' - Normalization
#' - Dimensional Reduction
#' - Clustering
#' - Visualization
#' - Differential Gene Expression
#'
#' @inheritParams read_10xdata
#' @inheritParams Seurat::Read10X
#' @inheritParams Seurat::CreateSeuratObject
#' @inheritParams Seurat::SCTransform
#' @inheritParams Seurat::RunPCA
#' @inheritParams Seurat::FindClusters
#' @inheritParams Seurat::RunUMAP
#' @inheritParams Seurat::DimPlot
#' @inheritParams Seurat::FindAllMarkers
#' @param output.dir Path to the output directory - where plots and markers file will be saved
#' @param save.rds Save the final Seurat object in RDS format (default set to True)
#' @param integrated.assay To run pipeline on the integrated assay
#'
#'
#' @return A processed Seurat Object along with QC and visualization plots
#'
#' @importFrom patchwork wrap_plots
#' @importFrom utils write.csv
#' @importFrom stringr str_to_title
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#' @export
#'
#' @import Seurat
#' @import clustree
#' @import MAST
#' @import sctransform

run_seurat <-
  function(data.dir = getwd(),
           output.dir = getwd(),
           min.cells = 3,
           min.features = 200,
           mt.pattern = '^mt-',
           project.name = 'SeuratProject',
           max.percent.mt = 15,
           max.features = NULL,
           max.nCount = NULL,
           sctransform = F,
           vars.to.regress = c('percent.mt', 'nCount_RNA'),
           ndims = 50,
           dims = 1:35,
           resolution = seq(0.1, 1, 0.1),
           cellcycle = TRUE,
           genes.FeaturePlot = NULL,
           genes.DotPlot = NULL,
           find.all.markers = F,
           only.pos = TRUE,
           min.pct = 0.25,
           logfc.threshold = 0.25,
           test.use = 'MAST',
           save.rds = TRUE,
           integrated.assay = FALSE,...) {

    options(warn = -1)

    obj <-
      read_10Xdata(
        data.dir = data.dir,
        output.dir = output.dir,
        min.cells = min.cells,
        min.features = min.features,
        mt.pattern = mt.pattern,
        project.name = project.name
      )

    message('QC filtering')

    cells_nb_pre_qc <- length(colnames(obj))

    if (is.null(max.nCount)) {
      max.nCount <- quantile(obj$nCount_RNA, 0.99, names = F)
    }
    if (is.null(max.features)) {
      max.features <- quantile(obj$nFeature_RNA, 0.99, names = F)
    }

    save_plot(output.dir = output.dir, filename = 'QC_percent.mt', Seurat::FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")+geom_hline(yintercept = max.percent.mt, linetype='dashed')+ggtitle(paste('QC plot \nFiltering Low-quality / dying cells \nMitochondrial counts over ', max.percent.mt,'%')))
    save_plot(output.dir = output.dir, filename = 'QC_nCount',Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+geom_hline(yintercept = max.features, linetype='dashed')+geom_vline(xintercept = max.nCount, linetype='dashed')+ggtitle(paste('QC plot \n','filter cells that have unique feature counts over', round(max.features),'and high gene count over ', round(max.nCount))))


    obj <-
      subset(
        obj,
        subset = nFeature_RNA > min.features &
          nFeature_RNA < max.features &
          percent.mt < max.percent.mt & nCount_RNA < max.nCount
      )

    cells_nb_post_qc <- length(colnames(obj))

    save_plot(
      output.dir=output.dir, 'QC_filtering',
      Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle(paste('QC plot after filtering (',cells_nb_pre_qc-cells_nb_post_qc,' cells removed)' ))
    )

    #NORMALIZATION
    if (sctransform) {
      message('Normalization method: SCTransform')
      obj <-
        Seurat::SCTransform(obj, vars.to.regress = vars.to.regress, verbose = F)
    }

    else {

      message('SCTransform set to FALSE')
      message('Normalization method: Log-Nomalization')
      obj <-
        Seurat::NormalizeData(
          obj,
          normalization.method = "LogNormalize",
          scale.factor = 10000,
          verbose = F
        )
      obj <-
        Seurat::FindVariableFeatures(
          obj,
          selection.method = "vst",
          nfeatures = 2000,
          verbose = F
        )
      obj <-
        Seurat::ScaleData(obj, features = rownames(obj), verbose = F)

    }

    #CLUSTERING
    message('Dimensional Reduction')
    obj <- Seurat::RunPCA(obj, verbose = FALSE)

    save_plot(output.dir = output.dir, filename = 'ElbowPlot', Seurat::ElbowPlot(obj, ndims = ndims))

    message('Clustering')
    obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = FALSE)
    obj <-
      Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)

    message('Running UMAP')
    obj <-
      Seurat::RunUMAP(obj, dims = dims, umap.method = "umap-learn", verbose=FALSE)
    #obj <- Seurat::FindClusters(obj, resolution = 0.3, verbose = FALSE)

    if (sctransform) {
      save_plot(output.dir = output.dir, filename = 'clustree', clustree(obj, prefix = 'SCT_snn_res.'))
    }
    else {
      save_plot(output.dir = output.dir, filename = 'clustree', clustree(obj, prefix = 'RNA_snn_res.'))
    }

    save_plot(
      output.dir = output.dir, filename = 'clustree_overlay',
      clustree_overlay(
        obj,
        red_dim = "umap",
        x_value = "umap1",
        y_value = "umap2"
      )
    )

    save_plot(output.dir=output.dir,filename = 'DimPlot_umap',
              Seurat::DimPlot(obj, reduction = 'umap', label = TRUE))
    save_plot(
      output.dir = output.dir, filename = 'Dimplot_umap_orig.ident',
      Seurat::DimPlot(obj, reduction = 'umap', group.by = 'orig.ident')
    )

    if (cellcycle) {
      message('Cell Cycling Scoring')
      s.genes <- str_to_title(tolower(cc.genes$s.genes))
      g2m.genes <- str_to_title(tolower(cc.genes$g2m.genes))
      obj <-
        Seurat::CellCycleScoring(
          obj,
          s.features = s.genes,
          g2m.features = g2m.genes,
          set.ident = FALSE
        )
      save_plot(output.dir=output.dir, filename = 'cellcycle',
                Seurat::DimPlot(obj, reduction = 'umap', label = TRUE, group.by = 'Phase'))
    }
    else{
      message('cellcycle set to FALSE')
    }

    if (!is.null(genes.FeaturePlot)) {
      save_plot(output.dir=output.dir, filename='FeaturePlot',
                Seurat::FeaturePlot(obj, features = genes.FeaturePlot))
    }
    else{
      message('No genes.FeaturePlot provided')
    }

    if (!is.null(genes.DotPlot)) {
      save_plot(output.dir=output.dir,filename='DotPlot',
                Seurat::DotPlot(obj, features = genes.DotPlot))
    }
    else{
      message('No genes.DotPlot provided')
    }

    if (find.all.markers) {
      message('Finding Markers')
      obj_markers <<-
        Seurat::FindAllMarkers(
          obj,
          only.pos = only.pos,
          min.pct = min.pct,
          logfc.threshold = logfc.threshold,
          test.use = test.use
        )
      write.csv(obj_markers, paste(output.dir, '/markers_obj.csv', sep=''))
    }
    else{
      message('find.all.markers set to FALSE')
    }

    if (save.rds){saveRDS(obj, file = output.dir)}

    message(paste('DONE ! All your plots and the final object are in :', output.dir))

    return(obj)

  }



