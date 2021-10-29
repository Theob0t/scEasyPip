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
#' @inheritParams Read10xData
#' @inheritParams RunTfPipeline
#' @inheritParams Seurat::Read10X
#' @inheritParams SeuratObject::CreateAssayObject
#' @inheritParams Seurat::SCTransform
#' @inheritParams Seurat::RunPCA
#' @inheritParams Seurat::FindClusters
#' @inheritParams Seurat::RunUMAP
#' @inheritParams Seurat::ElbowPlot
#' @inheritParams Seurat::DimPlot
#' @inheritParams Seurat::FindAllMarkers
#' @param object A Seurat object
#' @param output.dir Path to the destination folder of saved files
#' @param mt.pattern Regex pattern of the mitochondrial genes ('^MT-' or '^mt-')
#' @param save.rds Save final Seurat object in RDS format (default set to True)
#' @param integrated.assay If set, run the pipeline on the integrated assay
#' @param no.plot If set, run the pipeline without saving the plots
#' @param project.name Name of the project/object used for titles in plots
#' @param max.percent.mt Mitochondrial counts threshold (default set to 15)
#' @param max.features Maximum number of gene per cell (default 99-quantile)
#' @param max.nCount Maximum number of reads per cell (default 99-quantile)
#' @param sctransform If set, use SCTransform normalization
#' @param logtransform Run the default log-normalization from Seurat
#' @param cellcycle Run CellCycle Scoring from Seurat
#' @param genes.FeaturePlot A list of genes for Seurat FeaturePlot
#' @param genes.DotPlot A list of genes for Seurat DotPlot
#' @param find.all.markers Run Differential Gene Expression
#' @param tf.activity Run Seurat pipeline on TF activity
#'
#' @return A processed Seurat Object along with QC and visualization plots
#'
#' @importFrom patchwork wrap_plots
#' @importFrom utils write.csv
#' @importFrom stringr str_to_title
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#'
#'
#' @import Seurat
#' @import SeuratObject
#' @import umap
#' @import ggraph
#' @import clustree
#' @import MAST
#' @import sctransform
#'
#'
#' @export

RunSeurat <- function(data.dir = getwd(),
                      object = NULL,
                      output.dir = getwd(),
                      min.features = 200,
                      min.cells = 3,
                      project.name = "project_name",
                      mt.pattern = "^mt-",
                      max.percent.mt = 15,
                      max.features = NULL,
                      max.nCount = NULL,
                      sctransform = FALSE,
                      logtransform = TRUE,
                      no.plot = FALSE,
                      vars.to.regress = c("percent.mt", "nCount_RNA"),
                      ndims = 50,
                      dims = 1:35,
                      resolution = seq(0.1, 1, 0.1),
                      cellcycle = TRUE,
                      genes.FeaturePlot = NULL,
                      genes.DotPlot = NULL,
                      find.all.markers = TRUE,
                      only.pos = TRUE,
                      min.pct = 0.25,
                      logfc.threshold = 0.25,
                      test.use = "MAST",
                      save.rds = TRUE,
                      integrated.assay = FALSE,
                      tf.activity = FALSE,
                      species = "mouse",
                      dorothea.confidence = c("A", "B", "C"),
                      ...) {
    options(warn = -1)

    if (is.null(object)) {
        obj <- Read10xData(data.dir = data.dir, output.dir = output.dir, min.cells = min.cells,
            min.features = min.features, mt.pattern = mt.pattern, project.name = project.name,
            ...)
    } else {
        if (class(object)[1] != "Seurat") {
            stop("A 10x output directory or a Seurat object is needed")
        }
        obj <- object
    }

    if (tf.activity) {
        obj <- RunTfPipeline(object = obj, species = species, dorothea.confidence = dorothea.confidence)
    }

    message(DefaultAssay(obj))
    message("QC filtering")
    cells.nb.pre.qc <- length(colnames(obj))

    if (is.null(max.nCount)) {
        max.nCount <- quantile(obj$nCount_RNA, 0.99, names = FALSE)
    }
    if (is.null(max.features)) {
        max.features <- quantile(obj$nFeature_RNA, 0.99, names = FALSE)
    }

    SavePlot(output.dir = output.dir, filename = "QC_percent.mt", Seurat::FeatureScatter(obj,
        feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_hline(yintercept = max.percent.mt,
        linetype = "dashed") + ggtitle(paste("Filtering Low-quality / dying cells \nMitochondrial counts over ",
        max.percent.mt, "%")))

    SavePlot(output.dir = output.dir, filename = "QC_nCount", Seurat::FeatureScatter(obj,
        feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = max.features,
        linetype = "dashed") + geom_vline(xintercept = max.nCount, linetype = "dashed") +
        ggtitle(paste("QC plot \n", "filter out cells that have unique feature counts over",
            round(max.features), "and high gene count over", round(max.nCount))))

    obj <- subset(obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features &
        percent.mt < max.percent.mt & nCount_RNA < max.nCount)

    cells.nb.post.qc <- length(colnames(obj))

    SavePlot(output.dir = output.dir, "QC_filtering", Seurat::FeatureScatter(obj, feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA") + ggtitle(paste("QC plot after filtering (", cells.nb.pre.qc -
        cells.nb.post.qc, " cells removed)")))

    ## NORMALIZATION
    obj <- runNormalization(object = obj, sctransform = sctransform, logtransform = logtransform,
        tf.activity = tf.activity, vars.to.regress = vars.to.regress)
    ## CLUSTERING
    obj <- runClustering(object = obj, output.dir = output.dir, ndims = ndims, dims = dims,
        resolution = resolution, no.plot = no.plot, sctransform = sctransform)

    ## CellCycle
    if (cellcycle) {
        object <- runCellCycle(object = object, output.dir = output.dir, no.plot = no.plot, sctransform = sctransform)
    } else {
        message("cellcycle set to FALSE")
    }

    # Feature Plots
    if (!is.null(genes.FeaturePlot)) {
        for (g in genes.FeaturePlot) {
            SavePlot(output.dir = output.dir, filename = paste0(paste0("FeaturePlot_",
                g), ".tiff"), Seurat::FeaturePlot(obj, features = g))
        }
    } else {
        message("No genes.FeaturePlot provided")
    }

    # DotPlot
    if (!is.null(genes.DotPlot)) {
        SavePlot(output.dir = output.dir, filename = "DotPlot", Seurat::DotPlot(obj,
            features = genes.DotPlot))
    } else {
        message("No genes.DotPlot provided")
    }

    ## DGE
    if (find.all.markers) {
        message("Finding Markers")
        obj_markers <- Seurat::FindAllMarkers(obj, only.pos = only.pos, min.pct = min.pct,
            logfc.threshold = logfc.threshold, test.use = test.use, ...)
        write.csv(obj_markers, paste0(output.dir, "/markers_obj.csv"))
    } else {
        message("find.all.markers set to FALSE")
    }

    ## Saving final object
    if (save.rds) {
        message("saving rds file ...")
        saveRDS(obj, file = paste0(output.dir, "/seurat.obj.Rds"))
    }

    message(paste("DONE ! Plots and final object are in :", output.dir))

    return(obj)
}

#' Run Normalization using Seurat pipeline
#'
#' @inheritParams RunSeurat
#' @inheritParams Seurat::SCTransform
#' @inheritParams Seurat::RunPCA
#'
#' @return A normalized Seurat Object
#'
#' @import Seurat
#' @import sctransform
#'
#' @export
runNormalization <- function(object, sctransform, logtransform, tf.activity, vars.to.regress) {
    if (!tf.activity) {
        if (sctransform) {
            message("Normalization method: SCTransform")
            object <- Seurat::SCTransform(object, vars.to.regress = vars.to.regress,
                verbose = FALSE)
        }
        if (logtransform) {
            message("Normalization method: Log-Nomalization")
            object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize",
                scale.factor = 10000, verbose = FALSE, assay = "RNA")
            message("Finding Variable Features")
            object <- Seurat::FindVariableFeatures(object, selection.method = "vst",
                nfeatures = 3000, verbose = FALSE, assay = "RNA")
            message("Scaling Data")
            object <- Seurat::ScaleData(object, features = rownames(object), verbose = FALSE,
                assay = "RNA")

        }
        if (sctransform) {
            DefaultAssay(object) <- "SCT"
        }
    }
    return(object)

}


#' Perform clustering on Seurat object
#'
#'
#' @inheritParams RunSeurat
#' @inheritParams Seurat::RunPCA
#' @inheritParams Seurat::FindClusters
#' @inheritParams Seurat::RunUMAP
#'
#'
#' @return A Seurat object with cells clustered and some cluster visualizations
#'
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#'
#' @import Seurat
#' @import clustree
#'
#' @export
runClustering <- function(object, output.dir, ndims, dims, resolution, no.plot, sctransform) {
    message("Dimensional Reduction")
    object <- Seurat::RunPCA(object, verbose = FALSE)

    SavePlot(output.dir = output.dir, filename = "ElbowPlot", Seurat::ElbowPlot(object,
        ndims = ndims) + geom_vline(xintercept = max(dims)))

    message("Clustering")
    object <- Seurat::FindNeighbors(object, dims = dims, verbose = FALSE)
    object <- Seurat::FindClusters(object, resolution = resolution, verbose = FALSE)

    message("Running UMAP")
    object <- Seurat::RunUMAP(object, dims = dims, umap.method = "umap-learn", verbose = FALSE)

    if (!no.plot) {
        if (sctransform) {
            SavePlot(output.dir = output.dir, filename = "clustree", clustree::clustree(object,
                prefix = "SCT_snn_res."))
        }
        if (TRUE %in% grepl("^RNA", x = colnames(object@meta.data))) {
            SavePlot(output.dir = output.dir, filename = "clustree", clustree::clustree(object,
                prefix = "RNA_snn_res."))
        }
        if (TRUE %in% grepl("^dorothea", x = colnames(object@meta.data))) {
            SavePlot(output.dir = output.dir, filename = "clustree", clustree::clustree(object,
                prefix = "dorothea_snn_res."))
        }

        SavePlot(output.dir = output.dir, filename = "clustree_overlay", clustree_overlay(object,
            red_dim = "umap", x_value = "umap1", y_value = "umap2"))

        SavePlot(output.dir = output.dir, filename = "DimPlot_umap", Seurat::DimPlot(object,
            reduction = "umap", label = TRUE))
        SavePlot(output.dir = output.dir, filename = "Dimplot_umap_orig.ident", Seurat::DimPlot(object,
            reduction = "umap", group.by = "orig.ident"))

        SavePlot(output.dir = output.dir, filename = "Dimplot_res_0.5", Seurat::DimPlot(FindClusters(object,
            resolution = 0.5, verbose = FALSE), reduction = "umap") + ggtitle("Cluster resolution=0.5"))
    }
    return(object)
}


#' Perform cell cycle scoring
#'
#' @inheritParams RunSeurat
#'
#' @return A Seurat object with Cell Cycle phase for each cells
#'
#' @importFrom patchwork wrap_plots
#' @importFrom utils write.csv
#' @importFrom stringr str_to_title
#' @importFrom grDevices dev.off
#'
#' @import Seurat
#'
#' @export
runCellCycle <- function(object, output.dir, no.plot, sctransform) {
    assay <- DefaultAssay(object)
    if (sctransform) {
        DefaultAssay(object) <- "SCT"
    } else {
        DefaultAssay(object) <- "RNA"
    }

    message("Cell Cycling Scoring")
    s.genes <- str_to_title(tolower(cc.genes$s.genes))
    g2m.genes <- str_to_title(tolower(cc.genes$g2m.genes))
    object <- Seurat::CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes,
        set.ident = FALSE)

    if (!no.plot) {
        SavePlot(output.dir = output.dir, filename = "cellcycle", Seurat::DimPlot(object,
            reduction = "umap", label = TRUE, group.by = "Phase") + ggtitle("Cell Cycle Phases"))
    }
    DefaultAssay(object) <- assay
    return(object)
}
