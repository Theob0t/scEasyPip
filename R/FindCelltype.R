#' Identify cell types based on a user defined consensus markers
#'
#' Assign cell identities to clusters
#'
#' @inheritParams RunSeurat
#' @inheritParams Seurat::FindAllMarkers
#' @param markers.file Path to the markers file (FindAllMarkers() output)
#' @param consensus.file Path to the consensus file
#' @param label Boolean to display cluster labels on umap plots
#' @param filename Character string to name the saved umap plots
#'
#' @return A Seurat Object with cells identities assigned to clusters
#'
#' @importFrom patchwork wrap_plots
#' @importFrom utils write.csv
#' @importFrom utils read.csv
#' @importFrom stringr str_to_title
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#'
#'
#' @import Seurat
#' @import MAST
#' @import dplyr
#'
#' @export

FindCelltype <- function(markers.file, consensus.file, object, output.dir, filename = NULL,
    save.rds = TRUE, label = TRUE, ...) {
    ## 1 - Check inputs
    inputs <- checkInputs(markers.file = markers.file, consensus.file = consensus.file,
        object = object, filename = filename)

    markers <- inputs[[1]]
    consensus <- inputs[[2]]

    ## 2 - Run assignment function
    vectors <- getAssignmentsVectors(consensus = consensus, markers = markers)

    cluster.ident.fc <- vectors[[1]]
    cluster.ident.nb <- vectors[[2]]

    ## 3 - Save results
    names(cluster.ident.fc) <- levels(object)

    object$assignment_fc <- Idents(RenameIdents(object, cluster.ident.fc))

    names(cluster.ident.nb) <- levels(object)

    object$assignment_nb <- Idents(RenameIdents(object, cluster.ident.nb))

    SavePlot(output.dir = output.dir, filename = paste0("Dimplot_", filename), plot = Seurat::DimPlot(object,
        reduction = "umap", label = label) + ggtitle("Cluster Identity"))

    SavePlot(output.dir = output.dir, filename = paste0("DimPlot_fc_", filename), plot = Seurat::DimPlot(object,
        reduction = "umap", group.by = "assignment_fc", label = label) + ggtitle("Cluster assignment (FoldChange)"))

    SavePlot(output.dir = output.dir, filename = paste0("DimPlot_nb_", filename), Seurat::DimPlot(object,
        reduction = "umap", group.by = "assignment_nb", label = label) + ggtitle("Cluster assignment (# of markers)"))

    if (save.rds) {
        saveRDS(object, paste0(output.dir, "/assigned_object.Rds"))
    }

    return(object)
}


#' Check inputs for FindCelltypes function
#'
#' @inheritParams FindCelltype
#'
#' @return Return markers and consensus files as dataframes
#'
#' @importFrom utils write.csv
#' @importFrom utils read.csv
#'
#' @import Seurat
#' @import dplyr
#' @export
checkInputs <- function(markers.file, consensus.file, object, filename) {
    if (is.null(filename)) {
        filename <- object@project.name
    }

    ## check if data has been logtransformed
    if (dim(object@assays[["RNA"]]@scale.data)[1] == 0) {
        stop("RNA assay not Log-Normalized")
    }

    ## check if clustering has been ran
    if (is.null(object@meta.data[["seurat_clusters"]])) {
        stop("No cluster in object")
    }

    ## read files
    markers <- read.csv(markers.file)
    consensus <- read.csv(consensus.file)

    return(list(markers, consensus))
}



#' Assign clusters to cell identities from the consensus file provided
#'
#' @param markers Markers dataframe (FindAllMarkers() output)
#' @param consensus Consensus dataframe (2 columns per cell-type labelled
#' cell-type e.g., 'Stem-Cell', and FC_cell-type e.g., 'FC_Stem-Cell')
#' @return Return two vectors with cell assignments.
#' One using the max number of overlap, and one using max total FoldChange.
#'
#' @importFrom utils write.csv
#' @importFrom utils read.csv
#'
#' @import Seurat
#' @import dplyr
#' @export
getAssignmentsVectors <- function(consensus, markers) {
    celltypes <- colnames(consensus %>%
        select(-contains("FC_")))

    cluster.ident.fc <- vector()
    cluster.ident.nb <- vector()

    ## Get marker genes and scores of each cluster
    for (c in 0:max(markers$cluster)) {
        cluster <- subset(markers, cluster == c) %>%
            select(c("avg_log2FC", "gene"))

        fc.scores <- vector()
        nb.scores <- vector()

        ## Get intersection between cell-type markers and cluster markers
        i <- 0
        for (ctype in celltypes) {
            commun.genes <- intersect(cluster[, "gene"], consensus[, ctype])


            fc.consensus <- consensus[consensus[, ctype] %in% commun.genes,
                                      paste0("FC_", ctype)]

            nb.scores[i] <- length(fc.consensus)
            fc.scores[i] <- round(sum(fc.consensus), 1)

            i <- i + 1
        }
        cluster.ident.fc[c] <- celltypes[which.max(fc.scores)]
        cluster.ident.nb[c] <- celltypes[which.max(nb.scores)]
    }
    return(list(cluster.ident.fc, cluster.ident.nb))
}
