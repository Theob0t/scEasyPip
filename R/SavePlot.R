#' Save output plot in png format
#'
#' @param filename The filename under which to save the plot
#' @param output.dir Path to the destination folder of saved files
#' @param plot The plot to save
#'
#' @importFrom grDevices png
#'
#' @return  The plot saved as a png image
#'
#' @export
SavePlot <- function(filename, plot, output.dir) {
    png(filename = paste(output.dir, "/", filename, ".png", sep = ""), width = 720,
        height = 576)
    print(plot)
    dev.off()

}
