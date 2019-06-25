#' @useDynLib MeanShiftR
.onUnload <- function (libpath) {
  library.dynam.unload("MeanShiftR", libpath)

  invisible()
}
