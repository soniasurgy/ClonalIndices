#' Print method for index objects
#'
#' Prints the predicted genetic gains and selected genotypes.
#'
#' @param x An object of class \code{index}
#' @param ... Additional arguments passed to methods (currently unused)
#' @return Invisibly returns the input object \code{x}
#' @rdname index
#' @keywords internal
#' @export
print.index <- function(x, ...) {
  cat("Rank of genotypes for index", attr(x, "indice"), "\n")
  print(x$ranking)
  invisible(x)

}
