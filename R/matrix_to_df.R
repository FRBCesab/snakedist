#' Convert a distance matrix to a data.frame
#' 
#' @description
#' Converts a (distance) `matrix` to a `data.frame`.
#' 
#' @param x a `matrix`.
#' 
#' @param all a `logical` value. If `TRUE` (default), keep cells with `NA`.
#'
#' @return A `data.frame` with the following three columns: `from` (the 
#'   first site), `to` (the second site) and `distance` (the value of the link
#'   between two sites).
#' 
#' @export
#'
#' @examples
#' # Add an example ----

matrix_to_df <- function(x, all = FALSE) {
  
  ## Check 'x' argument ----
  
  if (missing(x)) {
    stop("Argument 'x' is required", call. = FALSE)
  }
  
  if (!is.matrix(x)) {
    stop("Argument 'x' must be a matrix", call. = FALSE)
  }
  
  if (!is.numeric(x)) {
    stop("Argument 'x' must be a numeric matrix", 
         call. = FALSE)
  }
  
  if (is.null(rownames(x))) {
    stop("Row names of 'x' must contain nodes labels", call. = FALSE)
  }

  
  ## Prepare for pivot ----
  
  x <- as.data.frame(x)
  x <- data.frame("from" = rownames(x), x)
  rownames(x) <- NULL
  colnames(x)[-1] <- x$"from"
  
  
  ## Pivot to longer format ----
  
  x <- tidyr::pivot_longer(x, cols = -1, names_to = "to", 
                           values_to = "distance")
  
  x <- as.data.frame(x)
  
  
  ## Remove 0/NA ----
  
  if (!all) {
    x <- x[which(!is.na(x$"distance") & x$"distance" != 0), ]
  }
  
  rownames(x) <- NULL
  x
}
