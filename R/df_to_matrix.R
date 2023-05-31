#' Convert a data.frame to a distance matrix
#' 
#' @description
#' Converts a `data.frame` to a distance `matrix`.
#' 
#' @param data a `data.frame` with the following three columns: `from` (the 
#'   first site), `to` (the second site) and `distance` (the value of the link
#'   between two sites). The output of the function [distance_along()].
#' 
#' @param lower a `logical` value. If `TRUE` (default), keep values in the 
#'   lower triangle of the matrix. Otherwise they will be replaced by `NA` 
#'   (or `0`).
#' 
#' @param upper a `logical` value. If `TRUE` (default), keep values in the 
#'   upper triangle of the matrix. Otherwise they will be replaced by `NA` 
#'   (or `0`).
#' 
#' @param diag a `logical` value. If `TRUE` (default), keep values in the 
#'   diagonal of the matrix. Otherwise they will be replaced by `NA` 
#'   (or `0`).
#'   
#' @param na_to_zero a `logical` value. If `TRUE` (default), missing edges are 
#'   coded as `0`. Otherwise they will be coded as `NA`.
#'
#' @return A matrix of dimensions `m x n`, where `m` is the number of sites 
#'   (`from`) and `n` is the number of sites (`to`).
#' 
#' @export
#'
#' @examples
#' # Add an example ----

df_to_matrix <- function(data, lower = TRUE, upper = TRUE, diag = TRUE, 
                         na_to_zero = TRUE) {
  
  # Check arguments ----
  
  check_dist_data(data)
  check_logical_value(lower)
  check_logical_value(upper)
  check_logical_value(diag)
  check_logical_value(na_to_zero)
  
  
  # Order sites ----
  
  sites_from <- sort(unique(data$"from"))
  sites_to   <- sort(unique(data$"to"))
  
  
  # Replace missing distances ----
  
  if (na_to_zero) {
    
    data$"distance" <- ifelse(is.na(data$"distance"), 0, data$"distance") 
    
  } else {
    
    data$"distance" <- ifelse(is.na(data$"distance") | data$"distance" == 0, 
                              NA, data$"distance")
  }
  
  
  # Pivot data frame ----
  
  mat <- tidyr::pivot_wider(data, names_from = "to", values_from = "distance", 
                            values_fn = ~.x)
  
  
  # Convert to matrix ----
  
  row_names <- mat[ , 1, drop = TRUE]
  mat <- data.matrix(mat[ , -1])
  rownames(mat) <- row_names
  
  
  # Order matrix ----
  
  mat <- mat[sites_from, sites_to]
  
  
  ## Apply filters ----
  
  if (!upper) mat[upper.tri(mat)] <- ifelse(na_to_zero, 0, NA)
  if (!lower) mat[lower.tri(mat)] <- ifelse(na_to_zero, 0, NA)
  if (!diag)  diag(mat)           <- ifelse(na_to_zero, 0, NA)
  
  mat
}
