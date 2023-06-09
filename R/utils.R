#' Convert a LINESTRING to POINT
#' 
#' @description
#' Converts a spatial line of type `LINESTRING` into `POINT` by sampling points
#' on a linear geometry using the function [sf::st_line_sample()]. This 
#' function is used by the function [distance_along()].
#' 
#' @param x an `sf` object of type `LINESTRING`.
#' 
#' @inheritParams distance_along
#' 
#' @noRd
#' 
#' @examples
#' \dontrun{
#' # Import Adour river (LINESTRING) ----
#' path_to_file <- system.file("extdata", "adour_lambert93.gpkg", 
#'                             package = "chessboard")
#' adour_river <- sf::st_read(path_to_file)
#' 
#' # Sample points of the Adour river ----
#' adour_river_pts <- line_to_points(adour_river)
#' }

line_to_points <- function(x, density = 0.01, type = "regular", ...) {
  
  ## Check 'x' argument ----
  
  if (missing(x)) {
    stop("Argument 'x' (spatial layer of linear shape) is required", 
         call. = FALSE)
  }
  
  if (!inherits(x, "sf")) {
    stop("The object 'x' must be an 'sf' object", 
         call. = FALSE)
  }
  
  if (nrow(x) != 1) {
    stop("Argument 'x' should have exactly one row", 
         call. = FALSE)
  }
  
  geom <- sf::st_geometry_type(x) |> as.character() |> unique()
  
  if (!("LINESTRING" %in% geom)) {
    stop("Linear shape geometry must be of type LINESTRING", call. = FALSE)
  }
  
  
  ## Check 'type' argument ----
  
  type <- tolower(type)
  
  if (!(type %in% c("regular", "random"))) {
    stop("Argument 'type' must either 'regular' or 'random'")
  }
  
  
  ## Check density argument ----
  
  if (!is.numeric(density) || length(density) != 1) {
    stop("Argument 'density' must be a numeric of length 1", call. = FALSE)
  }
  
  if (density <= 0) {
    stop("Argument 'density' must be > 0", call. = FALSE)
  }
  
  
  ## Convert LINESTRING to POINT ----
  
  sampled_points <- sf::st_line_sample(x, density = density, type = type, ...)
  
  if (sf::st_is_empty(sampled_points)) {
    stop("Unable to sample points along the linear shape. Please increase ", 
         "the value of 'density'", call. = FALSE)
  }
  
  sampled_points <- sampled_points |> 
    sf::st_sf() |> 
    sf::st_cast("POINT") |> 
    dplyr::mutate(id = 1:dplyr::n(), group = 1) |>
    dplyr::select(.data$id, .data$group)
}



#' Convert POINT to a LINESTRING
#' 
#' @description
#' Converts spatial points of type `POINT` into a `LINESTRING`. This function
#' is used by the function [distance_along()] and is the opposite of 
#' [line_to_points()].
#' 
#' @param points_sf an `sf` object of type `POINT`. Typically points sampled on
#'   a LINESTRING as the output of [line_to_points()].
#' 
#' @param from an `integer` of length 1. The start of the segment (row number 
#'   in `points_sf`).
#'   
#' @param to an `integer` of length 1. The end of the segment (row number in 
#'   `points_sf`).
#' 
#' @noRd
#' 
#' @examples
#' \dontrun{
#' # Import Adour river (LINESTRING) ----
#' path_to_file <- system.file("extdata", "adour_lambert93.gpkg", 
#'                             package = "chessboard")
#' adour_river <- sf::st_read(path_to_file)
#' 
#' # Sample points of the Adour river ----
#' adour_river_pts <- line_to_points(adour_river)
#' 
#' # Convert back points to LINESTRING ----
#' adour_river_all <- points_to_line(adour_river_pts, 1, nrow(adour_river_pts))
#' 
#' # Extract on segment on the river ----
#' adour_river_seg <- points_to_line(adour_river_pts, 100, 350)
#' 
#' # Map ----
#' plot(sf::st_geometry(adour_river))
#' plot(sf::st_geometry(adour_river_all), add = TRUE, col = "red")
#' plot(sf::st_geometry(adour_river_seg), add = TRUE, col = "blue", lwd = 2)
#' }

points_to_line <- function(points_sf, from, to) {
  
  ## Check 'points_sf' argument ----
  
  if (missing(points_sf)) {
    stop("Argument 'points_sf' is required", call. = FALSE)
  }
  
  if (!inherits(points_sf, "sf")) {
    stop("The object 'points_sf' must be an 'sf' object", call. = FALSE)
  }
  
  if (nrow(points_sf) < 2) {
    stop("Argument 'points_sf' should have at least two rows", call. = FALSE)
  }
  
  geom <- sf::st_geometry_type(points_sf) |> as.character() |> unique()
  
  if (length(geom) > 1) {
    stop("Argument 'points_sf' cannot contain different geometries", 
         call. = FALSE)
  }
  
  if (!("POINT" %in% geom)) {
    stop("Geometry of 'points_sf' must be of type POINT", call. = FALSE)
  }
  
  if (!("group" %in% colnames(points_sf))) {
    stop("The column 'group' is absent from 'points_sf'", call. = FALSE)
  }
  
  
  ## Check 'from' argument ----
  
  if (missing(from)) {
    stop("Argument 'from' is required", call. = FALSE)
  }
  
  if (!is.numeric(from) || length(from) != 1) {
    stop("Argument 'from' must be an integer of length 1", call. = FALSE)
  }
  
  if (!(from %in% 1:nrow(points_sf))) {
    stop("Argument 'from' must be between 1 and number of rows in 'points_sf'",
         call. = FALSE)
  }
  
  
  ## Check 'to' argument ----
  
  if (missing(to)) {
    stop("Argument 'to' is required", call. = FALSE)
  }
  
  if (!is.numeric(to) || length(to) != 1) {
    stop("Argument 'to' must be an integer of length 1", call. = FALSE)
  }
  
  if (!(to %in% 1:nrow(points_sf))) {
    stop("Argument 'to' must be between 1 and number of rows in 'points_sf'",
         call. = FALSE)
  }
  
  
  ## Extract segment ----
  
  points_sf[from:to, ] |> 
    dplyr::group_by(.data$group) |> 
    dplyr::summarise(do_union = FALSE) |> 
    sf::st_cast("LINESTRING") |> 
    dplyr::ungroup()
}



#' Check for boolean values
#' 
#' @param boolean a `logical` value of length 1.
#' 
#' @noRd

check_logical_value <- function(boolean) {
  
  if (is.null(boolean)) {
    stop("The argument '", deparse(substitute(boolean)), "' cannot be NULL", 
         call. = FALSE)
  }
  
  if (!is.logical(boolean) || length(boolean) != 1) {
    stop("The argument '", deparse(substitute(boolean)), "' must be a ", 
         "logical (TRUE or FALSE) of length 1", call. = FALSE)
  }
  
  if (is.na(boolean)) {
    stop("The argument '", deparse(substitute(boolean)), "' cannot be NA", 
         call. = FALSE)
  }
  
  invisible(NULL)
}



#' Check for the data.frame distance
#'
#' @param data a `data.frame` with the column `from`, `to` and `distance`. 
#'   The output of the function [distance_along()].
#'
#' @noRd

check_dist_data <- function(data) {
  
  if (missing(data)) {
    stop("Argument 'data' is required ", 
         "(output of the function distance_along())", call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data.frame ", 
         "(output of the function distance_along())", call. = FALSE)
  }
  
  if (!("from" %in% colnames(data))) {
    stop("The column 'from' is absent from the 'data' data.frame ", 
         "(output of the function distance_along())", call. = FALSE)
  }
  
  if (!("to" %in% colnames(data))) {
    stop("The column 'to' is absent from the 'data' data.frame ", 
         "(output of the function distance_along())", call. = FALSE)
  }
  
  if (!("distance" %in% colnames(data))) {
    stop("The column 'distance' is absent from the 'data' data.frame ", 
         "(output of the function distance_along())", call. = FALSE)
  }
  
  if (nrow(data) == 0) {
    stop("Argument 'data' must have at least one row", call. = FALSE)
  }
  
  if (!is.character(data$"from")) {
    stop("The column 'from' of the 'data' data.frame must be a character", 
         call. = FALSE)
  }
  
  if (!is.character(data$"to")) {
    stop("The column 'to' of the 'data' data.frame must be a character", 
         call. = FALSE)
  }
  
  invisible(NULL)  
}
