#' Compute the distance between two locations along a linear shape
#' 
#' @description
#' Computes distance between two sites along a linear shape (e.g. a river, a 
#' road, etc.). This function uses parallel computing to speed up the distance
#' computation.
#' 
#' @param sites an `sf` object of type `POINT`. A spatial object containing 
#'   coordinates of sites. Note that the first column must be the 
#'   site identifier.
#' 
#' @param along an `sf` object of type `LINESTRING`. A spatial object 
#'   containing coordinates of the linear shape (e.g. a river, a road, etc.)
#'   to follow while computing distances between the two sites.
#' 
#' @param density a `numeric` of length 1. The density of points to sample on 
#'   the linear structure. See [sf::st_line_sample()] for further detail. 
#'   Default is `1/100`.
#' 
#' @param type a `character` of length 1. The method to sample points on the 
#'   linear shape. Either `regular` (default) or `random`.
#'   
#' @param mc.cores an `integer` of length 1. The number of cores to use 
#'   (must be lesser than the number of cores available on the machine obtained
#'   with [parallel::detectCores()]).
#' 
#' @param ... other argument to pass to  [sf::st_line_sample()].
#' 
#' @return A three-column `data.frame` with:
#'   - `from`, the first site
#'   - `to`, the second site
#'   - `weight`, the distance between the two sites along the linear shape
#' 
#' @export
#' @import rlang
#'
#' @examples
#' ## Add an example ----

distance_along <- function(sites, along, density = 0.01, type = "regular", 
                           mc.cores = parallel::detectCores() - 1, ...) {

  ## Check 'sites' argument ----
  
  if (missing(sites)) {
    stop("Argument 'sites' (spatial layer of sites) is required", 
         call. = FALSE)
  }
  
  if (!inherits(sites, "sf")) {
    stop("The object 'sites' must be an 'sf' object", 
         call. = FALSE)
  }
  
  if (nrow(sites) < 2) {
    stop("Argument 'sites' should have at least two rows (sites)", 
         call. = FALSE)
  }
  
  if (ncol(sites) < 2) {
    stop("Argument 'sites' should have at least two columns: site label and ", 
         "geometry", call. = FALSE)
  }
  
  geom <- sf::st_geometry_type(sites) |> 
    as.character() |> 
    unique()
  
  if (length(geom) > 1) {
    stop("Argument 'sites' (spatial layer of sites) cannot contain different ", 
         "geometries", call. = FALSE)
  }
  
  if (!("POINT" %in% geom)) {
    stop("Sites geometry must be of type POINT", call. = FALSE)
  }
  
  if (any(duplicated(sites[ , 1, drop = TRUE]))) {
    stop("The argument 'sites' cannot contain duplicates", call. = FALSE)
  }
  
  if (is.na(sf::st_crs(sites))) {
    stop("The 'sites' layer has not a valid CRS", call. = FALSE)
  }
  
  
  ## Check 'along' argument ----
  
  if (missing(along)) {
    stop("Argument 'along' (spatial layer of linear shape) is required", 
         call. = FALSE)
  }
  
  if (!inherits(along, "sf")) {
    stop("The object 'along' must be an 'sf' object", 
         call. = FALSE)
  }
  
  if (nrow(along) != 1) {
    stop("Argument 'along' (linear shape) should have exactly one row", 
         call. = FALSE)
  }
  
  geom <- sf::st_geometry_type(along) |> 
    as.character() |> 
    unique()
  
  if (!("LINESTRING" %in% geom)) {
    stop("Linear shape geometry must be of type LINESTRING", call. = FALSE)
  }
  
  if (is.na(sf::st_crs(along))) {
    stop("The 'along' layer has not a valid CRS", call. = FALSE)
  }
  
  
  ## Check 'type' argument ----
  
  type <- tolower(type)
  
  if (!(type %in% c("regular", "random"))) {
    stop("Argument 'type' must either 'regular' or 'random'", call. = FALSE)
  }
  
  
  ## Check density argument ----
  
  if (!is.numeric(density) || length(density) != 1) {
    stop("Argument 'density' must be a numeric of length 1", call. = FALSE)
  }
  
  if (density <= 0) {
    stop("Argument 'density' must be > 0", call. = FALSE)
  }
  
  
  ## Check for different CRS ----
  
  if (sf::st_crs(sites) != sf::st_crs(along)) {
    stop("Layers 'sites' and 'along' have different CRS", call. = FALSE)
  }
  
  
  ## Convert LINESTRING to POINTS ----
  
  message("- Sampling points on the LINESTRING...")
  
  sampled_points <- line_to_points(along, density, type, ...)
  
  
  ## Find nearest points on line to each site ----
  
  message("- Finding nearest points on LINESTRING for each point...")
  
  nearest_points <- unlist(parallel::mclapply(1:nrow(sites), function(i) {
    which.min(sf::st_distance(sites[i, ], sampled_points))
  }, mc.cores = mc.cores))
  
  
  ## Create correspondence table ----
  
  message("- Creating point by point matrix...")
  
  nearest_points <- data.frame("id"   = nearest_points, 
                               "site" = sites[ , 1, drop = TRUE])
  
  
  ## Get all combinations with two sites ----
  
  pairs_of_sites <- expand.grid("from" = nearest_points$"site", 
                                "to"   = nearest_points$"site", 
                                stringsAsFactors = FALSE)

  
  # Select the upper triangle ----
  
  pairs_of_sites$"distance" <- 1
  
  pairs_of_sites <- df_to_matrix(pairs_of_sites, lower = FALSE, diag = FALSE)
  pairs_of_sites <- matrix_to_df(pairs_of_sites)[ , 1:2]
  
  
  ## Create spatial segment between all sites pairs ----
  
  message("- Creating segment between pairs of points...")
  
  along_segments <- do.call(rbind.data.frame, 
                            parallel::mclapply(1:nrow(pairs_of_sites), 
                                               function(i) {
                                                 
    site_from <- which(nearest_points$"site" == pairs_of_sites[i, "from"])
    site_to   <- which(nearest_points$"site" == pairs_of_sites[i, "to"])
    
    points_to_line(points_sf = sampled_points,
                   from      = nearest_points[site_from, "id"],
                   to        = nearest_points[site_to, "id"]) |> 
      
    dplyr::mutate(from = pairs_of_sites[i, "from"],
                  to   = pairs_of_sites[i, "to"]) |> 
      
    dplyr::select(.data$from, .data$to)
    
  }, mc.cores = mc.cores))

  
  ## Compute distance between all sites pairs ----
  
  message("- Computing distance (segment length)...")
  
  distance_along <- sf::st_length(along_segments) |> as.numeric()
  
  
  ## Export final table ----
  
  message("- Cleaning output...")
  
  along_segments <- along_segments |> 
    dplyr::mutate(distance = distance_along) |> 
    dplyr::select(1, 2, 4) |> 
    sf::st_drop_geometry() |> 
    as.data.frame()
  
  # Add lower triangle ----
  
  to_add <- along_segments
  colnames(to_add)[1:2] <- c("to", "from")
  to_add <- to_add[ , c("from", "to", "distance")]
  
  along_segments <- rbind(along_segments, to_add)
  
  
  # Add diagonal ----
  
  sites <- unique(along_segments$"from")
  to_add <- data.frame("from" = sites, "to" = sites, "distance" = 0)
  
  along_segments <- rbind(along_segments, to_add)
  
  
  # Order table ----
  
  along_segments <- along_segments[with(along_segments, order(from, to)), ]
  rownames(along_segments) <- NULL
  
  along_segments
}
