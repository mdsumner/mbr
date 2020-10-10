#' Minimum bounding rectangle
#'
#' Find the minimum-area-rectangle for a given set of points.
#'
#' Input is an array of point coordinates. Output (the value of mbr) is an array
#' of the vertices of the minimum bounding rectangle (with the first one
#' repeated to close it). Note the complete absence of any trigonometric
#' calculations.
#'
#' Timing is limited by the speed of the convex hull algorithm, because the
#' number of vertices in the hull is almost always much less than the total.
#' Most convex hull algorithms are asymptotically O(n*log(n)) for n points: you
#' can compute almost as fast as you can read the coordinates.
#'
#' `MBR()` is an alias of `mbr()`.
#'
#' For an sf wrapper see https://github.com/mdsumner/mbr/issues/1
#' @author Bill Huber
#' @references
#' https://gis.stackexchange.com/questions/22895/finding-minimum-area-rectangle-for-given-points/22934?stw=2#22934
#' @param p input points (matrix of 2 columns)
#'
#' @return matrix of rectangle coordinates (5 points)
#' @export
#' @name mbr
#' @importFrom grDevices chull
#' @aliases MBR
#' @examples
#' xy <- quakes[, 1:2]
#' plot(mb0 <- mbr(xy)); lines(mb0)
#' points(xy)
mbr <- function(p) {
  p <- as.matrix(p)[, 1:2, drop = FALSE]
  # Analyze the convex hull edges
  a <- chull(p)                                   # Indexes of extremal points
  a <- c(a, a[1])                                 # Close the loop
  e <- p[a[-1],] - p[a[-length(a)], ]             # Edge directions
  norms <- sqrt(rowSums(e^2))                     # Edge lengths
  v <- e / norms                                  # Unit edge directions
  w <- cbind(-v[,2], v[,1])                       # Normal directions to the edges

  # Find the MBR
  vertices <- p[a, ]                              # Convex hull vertices
  x <- apply(vertices %*% t(v), 2, range)         # Extremes along edges
  y <- apply(vertices %*% t(w), 2, range)         # Extremes normal to edges
  areas <- (y[1,]-y[2,])*(x[1,]-x[2,])            # Areas
  k <- which.min(areas)                           # Index of the best edge (smallest area)

  # Form a rectangle from the extremes of the best edge
  cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])
}

#' @name mbr
#' @export
MBR <- mbr
