#' ODMeans
#'
#' @param data A data frame with four columns:
#' Initial Latitude | Initial Longitude | Final Latitude | Final Longitude
#' @param numK Initial number of clusters in the first call of K-Means.
#' @param limitsSeparation Range to determine if a drastic change has happened between a cluster and its separation.
#' A bigger value makes more difficult to separate a cluster.
#' @param maxDist Maximum distance to join two points. This is based on the euclidean distance.
#' @param distHierarchical Maximum distance to create a new hierarchy per cluster
#'
#' @return Returns a structure that contains the final centers, clusters, sizes and hierarchy
#' @export
#'
#' @examples
#' data(ODMeansSampleData)
#' od_means(ODMeansSampleData, 5, 200, 2500, 500)
od_means <- function(data, numK, limitsSeparation, maxDist, distHierarchical){
  kcluster = dinamic_clusters(data, numK, limitsSeparation, maxDist)
  return(hierarchical_clusters(data, kcluster, distHierarchical))
}
