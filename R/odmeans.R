#' ODMeans Function
#'
#' @param data A data frame with four columns:\cr
#' Initial Latitude | Initial Longitude | Final Latitude | Final Longitude
#' @param numKGlobal Initial number of clusters in the first call of k-means in the global hierarchy.
#' @param limitSeparationGlobal Within cluster distance threshold to determine if a global cluster must be separated into two new clusters.
#' @param maxDistGlobal Meter distance threshold used to re-estimate centroids in global hierarchy.
#' @param distHierarchical Meter distance threshold between origin and destination to generate new local clusters from a first layer cluster
#' @param numKLocal Initial number of clusters in the first call of k-means in the local hierarchy.
#' @param limitSeparationLocal Within cluster distance threshold to determine if a local cluster must be separated into two new clusters.
#' @param maxDistLocal Meter distance threshold used to re-estimate centroids in local hierarchy.
#' @param kmeans_pp Boolean value, if TRUE it initialize centroids using kmeans++.
#'
#' @return Returns an S3 class object similar to kmeans S3 Class, with eight properties.
#' @export
#'
#' @examples
#' data(ODMeansTaxiData)
#' odmeans_data = odmeans(ODMeansTaxiData, 10, 300, 1000, 2200, 3, 50, 100)
odmeans <- function(data, numKGlobal, limitSeparationGlobal, maxDistGlobal, distHierarchical, numKLocal, limitSeparationLocal, maxDistLocal, kmeans_pp=FALSE){

  global_cluster = first_hierarchy(data, numKGlobal, limitSeparationGlobal, maxDistGlobal, kmeans_pp)
  local_cluster = second_hierarchy(data, global_cluster, distHierarchical, numKLocal, limitSeparationLocal, maxDistLocal)
  return(local_cluster)
}


