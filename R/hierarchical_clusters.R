#' Hierarchical Clusters
#'
#' @param data A data frame with four columns:\cr
#' Initial Latitude | Initial Longitude | Final Latitude | Final Longitude
#' @param Kcluster An ODMeans structure, result of function dinamic_clusters.
#' @param distHierarchical Maximum distance to create a new hierarchy per cluster.
#'
#' @return Hierarchical Clusters returns an object similar of class "kmeans". It is a list with at least the following components:
#'
#'cluster A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
#'centers A matrix of cluster centres.
#'totss The total sum of squares.
#'withinss Vector of within-cluster sum of squares, one component per cluster.
#'tot.withinss Total within-cluster sum of squares, i.e. sum(withinss).
#'betweenss The between-cluster sum of squares, i.e. totss-tot.withinss.
#'size The number of points in each cluster.
#'level_hierarchy Corresponds of the hierarchy level of the cluster, can be "Global" or "Local"
#' @export
#'
#' @examples
#' data(ODMeansSampleData)
#' hierarchical_clusters(ODMeansSampleData, dinamic_clusters(ODMeansSampleData, 5, 200, 2500), 500)
hierarchical_clusters <- function(data,Kcluster,distHierarchical) {
  #Structure to return
  finalCluster=Kcluster

  #Check if data is valid#
  #Latitude must be a value between -90 and 90
  if (!(all(data[,1]>=-90) && all(data[,1]<=90)) || !(all(data[,3]>=-90) && all(data[,3]<=90))){stop("Invalid Latitude. Must be a value between -90 and 90")}
  #Longitude must be a value between -180 and 180
  if (!(all(data[,2]>=-180) && all(data[,2]<=180)) || !(all(data[,4]>=-180) && all(data[,4]<=180))){stop("Invalid Longitude. Must be a value between -180 and 180")}
  odDataframe=data[,c(1:4)]

  #Check if parameters are numbers
  if (!(is.numeric(distHierarchical))){stop("distHierarchical must be a number")}

  #Analyzing the distances between origin and destiny per cluster
  distOD=diag(geosphere::distm(Kcluster$center[,c(2,1)],Kcluster$center[,c(4,3)]))

  for (i in c(1:length(distOD))){
    #If a distance is smaller than distHierarchical, a new hierarchy is created
    if (distOD[i]<distHierarchical){
      tempData=data[Kcluster$cluster==i,]
      #Creating the subcluster
      newSubCluster=dinamic_clusters(tempData,3,5,100)

      #Adding the subcluster to the new component finalCluster
      #It maps the previous points with the new points
      indexSC=as.vector(which(Kcluster$cluster==i))

      #In newSubCluster the indexes from 2 to K must be reemployed by length(KforAllCluster) to length(KforAllCluster)+k-1
      tempCluster=as.vector(newSubCluster$cluster)+max(finalCluster$cluster)-1
      #In newSubCluster the index 1 must be reemployed by i
      tempCluster[tempCluster==min(tempCluster)]=i
      #Changing the clusters
      finalCluster$cluster[indexSC]=tempCluster

      #Adding the new center to finalCluster
      finalCluster$centers[i,]=newSubCluster$centers[1,]
      finalCluster$centers=rbind(finalCluster$centers, newSubCluster$centers[2:nrow(newSubCluster$centers),])

      #Adding the new size to finalCluster
      finalCluster$size[i]=newSubCluster$size[1]
      finalCluster$size=c(finalCluster$size, newSubCluster$size[2:nrow(newSubCluster$centers)])

      #Adding the hierarchy to finalCluster
      finalCluster$level_hierarchy[i]="Local"
      finalCluster$level_hierarchy=c(finalCluster$level_hierarchy, rep("Local",nrow(newSubCluster$centers[2:nrow(newSubCluster$centers),])))
    }
  }
  ###Overwriting measures###
  #Calculation of vector of within-cluster sum of squares, one component per cluster.
  withinss = numeric()
  for (center in c(1:nrow(finalCluster$centers))){
    withinss = c(withinss, sum(rowSums(sweep(odDataframe[finalCluster$cluster==center,1:4], 2, finalCluster$centers[center,], `-`)^2)))
  }
  finalCluster$withinss = withinss
  #The total sum of squares.
  totss = sum(scale(as.matrix(odDataframe), scale = FALSE)^2)
  finalCluster$totss= totss
  #Total within-cluster sum of squares
  tot.withinss = sum(withinss)
  finalCluster$tot.withinss = tot.withinss
  #Betweeness
  finalCluster$betweenss=(totss-tot.withinss)

  return(finalCluster)
}
