#' First Hierarchy Function
#'
#' @param data A data frame with four columns:\cr
#' Initial Latitude | Initial Longitude | Final Latitude | Final Longitude
#' @param numK Initial number of clusters in the first call of k-means in the global hierarchy.
#' @param limitsSeparation Within cluster distance threshold to determine if a global cluster must be separated into two new clusters.
#' @param maxDist Meter distance threshold used to re-estimate centroids in global hierarchy.
#'
#' @return Returns an S3 class object similar to kmeans S3 Class, with eight properties.
#' @export
#'
#' @examples
#' data(ODMeansTaxiData)
#' first_hierarchy_data = first_hierarchy(ODMeansTaxiData, 10, 300, 1000)
first_hierarchy <- function(data, numK, limitsSeparation, maxDist) {

  #ejemplo sintetico, synthetic_data

  #Check if data is valid#
  #Latitude must be a value between -90 and 90
  if (!(all(data[,1]>=-90) && all(data[,1]<=90)) || !(all(data[,3]>=-90) && all(data[,3]<=90))){stop("Invalid Latitude. Must be a value between -90 and 90")}
  #Longitude must be a value between -180 and 180
  if (!(all(data[,2]>=-180) && all(data[,2]<=180)) || !(all(data[,4]>=-180) && all(data[,4]<=180))){stop("Invalid Longitude. Must be a value between -180 and 180")}

  #Check if parameters are numbers
  if (!(is.numeric(numK))){stop("numK must be a number")}
  if (!(is.numeric(limitsSeparation))){stop("limitsSeparation must be a number")}
  if (!(is.numeric(maxDist))){stop("maxDist must be a number")}

  #Iterative process to determine optimal number of clusters based on WithinClusterDistance.
  odDataframe=data[,c(1:4)]
  checkClusters=T
  while(checkClusters){
    checkClusters=F

    #Creation of clustered data based on centers = numK
    clusteredData=stats::kmeans(odDataframe, centers = numK, iter.max= 200, algorithm="Lloyd")
    newCenters=clusteredData$centers
    numCenters=nrow(newCenters)

    #Analyzing the behavior of each cluster
    for (i in c(1:numCenters)){
      if (nrow(odDataframe[clusteredData$cluster==i,])>1){ #At least two rows to analyze.
        #The center is separated into two. If the change is significant clusteredData is overwritten with the new clusters.
        subCluster=stats::kmeans(odDataframe[clusteredData$cluster==i,], centers = 2, iter.max= 200, algorithm="Lloyd")
        if (subCluster$tot.withinss<clusteredData$withinss[i]-limitsSeparation){
          checkClusters=T
          newCenters[i,]=subCluster$centers[1,]
          newCenters=rbind(newCenters,subCluster$centers[2,])
        }
      }
    }
    numK=newCenters #numK now is a data frame
  }

  #Joining process
  for (indexCol in c(1:ncol(newCenters))){

    #For each axis calculates the distance between points. If there are distances shorter than maxDist, start the joining process.
    if ((indexCol==1) | (indexCol==3)){
      tempLonLat=cbind(0, clusteredData$centers[,indexCol])
    } else {
      tempLonLat=cbind(clusteredData$centers[,indexCol],0)
    }

    #Calculation of distance matrix
    distanceMatrix=geosphere::distm(tempLonLat, fun = geosphere::distHaversine)
    distanceMatrix[lower.tri(distanceMatrix,T)]=maxDist
    #Obtaining close indexes
    indexToJoin=as.data.frame(which(distanceMatrix<maxDist,T))
    JointFlag = any(nrow(indexToJoin)>0)

    #Joining process
    while (JointFlag){
      #Vector pointing out what indexes have to be join
      tempIndex=rep(T,nrow(indexToJoin))
      #The first two element to join corresponds to the element of the first row
      subsetJoin=as.double(indexToJoin[1,])
      tempIndex[1]=F

      flagIndex=T #First loop
      while(flagIndex){
        flagIndex=F
        #Joining process of the subset
        for (j in c(1:nrow(indexToJoin))){
          ##Just in case three rows are connected through the last row, for example [1, 14][22, 10][14, 10]
          if (tempIndex[j] & ((sum(indexToJoin[j,1]==subsetJoin)+sum(indexToJoin[j,2]==subsetJoin))>0)){
            subsetJoin=union(subsetJoin,as.double(indexToJoin[j,]))
            tempIndex[j]=F
            flagIndex=T
          }
        }
      }
      #The new center consists in the mean
      clusteredData$centers[subsetJoin,indexCol]=mean(clusteredData$centers[subsetJoin,indexCol])
      #It is updated indexToJoin (erasing the indexes that already been joined)
      indexToJoin=as.data.frame(indexToJoin[tempIndex,])

      #If there is no index to join finish the loop
      if (sum(tempIndex)==0){
        JointFlag=F
      }
    }
  }
  #Calculation of vector of within-cluster sum of squares, one component per cluster.
  withinss = numeric()
  for (center in c(1:numCenters)){
    withinss = c(withinss, sum(rowSums(sweep(odDataframe[clusteredData$cluster==center,1:4], 2, clusteredData$centers[center,], `-`)^2)))
  }
  #The total sum of squares.
  totss = sum(scale(as.matrix(odDataframe), scale = FALSE)^2)
  #Total within-cluster sum of squares
  tot.withinss = sum(withinss)

  #Structure to return
  finalCluster=structure(list(cluster=clusteredData$cluster,
                              centers=clusteredData$centers,
                              totss = totss,
                              withinss = withinss,
                              tot.withinss = tot.withinss,
                              betweenss=(totss-tot.withinss),
                              size=clusteredData$size,
                              level_hierarchy = rep("Global",nrow(clusteredData$centers))))
  class(finalCluster) <- "ODMeans"
  return(finalCluster)
}
