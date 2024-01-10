#' Graph ODMeans Function
#'
#' @param odmeans_data It receives an object from S3 ODMeans class. However, it can also work with objects from similar classes like S3 k-Means
#' @param title It receives an string, and corresponds to the title of the plot.
#' @param maptype It receives a string with the type of the map. Depending on the map selected, it will change the background of it. The possible values are: “terrain”, “satellite”, “roadmap”, “hybrid”.
#' @param zoom An integer from 3 (continent) to 21 (building), which controls the level of zoom applied to the map.
#' @param add_cluster Receives TRUE or FALSE value. When True, it adds the number of the cluster to the arrows.
#'
#' @return A ggplot graph showing a map with the centers of the clusters.
#' @export
#'
#' @examples
#' data(ODMeansTaxiData)
#' odmeans_data = odmeans(ODMeansTaxiData, 10, 300, 1000, 2200, 3, 50, 100)
#' odmeans_plot = odmeans_graph(odmeans_data, "ODMeans Taxi Graph", "roadmap", 11, FALSE)
#'
odmeans_graph <- function(odmeans_data, title="ODMeans Graph", maptype="roadmap", zoom=4, add_cluster=TRUE){

  # Obtain all information to graph
  dfODMeans = data.frame(odmeans_data[["centers"]])
  colnames(dfODMeans) <- c("OriginLatitude", "OriginLongitude", "DestinationLatitude", "DestinationLongitude")

  clusters = c(1:nrow(dfODMeans))

  if (is.null(odmeans_data[["level_hierarchy"]]) != TRUE)
    Hierarchy = data.frame("Hierarchy" = odmeans_data[["level_hierarchy"]])
  else
    Hierarchy = rep("no_hierarchy", nrow(dfODMeans))

  dfODMeans = cbind(dfODMeans, clusters, Hierarchy)


  # Coordinates of the map
  lowerleftlon = min(dfODMeans["OriginLongitude"], dfODMeans["DestinationLongitude"])
  upperrightlon = max(dfODMeans["OriginLongitude"], dfODMeans["DestinationLongitude"])
  lowerleftlat = min(dfODMeans["OriginLatitude"], dfODMeans["DestinationLatitude"])
  upperrightlat = max(dfODMeans["OriginLatitude"], dfODMeans["DestinationLatitude"])

  center =  c(lon = mean(c(lowerleftlon, upperrightlon)), lat = mean(c(lowerleftlat, upperrightlat)))


  coordinates = c(lowerleftlon, lowerleftlat, upperrightlon, upperrightlat )

  map <- tryCatch(
    {
      ggmap::get_googlemap(center = center, zoom = zoom, maptype = maptype)
    },
    error = function(e) {
      message("An error occurred while getting the map, try changing zoom, maptype or data.")
      message(e)
      return(NULL)  # Return NULL or an appropriate value if there's an error
    }
  )

  # Check if map is NULL and return a NULL map
  if (is.null(map)) {
    return(NULL)
  }

  ## Add cluster include number of clusters
  if (add_cluster == TRUE)
    final_map = ggmap::ggmap(map)+
    ggplot2::geom_segment(data=dfODMeans,
                 ggplot2::aes_string(y= "OriginLatitude",
                     yend= "DestinationLatitude",
                     x= "OriginLongitude",
                     xend= "DestinationLongitude",
                     color = "Hierarchy"),
                 linewidth= 1,
                 arrow= ggplot2::arrow(length = ggplot2::unit(0.4, "cm"))) +
    ggplot2::scale_fill_manual(aesthetics = "colour", values = c("Local" = "royalblue","Global" = "red")) +
    ggrepel::geom_label_repel(data=dfODMeans,
                     ggplot2::aes_string(label="clusters", y="OriginLatitude", x="OriginLongitude"),
                     alpha= 0.7,
                     max.overlaps = 26,
                     na.rm=TRUE,
                     xlim = c(lowerleftlon, upperrightlon),
                     ylim = c (lowerleftlat, upperrightlat)) +
    ggplot2::theme_minimal() +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude") +
    ggplot2::labs(title = paste(title))

  else {
    final_map = ggmap::ggmap(map)+
      ggplot2::geom_segment(data=dfODMeans,
                  ggplot2::aes_string(y= "OriginLatitude",
                       yend= "DestinationLatitude",
                       x= "OriginLongitude",
                       xend= "DestinationLongitude",
                       color = "Hierarchy"),
                   linewidth= 1,
                   arrow= ggplot2::arrow(length = ggplot2::unit(0.4, "cm"))) +
      ggplot2::scale_fill_manual(aesthetics = "colour", values = c("Local" = "royalblue","Global" = "red")) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("Longitude") +
      ggplot2::ylab("Latitude") +
      ggplot2::labs(title = paste(title))
  }

  return(final_map)
}
