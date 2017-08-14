# getDendrogram.R

#' Draw a dendrogram from a square matrix.
#'
#' @param squareMatrix The square matrix to be plotted as a dendrogram.
#'
#' @return A dendrogram object.
#'
#' @examples getDendrogram(pFstMatrix.alt.roa)
#' 
getDendrogram <- function(squareMatrix, numGroups){
  # Cluster.
  clustered = hclust(as.dist(squareMatrix))
  groupColors = c("blue", "red", "purple", "darkgreen", "brown")
  clusterGroups = cutree(clustered, numGroups)
  
  # Function to color leafs based on the groups.
  colorGroups <- function(n) {
    if (is.leaf(n)) {
      a = attributes(n)
      labCol <- groupColors[clusterGroups[which(names(clusterGroups) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    return(n)
  }
  # Make a dendrogram
  dendro = as.dendrogram(clustered)
  # Color different groups.
  clusDendro = dendrapply(dendro, colorGroups)
  
  return(clusDendro)
}
