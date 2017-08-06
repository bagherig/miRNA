# getDendrogram.R

#' Draw a dendrogram from a square matrix.
#'
#' @param squareMatrix The square matrix to be plotted as a dendrogram.
#'
#' @return A dendrogram object.
#'
#' @examples getDendrogram(pFstMatrix.alt.roa)
#' 
getDendrogram <- function(squareMatrix){
  # Cluster.
  clustered = hclust(as.dist(squareMatrix))
  groupColors = c("blue", "red", "purple", "darkgreen", "brown")
  clusterGroups = cutree(clustered, 6)
  
  # Function to color leafs based on the groups.
  colorGroups <- function(n) {
    if (is.leaf(n)) {
      a = attributes(n)
      labCol <- groupColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    return(n)
  }
  # Make a dendrogram
  dendro = as.dendrogram(clustered)
  # Color different groups.
  clusDendro = dendrapply(hcd, colorGroups)
  
  return(clusDendro)
}
