# Function obtained from http://www.phaget4.org/R/image_matrix.html

myImagePlot <- function(x, border=TRUE, min=NULL, max=NULL, 
                        cex=0.7, lwd=0.3, sub=FALSE, ...){
  superl = c("African", "American", "European", "East Asian", "South Asian")
  supers = c(0, 7, 11, 16, 21, 26)
  if (is.null(min)) min <- min(x)
  if (is.null(max)) max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(1,0,length=256),  # Red
                    seq(1,0,length=256),  # Green
                    seq(1,1,length=256))  # Blue
  ColorRamp <- c(ColorRamp, rgb( seq(0,1,length=256),  # Red
                                 seq(0,0,length=256),  # Green
                                 seq(1,0,length=256)))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(6,6,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  if (is.null(cex)){cex = 1 - (length(yLabels) / 200)}
  if (cex < 0.3) cex = 0.3
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, las=2, cex.axis=cex)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las=2,
       cex.axis=cex)

  if (border == TRUE){
    for (i in 0:length(xLabels)){
      yLen = length(yLabels)
      if (sub){
        if (i %in% supers){
          axis(2, pos = i+0.5, labels = FALSE, at = 0.5:(yLen+0.5),
               lwd.ticks = 0, lwd = 1.5)
        }
      }
      axis(2, pos = i+0.5, labels = FALSE, at = 0.5:(yLen+0.5),
           lwd.ticks = 0, lwd = lwd)
    }
    
    for (i in 0:length(yLabels)){
      xLen = length(xLabels)
      if (sub){
        if (i %in% (26-supers)){
          axis(1, pos = i+0.5, labels = FALSE, at = 0.5:(xLen+0.5), 
               lwd.ticks = 0, lwd = 1.5)
        }
      }
      axis(1, pos = i+0.5, labels = FALSE, at = 0.5:(xLen+0.5), 
           lwd.ticks = 0, lwd = lwd)
    }
  }
  # Color Scale
  par(mar = c(6,2.5,2.5,1))
  image(1, y=ColorLevels,
        z=matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n", las=1)
  layout(1)
}
