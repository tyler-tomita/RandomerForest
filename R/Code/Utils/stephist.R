library(ggplot2)
library(proto)

geom_stephist <- function(mapping = NULL, data = NULL, stat="bin", position="identity", ...) {
  GeomStepHist$new(mapping=mapping, data=data, stat=stat, position=position, ...)
}

GeomStepHist <- proto(ggplot2:::Geom, {
  objname <- "stephist"
  
  default_stat <- function(.) StatBin
  default_aes <- function(.) aes(colour="black", size=0.5, linetype=1, alpha = NA)
  
  reparameterise <- function(., df, params) {
    transform(df,
              ymin = pmin(y, 0), ymax = pmax(y, 0),
              xmin = x - width / 2, xmax = x + width / 2, width = NULL
    )
  }
  
  draw <- function(., data, scales, coordinates, ...) {
    data <- as.data.frame(data)[order(data$x), ]
    
    n <- nrow(data)
    i <- rep(1:n, each=2)
    newdata <- rbind(
      transform(data[1, ], x=xmin, y=0),
      transform(data[i, ], x=c(rbind(data$xmin, data$xmax))),
      transform(data[n, ], x=xmax, y=0)
    )
    rownames(newdata) <- NULL
    
    GeomPath$draw(newdata, scales, coordinates, ...)
  }
  guide_geom <- function(.) "path"
})