display.pointset <- function(file){
  X <- read.csv(file, header=FALSE)
  dev.new()
  index <- which( X[,1]=="Fixed" | X[,1]=="Moving" | X[,1]=="FixedTransformed" )
  plot(X[index, 2:3], col=X[index,1])
  title(file)
}

display.pointset("jensen-points.csv")
display.pointset("euclidean-points.csv")
display.pointset("trimmed-euclidean-points.csv")
display.pointset("weighted-euclidean-points.csv")

