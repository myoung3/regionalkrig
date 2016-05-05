
poly.int2 <- function(xpoly,xtest) {
  # xtest is the mx2 array of points to be tested.
  # xpoly is an nx2 array of points defining the polygon,
  # assuming the nth point is a replicate of the first.
  m <- nrow(xtest)
  n <- nrow(xpoly)
  xp <- xpoly[,1]; yp <- xpoly[,2]
  c <- rep(F,m)
  for (k in 1:m) {
   x <- xtest[k,1]; y <- xtest[k,2]
   crossings <- 0
  for (i in 1:(n-1)) {
   if ( ( xp[i]<x & x<xp[i+1] ) | ( xp[i]>x & x>xp[i+1] )) {
    t <- (x-xp[i+1])/(xp[i]-xp[i+1])
    cy <- t*yp[i] + (1-t)*yp[i+1]
    #if (y == cy) #boundary#
    if (y>cy) crossings <- crossings+1
   }
  }
  if ((crossings/2 - round(crossings/2))!=0) c[k] <- T
 }
 return(c)
}

