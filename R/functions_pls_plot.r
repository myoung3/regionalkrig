
library(reshape2)

#=======================================#
#=======================================#

#' Parse Loading Types
#'
#' A function that labels the loadings matrix with variable types
#' @param x See rawdata under PLSK.full
#' @return data.table with loading components labeled
#' @keywords 
#' @export
#' @examples 

parse_loadingtype <- function(x){ 
  new.x         <- x$model.obj$PLS$loadings[,1:2]
  out <- data.table(x = new.x, type = rep(as.character(NA), dim(new.x)[1]), subtype = rep(as.character(NA), dim(new.x)[1]))
  out$varname <- rownames(new.x)
  level.order <- rev(c("roadway proximity", 
                       "roadway length buffer",
                       "intersection buffer", 
                       "truck length buffer",
                       "population buffer", 
                       "land use (2006) buffer",
                       "land use (historical) buffer", 
                       "ndvi buffer",
                       "impervious surface buffer", 
                       "elevation", 
                       "emissions buffer", 
                       "proximity other",
                       "satellite NO2",
                       "other"))
  out[grep("^pop", out$varname),]$type <- "population buffer"
  out[grep("^lu",  out$varname),]$type <- "land use (historical) buffer"
  out[grep("^rlu", out$varname),]$type <- "land use (2006) buffer"
  out[grep("^ndvi",out$varname),]$type <- "ndvi buffer"
  out[grep("^imp", out$varname),]$type <- "impervious surface buffer"
  out[grep("^elev",out$varname),]$type <- "elevation"
  out[grep("^tl",  out$varname),]$type <- "truck length buffer"
  out[grep("^intersect", out$varname),]$type <- "intersection buffer"
  out[grep("^ll",  out$varname),]$type <- "roadway length buffer"
  out[grep("^em",  out$varname),]$type <- "emissions buffer"
  out[grep("^no2",  out$varname),]$type <- "satellite NO2"
  out[grepl("^log10", out$varname) &  ( grepl("a1", out$varname) | grepl("a2", out$varname) | 
                                            grepl("a3", out$varname)) ,]$type <- "roadway proximity"
  out[grepl("^log10", out$varname) & !( grepl("a1", out$varname) | grepl("a2", out$varname) | 
                                            grepl("a3", out$varname)) ,]$type <- "proximity other"
  out[is.na(out$type),]$type <- "other" 
  out$type <- factor(out$type, levels = level.order)
  stopifnot(sum(is.na(out$type))==0)
  setkey(out, "varname")
  out
}

#=======================================#
#=======================================#

#' Plot PLS: Plots PLS loadings
#'
#' A function that 
#' @param x See rawdata under PLSK.full
#' @return NA, produces plot
#' @keywords 
#' @export
#' @examples 

plot_pls <- function(x, ...){
  loadings_yrs  <- parse_loadingtype(x)
  names(loadings_yrs) <- c("Comp1","Comp2","type","subtype","varname")

  num.vars <- length(levels(loadings_yrs$type))
  pt.label <- levels(loadings_yrs$type)
  dummyx  <- rep(range(loadings_yrs$Comp1),2)
  dummyy  <- c(1, 1, num.vars, num.vars)

  dev.new(height = 6, width = 7.5)
  par(mar=c(5,15,4,2))
  
  set.ticks <- nice_x_axis(loadings_yrs$Comp1)
  plot(dummyy~dummyx, axes = FALSE, frame.plot=TRUE, type="n", ylab = "", ...)
  axis(2, las = 2, labels = pt.label, at=c(1:num.vars))
  axis(1, las = 2, labels = set.ticks, at=set.ticks)

  j<-1
  for (i in c(1:num.vars)){
    current <- loadings_yrs$Comp1[loadings_yrs$type == levels(loadings_yrs$type)[i]]
    points(current,rep(i, length(current)) )
    j <-  j+1
  }
  abline(v=0,col="red")
}

#=======================================#
#=======================================#

#' Nice X Axis
#'
#' A function that 
#' @param x numeric vector
#' @return numeric vector, for x-axis locations
#' @keywords 
#' @export
#' @examples 

nice_x_axis <- function(x){
  max.x   <- max(abs(x), na.rm=TRUE)
  x.range <- range(x)
  x.inc   <- (x.range[2] - x.range[1])/11 

  if (x.inc  > 1){ scale.x <-  10 }
  if (x.inc  < 1){ scale.x <- 0.1 }
  if (x.inc == 1){ scale.x <-   1 }
  n.dig <- 0
  while((x.inc/scale.x) > 10 | (x.inc/scale.x) < 1 ){
    scale.x <- scale.x*scale.x
    if (scale.x < 1) n.dig <- n.dig + 1
  }
  n.dig <- n.dig + 1
  new.min <- round(x.range[1]-x.inc, n.dig)
  round(seq(from = new.min, to = max.x, by = scale.x), n.dig)
}


