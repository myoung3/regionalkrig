
library(reshape2)

#=======================================#
#=======================================#

#' Parse Loading Types
#'
#' A function that 
#' @param x See rawdata under PLSK.full
#' @return String vector listing labels for groups of variables
#' @keywords 
#' @export
#' @examples 

parse_loadingtype <- function(x){ 
  new.x         <- x$model.obj$PLS$loadings[,1:2]
  out <- data.table(x = new.x, type = rep(as.character(NA), dim(new.x)[1]), subtype = rep(as.character(NA), dim(new.x)[1]))
  out[grep("^pop", rownames(new.x)),"type"] <- "population buffer"
  out[grep("^lu",  rownames(new.x)),"type"] <- "land use (historical) buffer"
  out[grep("^rlu", rownames(new.x)),"type"] <- "land use (2006) buffer"
  out[grep("^ndvi",rownames(new.x)),"type"] <- "ndvi buffer"
  out[grep("^imp", rownames(new.x)),"type"] <- "impervious surface buffer"
  out[grep("^elev",rownames(new.x)),"type"] <- "elevation"
  out[grep("^tl",  rownames(new.x)),"type"] <- "truck length buffer"
  out[grep("^intersect", rownames(new.x)),"type"] <- "intersection buffer"
  out[grep("^ll",  rownames(new.x)),"type"] <- "roadway length buffer"
  out[grep("^em",  rownames(new.x)),"type"] <- "emissions buffer"
  out[grepl("^log10", rownames(new.x)) &  ( grepl("a1", rownames(new.x)) | grepl("a2", rownames(new.x)) | 
                                            grepl("a3", rownames(new.x))) ,"type"] <- "roadway proximity"
  out[grepl("^log10", rownames(new.x)) & !( grepl("a1", rownames(new.x)) | grepl("a2", rownames(new.x)) | 
                                            grepl("a3", rownames(new.x))) ,"type"] <- "proximity other"
  out[is.na(out$type), "type"] <- "other" 

  out[, type:=as.factor(type)]
  levels(out$type) <- rev(c("roadway proximity", 
                            "roadway length buffer",
                            "intersection buffer", 
                            "truck length buffer",
                            "population buffer", 
                            "land use (1990s) buffer",
                            "land use (historical) buffer", 
                            "ndvi buffer",
                            "impervious surface buffer", 
                            "elevation", 
                            "emission buffer", 
                            "proximity other",
                            "other"))
  stopifnot(sum(is.na(out$type))==0)
  out$varname <- rownames(new.x)
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

plot_pls <- function(x){
  loadings_yrs  <- parse_loadingtype(x)
  names(loadings_yrs) <- c("Comp1","Comp2","type","subtype","varname")

  num.vars <- length(levels(loadings_yrs$type))
  pt.label <- levels(loadings_yrs$type)
  dummyx  <- rep(range(loadings_yrs$Comp1),2)
  dummyy  <- c(1, 1, num.vars, num.vars)

  par(mar=c(5,15,4,2))
  
  set.ticks <- nice_x_axis(loadings_yrs$Comp1)
  plot(dummyy~dummyx, axes = FALSE, frame.plot=TRUE, type="n", ylab = "")
  axis(2, las = 2, labels = pt.label, at=seq(num.vars,1,by=-1))
  axis(1, las = 2, labels = set.ticks, at=set.ticks)

  j<-1
  for (i in seq(num.vars,1,by=-1)){
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
  while((x.inc/scale.x) > 10 | (x.inc/scale.x) < 0.1 ){
    scale.x <- scale.x*scale.x
    if (scale.x < 1) n.dig <- n.dig + 1
  }
  new.min <- round(x.range[1]-x.inc, n.dig)
  seq(from = new.min, by = scale.x, length.out = 10)

}


