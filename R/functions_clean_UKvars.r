

#unfinished--consider the functions that USE this object (ie, the process of prediction)
##don't forget it's model.obj$PLS not model.obj$pls now
make_predict_object <- function(rawdata, model.obj, pls.comps)
{
  ppts <- c()
  ppts$covars.pls <- scale(rawdata[, colnames(model.obj$covars.pls),drop=FALSE],
                       center=attr(model.obj$covars.pls, 'scaled:center'),
                       scale=attr(model.obj$covars.pls, 'scaled:scale'))
  ppts$PLS$scores <- predict(model.obj$PLS, newdata=ppts, comps=1:pls.comps, type='scores')  
  ppts$UK.covars <- as.matrix(rawdata[, colnames(model.obj$UK.covars),drop=FALSE])
  ppts$monitors <- rawdata$location_id
  ppts$obs <- length(ppts$monitors)
  ppts$coords <- rawdata[, c('lambert_x','lambert_y')]/1000
  ppts$gps <- rawdata[, c('longitude','latitude')]
  ppts$weco.monitors <- ppts$monitors[poly.int2(weco.polygon, ppts$gps)]
  ppts$east.monitors <- ppts$monitors[poly.int2(east.polygon, ppts$gps)]
  if (length(c(ppts$weco.monitors, ppts$east.monitors)) == 0)
    ppts$west.monitors <- ppts$monitors
  else
    if (length(ppts$weco.monitors) == 0)
      ppts$west.monitors <- ppts$monitors[-c(match(ppts$east.monitors, ppts$monitors))]
    else if (length(ppts$east.monitors) == 0)
      ppts$west.monitors <- ppts$monitors[-c(match(ppts$weco.monitors, ppts$monitors))]
    else
      ppts$west.monitors <- ppts$monitors[-match(unlist(list(ppts$east.monitors, ppts$weco.monitors)), ppts$monitors)]
  rownames(ppts$coords) <- ppts$monitors
  rownames(ppts$gps) <- ppts$monitors
  return(ppts)
}







#UK.varnames is a character vector corresponding to the names of variables to be included directly in the UK model (in addition to PLS components)
##new components of model.obj:
##covars.pls are the (scaled) covariates that go into the PLS models--covars except exclude and UK.vars\
##PLS is the pls output, comparable to "pls" in the original version of this funcion (renamed to spit out an error any time it comes across code that hasn't been changed)
##therefore the variables to be used in UK are cbind(1, model$PLS$scores[, 1:ncomp], model.obj$UK.covars)
##where ncomp is the number of pls.components desired

make_modeling_object <- function(rawdata, exclude.vars,UK.varnames, pls.comps=5)
{
  if("data.table" %in% class(rawdata)) rawdata <- as.data.frame(rawdata)

  covars.pls.names <- setdiff(setdiff(colnames(rawdata), exclude.vars), UK.varnames)
  model.obj <- c()
  model.obj$y <- sqrt(rawdata[, 'pollutant_conc'])
  model.obj$UK.covars <- as.matrix(rawdata[, UK.varnames, drop=FALSE])
  model.obj$covars.pls <- scale(as.matrix(rawdata[,covars.pls.names,drop=FALSE]))
  model.obj$PLS <- pls::plsr(y ~ covars.pls, ncomp=pls.comps, validation='none', data=model.obj)
  model.obj$obs <- nrow(rawdata)
  model.obj$coords <- rawdata[, c('lambert_x', 'lambert_y')]/1000
  model.obj$gps <- rawdata[, c('longitude', 'latitude')]
  model.obj$monitors <- rawdata$native_id
  model.obj$weco.monitors <- model.obj$monitors[poly.int2(weco.polygon, model.obj$gps)]
  model.obj$east.monitors <- model.obj$monitors[poly.int2(east.polygon, model.obj$gps)]
  model.obj$west.monitors <- model.obj$monitors[-match(unlist(list( model.obj$east.monitors,  model.obj$weco.monitors)), model.obj$monitors)]
  rownames(model.obj$coords) <- model.obj$monitors
  rownames(model.obj$gps) <- model.obj$monitors
  names(model.obj$y) <- model.obj$monitors
  return (model.obj)
}



#============================================#
# transforms all distance variables to log10 #
# and lower bounds them to 10 meters         #
#============================================#
log_transform_distances <- function(my.data)
{
  distance.vars <- grep("^m_to_", colnames(my.data))
  new.varnames <- c()
  for (i in distance.vars)
  {
    newcol.index <- 1 + length(colnames(my.data))
    # throw out records with covariate values of -1
    my.data <- my.data[!my.data[, i] == -1, ]
    my.data[, newcol.index] <- log10(sapply(my.data[, i], function(Col) { max(10, Col)}))
    colnames(my.data)[newcol.index] <- paste('log10_', colnames(my.data)[i], sep='')
  }
  my.data <- my.data[, -distance.vars]
  return (my.data)
}


log_transform_distances_dt <- function(my.data.table, desc.vars)
{
  stopifnot("data.table" %in% class(my.data.table))
  distance.vars <- grep("^m_to_", colnames(my.data.table), v=TRUE)
  if(any(distance.vars %in% desc.vars)){stop("distance variables in desc.vars. \n
	to avoid errors in the cleaning function, remove distance variables prior to modeling")}
  if(!length(grep("^m_to_", colnames(my.data.table)))) stop("no m_to_ variables in dataset. Have you already used log_transform_distances?")
  # throw out records with covariate values of -1
  my.data.table <- my.data.table[eval(parse(text=paste(distance.vars, "!= -1", collapse= " & "))),]
  for (i in distance.vars)
  {
    #my.data[, newcol.index] <- log10(sapply(my.data[, i], function(Col) { max(10, Col)}))
   my.data.table[, eval(paste('log10_', i, sep='')) := log10( sapply(my.data.table[[i]], function(Col){max(10,Col)}) )]
  }
   my.data.table[, eval(distance.vars) := NULL]
	return(my.data.table)

}


backtransform_lognames <- function(x){
	transformed <- grep("log10", x, v=TRUE)
	suffixes <- sapply(transformed, function(.v) strsplit(.v, '_')[[1]][-1])
	original <- sapply(suffixes, paste, collapse="_")
	union(setdiff(x, transformed), original)
}




#=============================================#
# makes new variables that are the sum of a3  #
# and a3 roads in buffers, removes individual #
# covariates for a2 and a3 roads in buffers   #
#=============================================#
combine_a23_ll <- function(my.data)
{
  a2.vars <- grep("ll_a2", colnames(my.data))
  for (i in a2.vars)
  {
    newcol.index <- 1 + length(colnames(my.data))
    a3.var <- grep(gsub('a2','a3',colnames(my.data)[i]), colnames(my.data))
    my.data[, newcol.index] <- my.data[, i] + my.data[, a3.var]
    colnames(my.data)[newcol.index] <- paste("ll_a23_", strsplit(colnames(my.data)[a3.var], '_')[[1]][3], sep="")
  }
  ll.vars <- grep("ll_a[^1]_s", colnames(my.data))
  my.data <- my.data[, -ll.vars]
  return(my.data)
}


combine_a23_ll_dt <- function(my.data.table, desc.vars)
{
  stopifnot("data.table" %in% class(my.data.table))
  a2.vars <- grep("ll_a2", colnames(my.data.table),v=TRUE)

  if(any(a2.vars %in% desc.vars)){stop("a2 variables in desc.vars. \n
	to avoid errors in the cleaning function, remove a2 variables prior to modeling")}

  if(length(grep("ll_a23", colnames(my.data.table)))) stop("ll_a23 variables already exist. Have you already used combine_a23_ll_dt?")
  for (i in a2.vars)
  {
    newcol.index <- 1 + length(colnames(my.data.table))
    a3.var <- grep(gsub('a2','a3',i), colnames(my.data.table),v=TRUE)
    my.data.table[, eval(paste("ll_a23_", strsplit(a3.var, '_')[[1]][3], sep="")) := my.data.table[[i]] + my.data.table[[a3.var]] ]
  }
  ll.vars <- grep("ll_a[^1]_s", colnames(my.data.table),v=TRUE)
  my.data.table <- my.data.table[, eval(ll.vars):=NULL]

  return(my.data.table)
}


backtransform_llnames <- function(x){
	transformed <- grep("ll_a23", x, v=TRUE)
	suffixes <- sapply(transformed, function(.v) strsplit(.v, '_')[[1]][3])
	original <- c(
		paste("ll_a2", suffixes, sep='_'),
		paste("ll_a3", suffixes, sep='_')
	)
	union(setdiff(x, transformed), original)
}


#========================================#
# returns random cross validation groups #
# with equal proportion per region       #
#========================================#
get_cv_groups <- function(national, my.segments = 10)
{
  rand.cv.groups <- rep(NA, national$obs)
  names(rand.cv.groups) <- national$monitors
  for (mons in list(national$weco.monitors, national$west.monitors, national$east.monitors))
  {
    obs <- length(mons)
    cv.groups <- trunc(my.segments * 1:obs / (obs+1)) + 1
    rand.cv.groups[as.character(mons)] <- cv.groups[order(runif(obs, 0, 1))]
  }
  return (rand.cv.groups)
}


#==================================================#
# the kriging code requires targets in each region #
# this function pads the X matrix with needed data #
#==================================================#
fill_X <- function(X, regions.missing, all.regions=c('east','weco','west'), pls.comps,uk.vars)
{
 pls.comps <- pls.comps + uk.vars
  if (!length(regions.missing))
    return (X)
  for (i in match(regions.missing, all.regions))
  {
    #X <- rbind(X, rep(0, length(all.regions)*pls.comps))
    X <- rbind(X, rep(0, length(all.regions)*(pls.comps+1)))
    #X[nrow(X), (i*pls.comps-(pls.comps-1)):(i*pls.comps)] = rep(.2, pls.comps)
    X[nrow(X), (i*(pls.comps+1)-pls.comps):(i*(pls.comps+1))] = c(1, rep(.2, pls.comps))
  }
  return (X)
}

fill_X_old <- function(X, regions.missing, all.regions=c('east','weco','west'), pls.comps)
{
  if (!length(regions.missing))
    return (X)
  for (i in match(regions.missing, all.regions))
  {
    #X <- rbind(X, rep(0, length(all.regions)*pls.comps))
    X <- rbind(X, rep(0, length(all.regions)*(pls.comps+1)))
    #X[nrow(X), (i*pls.comps-(pls.comps-1)):(i*pls.comps)] = rep(.2, pls.comps)
    X[nrow(X), (i*(pls.comps+1)-pls.comps):(i*(pls.comps+1))] = c(1, rep(.2, pls.comps))
  }
  return (X)
}



#==================================================#
# the kriging code requires targets in each region #
# this function pads the coords with needed data   #
#==================================================#
fill_coords <- function(coords, regions.missing, all.regions=c('east','weco','west'))
{
  if (!length(regions.missing))
    return (coords)
  for (i in match(regions.missing, all.regions))
    coords <- rbind(coords, c(-2000, 10))
  return (coords)
}


#====================================================#
# the kriging code requires targets in each region   #
# this function pads the region vec with needed data #
#====================================================#
fill_regional_vector <- function(vec, regions.missing, all.regions=c('east','weco','west'))
{
  if (!length(regions.missing))
    return (vec)
  return (c(vec, regions.missing))
}


