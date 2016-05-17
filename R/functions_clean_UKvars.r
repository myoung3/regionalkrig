

#============================================#
# transforms all distance variables to log10 #
# and lower bounds them to 10 meters         #
#============================================#


log_transform_distances_dt <- function(my.data.table, desc.vars)
{
  stopifnot("data.table" %in% class(my.data.table))
  distance.vars <- grep("^m_to_", colnames(my.data.table), v=TRUE)
  if(any(distance.vars %in% desc.vars)){
      stop("distance variables in desc.vars. \n
	to avoid errors in the cleaning function, remove distance variables prior to modeling")
  }
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


