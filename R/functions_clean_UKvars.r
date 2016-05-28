

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


