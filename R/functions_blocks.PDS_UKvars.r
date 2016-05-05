
make_reg_design <- function(my.object, region.vector, my.pls.comps, all.regions=c('east', 'weco', 'west'))
{
  my.pls.comps <- as.integer(my.pls.comps)
  total.vars <- my.pls.comps + ncol(my.object$UK.covars)
  n.regions <- length(all.regions)
  #design.matrix <- matrix(0, nrow=my.object$obs, ncol=n.regions * my.pls.comps)
  design.matrix <- matrix(0, nrow=my.object$obs, ncol=n.regions * (total.vars+1))
  j <- 1
  for (r in 1:n.regions)
  {
    indices <- region.vector == all.regions[r]
    #design.matrix[indices, j:(j+my.pls.comps-1)] <- my.object$pls$scores[indices, 1:my.pls.comps]
    if(sum(indices)>1){ 
	    design.matrix[indices, j:(j+total.vars)] <- cbind(1,my.object$PLS$scores[indices, 1:my.pls.comps],my.object$UK.covars[indices,])
    }else{
 	    if(sum(indices)==1){
              design.matrix[indices, j:(j+total.vars)] <- c(1,my.object$PLS$scores[indices, 1:my.pls.comps], my.object$UK.covars[indices,])
          }
    }
    #j <- j + my.pls.comps
    j <- j + total.vars+1
  }
  #design.matrix <- matrix(design.matrix, ncol=n.regions*my.pls.comps)
  rownames(design.matrix) <- my.object$monitors
  return (design.matrix)
}

#====================================================#
# returns a design matrix with national coefficients #
#====================================================#
make_nat_design_old <- function(my.object, region.vector, my.pls.comps)
{
 stop()
  design.matrix <- matrix(cbind(1,my.object$pls$scores[, 1:my.pls.comps]), ncol=(my.pls.comps+1))
  rownames(design.matrix) <- my.object$monitors
  return (design.matrix)
}


make_nat_design <- function(my.object, region.vector, my.pls.comps)
{
  design.matrix <- matrix(cbind(1,my.object$PLS$scores[, 1:my.pls.comps], my.object$UK.covars), ncol=(my.pls.comps+1+ncol(my.object$UK.covars)))
  rownames(design.matrix) <- my.object$monitors
  return (design.matrix)
}


#==============================================#
# returns a vector with elements corresponding #
# to which region each monitor is in           #
# aka v.hash                                   #
#==============================================#
make_reg_vector <- function(my.object)
{
  the.regions <- c(rep(NA, my.object$obs))
  the.regions[match(my.object$weco.monitors, my.object$monitors)] <- "weco"
  the.regions[match(my.object$west.monitors, my.object$monitors)] <- "west"
  the.regions[match(my.object$east.monitors, my.object$monitors)] <- "east"
  names(the.regions) <- my.object$monitors
  return (the.regions)
}


#====================================#
# returns a vector with all elements #
# set to be one region, namely 'all' #
#====================================#
make_nat_vector <- function(my.object)
{
  the.regions <- rep('all', my.object$obs)
  names(the.regions) <- my.object$monitors
  return (the.regions)
}


