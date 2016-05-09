

#=======================================#
# returns a list of variables whose     #
# 10TH and 90TH quantiles are different #
#=======================================#

#' Variable Screening: Sufficient variability
#'
#' A function that screens out covariates with too little variability (10th percentile == 90th percentile)
#' @param my.data 
#' @param my.all.vars
#' @keywords 
#' @export
#' @examples 

fail_quantile_check <- function(my.data, my.all.vars)
{
  quantile.vars <- sapply(my.all.vars, function(Col) return (quantile(my.data[, Col], .10) == quantile(my.data[, Col], .90)))
  fail <- as.vector(sapply(names(quantile.vars[quantile.vars == TRUE]), function(Str) gsub('.10%', '', Str)))
  return (fail)
}


#==================================#
# returns a list of variables with #
# an absurdly skewed distrobution  #
#==================================#

#' Variable Screening: Skewness
#'
#' A function that screens out covariates that have extremely skewed distributions (maximum normalized value > 10).  Note: -10, too?
#' @param my.data 
#' @param my.all.vars
#' @keywords 
#' @export
#' @examples 

fail_skewed_distro <- function(my.data, my.all.vars)
{
  fail <- sapply(my.all.vars, function(Col) return (max(scale(my.data[, Col])) > 10))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}


#========================================#
# returns a list of land use variables 3 #
# whose max percentage is < 10 percent   #
#========================================#

#' Variable Screening: Rare land use types
#'
#' A function that screens out land use that never makes up more than 10% of the buffer area
#' @param my.data 
#' @param my.all.vars
#' @keywords 
#' @export
#' @examples 

fail_low_landuse <- function(my.data, my.all.vars)
{
  lu.vars <- grep("^lu_", my.all.vars, valu=T)
  fail <- sapply(lu.vars, function(Col) return (max(my.data[, Col]) < 10))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}


#==================================#
# returns a list of variabes with  #
# only 2 or less non-zero values   #
#==================================#

#' Variable Screening: All values are zero
#'
#' A function that screens out variables that are 0 for all points
#' @param my.data 
#' @param my.all.vars
#' @keywords 
#' @export
#' @examples 

fail_non_zero <- function(my.data, my.all.vars)
{
  fail <- sapply(my.all.vars, function(Col) return (sum(my.data[, Col] > 0) < 3))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}


