

#=======================================#
# returns a list of variables whose     #
# 10TH and 90TH quantiles are different #
#=======================================#
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
fail_non_zero <- function(my.data, my.all.vars)
{
  fail <- sapply(my.all.vars, function(Col) return (sum(my.data[, Col] > 0) < 3))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}


