

#=======================================#
# returns a list of variables whose     #
# 10TH and 90TH quantiles are different #
#=======================================#

#' Variable Screening: Sufficient variability
#'
#' A function that screens out covariates with too little variability (10th percentile == 90th percentile)
#' @param my.data See rawdata under PLSK.full
#' @param my.all.vars String vector listing candidate variables; often, rownames(my.data)
#' @return String vector listing variables that failed the cleaning check
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
#' @param my.data See rawdata under PLSK.full
#' @param my.all.vars String vector listing candidate variables; often, rownames(my.data)
#' @return String vector listing variables that failed the cleaning check
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
#' A function that screens out land use that never makes up more than 10% of the buffer area.  Searches my.all.vars for 
#' variable names that start with "lu_".
#' @param my.data See rawdata under PLSK.full
#' @param my.all.vars String vector listing candidate variables; often, rownames(my.data)
#' @param lu.search A regular expression (string) that is passed to grep and identifies the land use variables.  Default is "^lu_".
#' @return String vector listing variables that failed the cleaning check
#' @keywords 
#' @export
#' @examples 

fail_low_landuse <- function(my.data, my.all.vars, lu.search = "^lu_")
{
  lu.vars <- grep(lu.search, my.all.vars, valu=T)
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
#' @param my.data See rawdata under PLSK.full
#' @param my.all.vars String vector listing candidate variables; often, rownames(my.data)
#' @return String vector listing variables that failed the cleaning check
#' @keywords 
#' @export
#' @examples 

fail_non_zero <- function(my.data, my.all.vars)
{
  fail <- sapply(my.all.vars, function(Col) return (sum(my.data[, Col] > 0) < 3))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}



#============================================#
# transforms all distance variables to log10 #
# and lower bounds them to 10 meters         #
#============================================#

#' Variable Transformation: Log Transform Distances
#'
#' Log10 transformation of variables identified by regular expression (search string)
#' @param my.data See rawdata under PLSK.full
#' @param my.all.vars String vector listing candidate variables; often, rownames(my.data)
#' @param dist.search A regular expression (string) that is passed to grep and identifies the distance variables.  Default is "^m_to_".
#' @return A data table with variables transformed and re-named.
#' @keywords 
#' @export
#' @examples 



log_transform_distances_dt <- function(my.data.table, desc.vars, dist.search = "^m_to_")
{
  distance.vars <- grep(dist.search, colnames(my.data.table), v=TRUE)

  # Error trap
  stopifnot("data.table" %in% class(my.data.table))
  if(any(distance.vars %in% desc.vars)){
      stop("distance variables in desc.vars. \n
	to avoid errors in the cleaning function, remove distance variables prior to modeling")
  }
  if(!length(grep(dist.search, colnames(my.data.table)))) stop("no m_to_ variables in dataset. Have you already used log_transform_distances?")

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


#=============================================#
# Data cleaning                               #
#=============================================#

clean_data <- function(my.data, desc.vars){
  my.data <- log_transform_distances_dt(my.data, desc.vars=desc.vars)
  my.data <- combine_a23_ll_dt(my.data, desc.vars=desc.vars)
  my.data <- as.data.frame(my.data)

  desc.vars <- union(desc.vars, 
     c('pollutant_conc','location_id','native_id','latitude','longitude', 'monitor_type')
  )
   
  # Variable screening
  all.vars <- colnames(my.data)[! colnames(my.data) %in% desc.vars]
  exclude.vars <- desc.vars
  exclude.vars <- c(exclude.vars, fail_quantile_check(my.data, all.vars))
  exclude.vars <- c(exclude.vars, fail_skewed_distro(my.data, all.vars))
  exclude.vars <- c(exclude.vars, fail_low_landuse(my.data, all.vars))
  exclude.vars <- c(exclude.vars, fail_non_zero(my.data, all.vars))
  exclude.vars <- c(exclude.vars,c("lambert_x","lambert_y"))
  exclude.vars <- unique(exclude.vars)
  exclude.vars <- exclude.vars[exclude.vars %in% colnames(my.data)]

  my.data <- my.data[, !(colnames(my.data) %in% exclude.vars)]
 
  # Return data and variables chosen for exclusion
  list("rawdata" = my.data, "exclude.vars" = exclude.vars)
}




