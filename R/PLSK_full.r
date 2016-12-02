

#' Function for Running a Full Model
#'
#' Function for running a full model. 
#' @param rawdata A data.frame of the rawdata.  See rawdata section for details. Failing to supply rawdata will terminate the function call.  
#' @param desc.vars A vector of ?
#' @param pls.comps An integer between 1 and 5 indicating the number of PLS components.
#' @param UK.varnames Defaults to NULL.  These are covariates to be used in universal kriging.
#' @param factr Defaults to 1e9.
#' @param verbose Defaults to FALSE.  This will?
#' @param regional Defaults to TRUE.  This will?
#' @param clean.data Boolean indicating whether to apply data cleaning functions.  Defaults to TRUE.
#' @details This is one of the workhorse functions of this package.
#' By default, distance variables are log-transformed.  Geographic variables are screened for sufficient variability, etc. and other stuff I'll specify later.
#' The function is expecting something like the following to be in the rawdata object:
#' \describe{
#'   \item{Location-specifying columns}{location_id, native_id, latitude, longitude, monitor_type, lambert_x, lambert_y, state_plane, state, county}
#'   \item{Pollution concentrations}{pollutant_conc}
#' }
#' @return A list containing a lot of stuff?
#' \describe{
#'   \item{ini.l.pars}{Vector of ?.  NOTE: this is hard-coded, can it be removed from the output?}
#'   \item{region.vec}{List of vector of named strings, indicating which monitors are in which regions}
#'   \item{b.region.vec}{List of vector of named strings, indicating which monitors are in which regions}
#'   \item{v.type}{Character indicating variogram type.  Options include "reg" for regional, or "nat" for national.}
#'   \item{b.type}{Character indicating beta type.  Options include "reg" for regional, or "nat" for national.}
#'   \item{v.hash}{List of 2 objects: reg:function and nat:function.  These are functions that have something to do with the variogram.}
#'   \item{b.hash}{List of 2 objects: reg:function and nat:function.  These are functions that have something to do with the betas.}
#'   \item{model.obj}{List of 13.  The modeling object, generated by make_modeling_object}
#'   \item{exclude.vars}{List of string vector listing covariates}
#'   \item{all.vars}{List of string vector listing covariates}
#'   \item{miss}{List of Boolean vector indicating ?}
#'   \item{year}{List of integer indicating year of model}
#'   \item{pollutant}{List of string indicating pollutant}
#'   \item{invars}{List of string vector, listing all variables passed to modeling process.  NOTE since this is just rownames(rawdata), can it be removed?}
#'   \item{rawdata}{List of data.frame inluding all original data}
#'   \item{desc.vars}{List of string vector, listing ?}
#'   \item{pls.comps}{List of integer indicating number of PLS components}
#'   \item{UK.varnames}{Stuff}
#'   \item{factr}{List of ?}
#'   \item{verbose}{List of boolean; indicates?}
#'   \item{regional}{List of boolean; if true then regional, if false then single region}
#' }
#' @section Modeling Object:
#' The modeling object (model.obj) is a list of 13
#' \describe{
#'   \item{y}{Vector of named numbers (concentrations, I assume)}
#'   \item{UK.covars}{Stuff}
#'   \item{covars.pls}{A matrix of something ?}
#'   \item{PLS}{List of 19}
#'   \item{obs}{An integer indicating the number of observations (I think)}
#'   \item{coords}{A data.frame with x and y coordinates (lambert projection)}
#'   \item{gps}{A data.frame with lat and long coordinates}
#'   \item{monitors}{String vector listing monitors}
#'   \item{weco.monitors}{String vector listing monitors in the weco region}
#'   \item{east.monitors}{String vector listing monitors in the east region}
#'   \item{west.monitors}{String vector listing monitors in the west region}
#'   \item{X}{A matrix with rows equal to the number of monitors and columns?}
#'   \item{parms}{List of 5, generated from my.likfit}
#' }
#' @section PLS:
#' The PLS object is a list of 19.  In general below, m is the number of monitors.
#' \describe{
#'   \item{coefficients}{List of 3-dimensional matrix (number of... monitors plus 32?, y, and 5 components?)}
#'   \item{scores}{Matrix that is... monitors by 5 components}
#'   \item{loadings}{Matrix that is number of monitors plus 32 by 5}
#'   \item{loading.weights}{a similar matrix?}
#'   \item{Yscores}{Matrix that is number of m x 5}
#'   \item{Yloadings}{Matrix 1 x 5 that?}
#'   \item{projection}{366 by 5 matrix}
#'   \item{Xmeans}{numeric vector}
#'   \item{Ymeans}{Number indicating?}
#'   \item{fitted.values}{334  by 5 matrix}
#'   \item{residuals}{34x1x5 matrix}
#'   \item{Xvar}{Vector of 5 named numbers}
#'   \item{Xtotvar}{Number indicating?}
#'   \item{fit.time}{Named number indicating?}
#'   \item{ncomp}{Integer.... 5}
#'   \item{method}{"kernelpls"}
#'   \item{call}{language plsr(formula = y ~ covars.pls, ncomp = pls.comps, data = model.obj, validation = "none")}
#'   \item{terms}{Classes 'terms', 'formula' length 3 y ~ covars.pls}
#'   \item{model}{data frame of y and covars.pls, which is 334 by 366?}
#' }
#' @section Parms:
#' The parms object is a list of 5, and is generated by my.likfit.  In general below, m is the number of monitors.
#' \describe{
#'   \item{log.cov.pars}{List of vector of named numbers; includes tau, sigma and rho for each region (variogram parms?)}
#'   \item{beta}{Matrix (instead of vector?) intercept and betas (b1 and b2?) for each region}
#'   \item{max.log.lik}{A list of a number}
#'   \item{hess.pd}{A list of a boolean indicating whether the hessian was positive definite (?)}
#'   \item{hessian}{A matrix (9x9) which is the Hessian.}
#' }
#' @section To Do List:
#' I'm adding this section as a place to track things we might want to do as improvements
#' \describe{
#'   \item{Hard-coded region parameters}{Re-work for flexible regions}
#' }
#' @keywords 
#' @importFrom ruk rlikfit
#' @examples
#' @export

PLSK.full <- function(rawdata, desc.vars, pls.comps, UK.varnames=NULL, factr=1e9, verbose=FALSE, regional=TRUE, clean.data=TRUE, debug.mode=FALSE){

  # Which variables came in with the raw data
  invars    <- names(rawdata)
  
  # Size of raw dataset 
  print(paste("dim(rawdata):", dim(rawdata)[1],dim(rawdata)[2]))

  # Stop if there are duplicated locations
  stopifnot(sum(duplicated(rawdata$location_id))==0)

  # Number of PLS components
  pls.comps <- as.integer(pls.comps)
  if(pls.comps > 5) stop("PLS components must be <= 5. this is hardcoded somewhere. can be changed if necessary")
  print(paste("n.comps:", pls.comps))

  # Pollutant and year of model
  pollutant <- attr(rawdata, "pollutant")
  year      <- attr(rawdata, "year")
  if(is.null(pollutant) | is.null(year)){
	warning("attr(rawdata,'pollutant') or attr(rawdata, 'year') is null. 
		If you set these values manually prior to modeling, 
		they will be stored for later reference in the model object as object$pollutant, object$year"
	)
  }

  # Descriptive variables to retain
  desc.vars <- union(desc.vars, 
                   c('pollutant_conc','location_id','native_id','latitude','longitude', 'monitor_type'))

  # Data pre-processing
  rawdata <- as.data.table(rawdata)
  miss <- rawdata[, sapply(.SD, is.na),.SDcols=setdiff(names(rawdata),desc.vars)]
  if(sum(miss)!=0){
	cat("rawdata contains", sum(as.logical(rowSums(miss))), "rows with missing values.\n")
	cat("rawdata contains", sum(as.logical(colSums(miss))), "columns with missing values.\n")
	stop("remove missing columns or rows prior to modeling")
  }

  ## DATA CLEANING ##
  
  rawdata <- log_transform_distances_dt(rawdata, desc.vars=desc.vars)
  rawdata <- combine_a23_ll_dt(rawdata,desc.vars=desc.vars)
  rawdata <- as.data.frame(rawdata)
  
  # Variable screening
  all.vars     <- colnames(rawdata)[! colnames(rawdata) %in% desc.vars]
  exclude.vars <- desc.vars
  if (clean.data){
    exclude.vars <- c(exclude.vars, fail_quantile_check(rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_skewed_distro(rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_low_landuse(rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_non_zero(rawdata, all.vars))
    exclude.vars <- c(exclude.vars,c("lambert_x","lambert_y"))
  }
  exclude.vars <- unique(exclude.vars)
  if (debug.mode){ print(exclude.vars) }
  exclude.vars <- exclude.vars[exclude.vars %in% colnames(rawdata)]

  if(!is.null(UK.varnames)) if(any(UK.varnames %in% exclude.vars)){
	stop("variable/s ", instersect(UK.varnames,exclude.vars), " have been excluded, probably due to failing cleaning checks")
  }
 
  ## END DATA CLEANING ##


  # create national R object out of rawdata object
  model.obj <- make_modeling_object(rawdata, exclude.vars=exclude.vars, UK.varnames=UK.varnames)

  # set up the design matrix and region vector hashes
  b.hash <- list(make_reg_design, make_nat_design)
  names(b.hash) <- c('reg', 'nat')
  v.hash <- list(make_reg_vector, make_nat_vector)
  names(v.hash) <- c('reg', 'nat')
  
  #
  # the betas are regional, the variograms are too!
  # NOTE: if betas and variograms are either both regional or national, why 2 params instead of 1?
    if(regional){
      b.type <- 'reg'
      v.type <- 'reg'
    }else{
      b.type <- 'nat'
      v.type <- 'nat'
    }
  #### Region vectors
  b.region.vec <- v.hash[[b.type]](model.obj)
  region.vec   <- v.hash[[v.type]](model.obj)
  
  # Design matrix
  # Generally make_reg_design(model.obj, b.region.vec, pls.comps)
  model.obj$X <- b.hash[[b.type]](my.object=model.obj, region.vector=b.region.vec, my.pls.comps=pls.comps) 
  rownames(model.obj$X) <- model.obj$monitors

  # Workhorse
  ini.l.pars <- c(0, log(.1), log(750), 0, log(.1), log(10), 0, log(.1), log(60))
  model.obj$parms <- rlikfit(y=model.obj$y,
  				     X=model.obj$X,
  				     coords=model.obj$coords,
  				     reg.ind=region.vec,
  				     init.pars= ini.l.pars[1:(3*length(unique(region.vec)))],
                             optim.args=list(control=list(trace= verbose,factr=factr)))
				
  if(regional){
    rownames(model.obj$parms$beta) <-    c("east_intercept",
      paste("east_b", c(1:pls.comps,UK.varnames), sep=''), 
      "weco_intercept",
      paste("weco_b", c(1:pls.comps,UK.varnames), sep=''),
      "west_intercept",
      paste("west_b", c(1:pls.comps,UK.varnames), sep=''))
    
    names(model.obj$parms$log.cov.pars) <- c('east_tau',
                                           'east_sigma',
                                           'east_rho',
                                           'weco_tau',
                                           'weco_sigma',
                                           'weco_rho',
                                           'west_tau',
                                           'west_sigma',
                                           'west_rho')
  }
  #rm(c("invars", "ini.l.pars"))
  as.list(environment())
}
