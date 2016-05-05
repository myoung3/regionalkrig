


PLSK.full <- function(rawdata, desc.vars, pls.comps,UK.varnames=NULL, factr=1e9, verbose=FALSE,regional=TRUE){
  
  invars <- names(rawdata) 
  pls.comps <-   as.integer(pls.comps)
  if(pls.comps > 5) stop("PLS components must be <= 5. this is hardcoded somewhere. can be changed if necessary")
  pollutant <- attr(rawdata, "pollutant")
  year <- attr(rawdata, "year")
  if(is.null(pollutant) | is.null(year)){
	warning("attr(rawdata,'pollutant') or attr(rawdata, 'year') is null. 
		If you set these values manually prior to modeling, 
		they will be stored for later reference in the model object as object$pollutant, object$year"
	)
  }
  print(paste("n.comps:", pls.comps))
  print(paste("dim(rawdata):", dim(rawdata)[1],dim(rawdata)[2]))

  stopifnot(sum(duplicated(rawdata$location_id))==0)
        rawdata <- as.data.table(rawdata)
	  miss <- rawdata[, sapply(.SD, is.na),.SDcols=setdiff(names(rawdata),desc.vars)]
	  if(sum(miss)!=0){
		cat("rawdata contains", sum(as.logical(rowSums(miss))), "rows with missing values.\n")
		cat("rawdata contains", sum(as.logical(colSums(miss))), "columns with missing values.\n")
		stop("remove missing columns or rows prior to modeling")
	  }
	
	  rawdata <- log_transform_distances_dt(rawdata, desc.vars=desc.vars)
	  rawdata <- combine_a23_ll_dt(rawdata,desc.vars=desc.vars)
	  rawdata <- as.data.frame(rawdata)

  desc.vars <- union(desc.vars, 
     c('pollutant_conc','location_id','native_id','latitude','longitude', 'monitor_type')
  )
    
  all.vars <- colnames(rawdata)[! colnames(rawdata) %in% desc.vars]
  exclude.vars <- desc.vars
  exclude.vars <- c(exclude.vars, fail_quantile_check(rawdata, all.vars))
  exclude.vars <- c(exclude.vars, fail_skewed_distro(rawdata, all.vars))
  exclude.vars <- c(exclude.vars, fail_low_landuse(rawdata, all.vars))
  exclude.vars <- c(exclude.vars, fail_non_zero(rawdata, all.vars))
  exclude.vars <- c(exclude.vars,c("lambert_x","lambert_y"))
  exclude.vars <- unique(exclude.vars)
  exclude.vars <- exclude.vars[exclude.vars %in% colnames(rawdata)]

  if(!is.null(UK.varnames)) if(any(UK.varnames %in% exclude.vars)){
	stop("variable/s ", instersect(UK.varnames,exclude.vars), " have been excluded, probably due to failing cleaning checks")
  }

  # create national R object out of rawdata object
  model.obj <- make_modeling_object(rawdata, exclude.vars=exclude.vars,UK.varnames=UK.varnames)


  # set up the design matrix and region vector hashes
  b.hash <- list(make_reg_design, make_nat_design)
  names(b.hash) <- c('reg', 'nat')
  v.hash <- list(make_reg_vector, make_nat_vector)
  names(v.hash) <- c('reg', 'nat')
  
  #
  # the betas are regional, the variograms are too!
    if(regional){
	    b.type <- 'reg'
	    v.type <- 'reg'
    }else{
	    b.type <- 'nat'
	    v.type <- 'nat'
    }
  b.region.vec <- v.hash[[b.type]](model.obj)
  region.vec <- v.hash[[v.type]](model.obj)
  
  model.obj$X <- b.hash[[b.type]](my.object=model.obj, region.vector=b.region.vec, my.pls.comps=pls.comps) #generally make_reg_design(model.obj, b.region.vec, pls.comps)
  rownames(model.obj$X) <- model.obj$monitors
  ini.l.pars <- c(0, log(.1), log(750), 0, log(.1), log(10), 0, log(.1), log(60))
  model.obj$parms <- my.likfit(ini.l.pars = ini.l.pars[1:(3*length(unique(region.vec)))],
                               X = model.obj$X,
                               coords = model.obj$coords,
                               reg.ind = region.vec,
                               data = model.obj$y,
								trace=verbose,
								factr=factr)
  if(regional){
	 rownames(model.obj$parms$beta) <-    c("east_intercept",
				paste("east_b",c(1:pls.comps,UK.varnames), sep=''), 
				"weco_intercept",
				paste("weco_b",c(1:pls.comps,UK.varnames),sep=''),
				"west_intercept",
				paste("west_b",c(1:pls.comps,UK.varnames), sep=''))
  
  
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
  
  	as.list(environment())
}
