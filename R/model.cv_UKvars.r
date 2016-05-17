
###note that additional effort has been taken to save components that may be later needed for graphing, etc
###specifically: parms, X, X.pred

#' Function for Running Cross-Validation
#'
#' Running cross-validation.
#' @param rawdata A matrix (?) of the rawdata, which must include the following elements: ???  Defaults to NULL.  This will terminate the function call.  
#' @param desc.vars A list of ?
#' @param pls.comps Integer indicating the number of PLS components.
#' @param UK.varnames A list of ?
#' @param manual.cv Defaults to NULL.  This will?
#' @param factr Defaults to 1e10.  What does it do?
#' @param verbose Boolean.  Defaults to TRUE.
#' @param my.segments Integer; number of segments?
#' @param regional Boolean indicating whether model is regional.  Defaults to TRUE.
#' @param kriging Boolean indicating whether to krig residuals.  Defaults to TRUE.
#' @return
#' @section To Do List:
#' I'm adding this section as a place to track things we might want to do as improvements
#' \describe{
#'   \item{}{}
#' }
#' @keywords cross-validation
#' @export
#' @examples 

PLSK.cv <- function(rawdata, desc.vars, pls.comps, UK.varnames=NULL, 
                    manual.cv=NULL, factr=1e9, verbose=FALSE, my.segments=10, regional=TRUE, kriging=TRUE){
  invars <- names(rawdata)
  pls.comps <-   as.integer(pls.comps)
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

  # Data pre-processing
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
  

  # set up the design matrix and region vector hashes
  b.hash <- list(make_reg_design, make_nat_design)
  names(b.hash) <- c('reg', 'nat')
  v.hash <- list(make_reg_vector, make_nat_vector)
  names(v.hash) <- c('reg', 'nat')
  
  # establish CV groups here
  weco.monitors <- poly.int2(weco.polygon, rawdata[, c('longitude','latitude')]) 
  east.monitors <- poly.int2(east.polygon, rawdata[, c('longitude','latitude')])
  west.monitors <- !east.monitors & !weco.monitors

  ###NOTE:one group can be entirely in one region but one region can't be made up of one group
  ###(if one region was made of one group,
  ###there would be no way to predict to that region when that group is left out)

  if(is.null(manual.cv)){
    cv.mx <- matrix(FALSE, nrow=nrow(rawdata), ncol=10)
    for (mons in list(weco.monitors, west.monitors, east.monitors)){
      obs <- sum(mons)
      cv.groups <- trunc(my.segments * 1:obs / (obs+1)) + 1
      cv.groups <- cv.groups[order(runif(obs, 0, 1))]
      for (i in 1:my.segments){
        cv.mx[which(mons)[cv.groups == i], i] = TRUE  #changed from grep(TRUE,mons) to which(mons) for readability
								    #pretty sure this could actually be written as mons & cv.groups==i
	}
    }
  }else{
    my.segments <- length(unique(manual.cv))
    cv.mx <- matrix(FALSE, nrow=nrow(rawdata), ncol=my.segments)
    for(i in 1:my.segments){
      cv.mx[manual.cv==i,i] <- TRUE
    }
    for (mons in list(weco.monitors, west.monitors, east.monitors)){
	stopifnot(sum(as.logical(colSums(cv.mx[mons,])))>1) ##check to see that every region contains multiple groups
    }
  }
  
  stopifnot(all(rowSums(cv.mx)==1))  
  
  cv.pred <- rep(NA, l=nrow(rawdata))
  cv.var <- rep(NA, l=nrow(rawdata))
  names(cv.pred) <- rawdata$native_id
  names(cv.var) <- rawdata$native_id
  
  parms <- list()
  X <- list()
  X.pred <- list()
 

  for (i in 1:my.segments)
  {
    print (i)
    #
    # prune out variables
    cv.rawdata <- rawdata[!cv.mx[,i],] 
    #cv.rawdata is the "90%" sample
    if(verbose){print(paste("dim(cv.rawdata):", dim(cv.rawdata)[1],dim(cv.rawdata)[2]))}

    all.vars <- colnames(cv.rawdata)[! colnames(cv.rawdata) %in% desc.vars]
    exclude.vars <- desc.vars
    exclude.vars <- c(exclude.vars, fail_quantile_check(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_skewed_distro(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_low_landuse(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_non_zero(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars,c("lambert_x","lambert_y"))
    exclude.vars <- unique(exclude.vars)
    exclude.vars <- exclude.vars[exclude.vars %in% colnames(cv.rawdata)]
    
    if(!is.null(UK.varnames)) if(any(UK.varnames %in% exclude.vars)){
      stop("variable/s ", instersect(UK.varnames,exclude.vars), " have been excluded, probably due to failing cleaning checks")
    }

    #
    # create national R object out of no2 object
    model.obj <- make_modeling_object(rawdata=cv.rawdata, exclude.vars=exclude.vars, UK.varnames=UK.varnames)
    
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
    region.vec   <- v.hash[[v.type]](model.obj)
    model.obj$X  <- b.hash[[b.type]](model.obj, b.region.vec, pls.comps)
    rownames(model.obj$X) <- model.obj$monitors
    X[[i]] <- model.obj$X
    
    if(kriging){
      ini.l.pars <- c(0, log(.1), log(750), 0, log(.1), log(10), 0, log(.1), log(60))
      model.obj$parms <- my.likfit(ini.l.pars = ini.l.pars[1:(3*length(unique(region.vec)))],
  		X = model.obj$X,
		coords = model.obj$coords,
 		reg.ind = region.vec,
		data = model.obj$y,trace=verbose,factr=factr)

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
    }else{ #end "if kriging
	    lmdata <- data.frame(y=model.obj$y, 	as.data.frame(model.obj$X))
	    lmRHS <- paste0(paste0("V", 1:ncol(model.obj$X),collapse="+"), "-1")
	    lmformula <- 	as.formula(c("y~", lmRHS ))

	    #model.obj$parms <- lm(model.obj$y~model.obj$X-1)
	    model.obj$parms <- lm(eval(lmformula), data=lmdata)
	    parms[[i]] <- model.obj$parms
    }


    parms[[i]] <- model.obj$parms

    # predict on the leftout group
    ppts <- make_predict_object(rawdata=rawdata[cv.mx[,i],], model.obj, pls.comps)
    b.region.vec    <- v.hash[[b.type]](ppts)
    ppts.region.vec <- v.hash[[v.type]](ppts)
    modl.region.vec <- v.hash[[v.type]](model.obj)
    ppts$X          <- b.hash[[b.type]](ppts, b.region.vec, pls.comps)

    X.pred[[i]] <- ppts$X
    rownames(X.pred[[i]]) <- rawdata[cv.mx[,i],"native_id"]

    if(kriging){
      mark <- c.exp(model.obj$parms,
                  model.obj$X, model.obj$coords, modl.region.vec, model.obj$y,
                  ppts$X, ppts$coords, ppts.region.vec)
      cv.pred[cv.mx[,i]] <- mark[, 1]
      cv.var[cv.mx[,i]] <- mark[, 2]
    }else{
      cv.pred[cv.mx[,i]] <- predict(model.obj$parms, as.data.frame(ppts$X))
    }

    #
    # save the R workspace
 
    print (i)
  } 
  as.list(environment())
}







PLSnoK.cv <- function(rawdata, desc.vars, pls.comps, UK.varnames=NULL, manual.cv=NULL, factr=1e9, verbose=FALSE, my.segments=10, regional=TRUE){
  invars <- names(rawdata)
  pls.comps <-   as.integer(pls.comps)
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
  



  
  # set up the design matrix and region vector hashes
  b.hash <- list(make_reg_design, make_nat_design)
  names(b.hash) <- c('reg', 'nat')
  v.hash <- list(make_reg_vector, make_nat_vector)
  names(v.hash) <- c('reg', 'nat')
  
  
  # establish CV groups here
  weco.monitors <- poly.int2(weco.polygon, rawdata[, c('longitude','latitude')]) 
  east.monitors <- poly.int2(east.polygon, rawdata[, c('longitude','latitude')])
  west.monitors <- !east.monitors & !weco.monitors

	###NOTE:one group can be entirely in one region but one region can't be made up of one group
	###(if one region was made of one group,
	###there would be no way to predict to that region when that group is left out)

  if(is.null(manual.cv)){
  cv.mx <- matrix(FALSE, nrow=nrow(rawdata), ncol=10)
  for (mons in list(weco.monitors, west.monitors, east.monitors))
  {
    obs <- sum(mons)
    cv.groups <- trunc(my.segments * 1:obs / (obs+1)) + 1
    cv.groups <- cv.groups[order(runif(obs, 0, 1))]
    for (i in 1:my.segments){
      cv.mx[which(mons)[cv.groups == i], i] = TRUE  #changed from grep(TRUE,mons) to which(mons) for readability
								    #pretty sure this could actually be written as mons & cv.groups==i
	}
  }
  }else{
      my.segments <- length(unique(manual.cv))
	cv.mx <- matrix(FALSE, nrow=nrow(rawdata), ncol=my.segments)
	for(i in 1:my.segments){
		cv.mx[manual.cv==i,i] <- TRUE
	}
	for (mons in list(weco.monitors, west.monitors, east.monitors)){
		stopifnot(sum(as.logical(colSums(cv.mx[mons,])))>1) ##check to see that every region contains multiple groups
	}
  }
  
  stopifnot(all(rowSums(cv.mx)==1))  
  
  cv.pred <- rep(NA, l=nrow(rawdata))
  cv.var <- rep(NA, l=nrow(rawdata))
  names(cv.pred) <- rawdata$native_id
  names(cv.var) <- rawdata$native_id
  
  parms <- list()
  X <- list()
  X.pred <- list()
 

  for (i in 1:my.segments)
  {
    print (i)
    #
    # prune out variables
    cv.rawdata <- rawdata[!cv.mx[,i],] 
    #cv.rawdata is the "90%" sample
    if(verbose){print(paste("dim(cv.rawdata):", dim(cv.rawdata)[1],dim(cv.rawdata)[2]))}

    all.vars <- colnames(cv.rawdata)[! colnames(cv.rawdata) %in% desc.vars]
    exclude.vars <- desc.vars
    exclude.vars <- c(exclude.vars, fail_quantile_check(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_skewed_distro(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_low_landuse(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars, fail_non_zero(cv.rawdata, all.vars))
    exclude.vars <- c(exclude.vars,c("lambert_x","lambert_y"))
    exclude.vars <- unique(exclude.vars)
    exclude.vars <- exclude.vars[exclude.vars %in% colnames(cv.rawdata)]
    
  if(!is.null(UK.varnames)) if(any(UK.varnames %in% exclude.vars)){
	stop("variable/s ", instersect(UK.varnames,exclude.vars), " have been excluded, probably due to failing cleaning checks")
  }

    #
    # create national R object out of no2 object
    model.obj <- make_modeling_object(rawdata=cv.rawdata, exclude.vars=exclude.vars, UK.varnames=UK.varnames)
    
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
    model.obj$X <- b.hash[[b.type]](model.obj, b.region.vec, pls.comps)
    rownames(model.obj$X) <- model.obj$monitors
    X[[i]] <- model.obj$X

    lmdata <- data.frame(y=model.obj$y, 	as.data.frame(model.obj$X))
     lmRHS <- paste0(paste0("V", 1:ncol(model.obj$X),collapse="+"), "-1")
    lmformula <- 	as.formula(c("y~", lmRHS ))

    #model.obj$parms <- lm(model.obj$y~model.obj$X-1)
    model.obj$parms <- lm(eval(lmformula), data=lmdata)
    parms[[i]] <- model.obj$parms



    # predict on the leftout group
    ppts <- make_predict_object(rawdata=rawdata[cv.mx[,i],], model.obj, pls.comps)
    b.region.vec <- v.hash[[b.type]](ppts)
    ppts.region.vec <- v.hash[[v.type]](ppts)
    modl.region.vec <- v.hash[[v.type]](model.obj)
    ppts$X <- b.hash[[b.type]](ppts, b.region.vec, pls.comps)

    X.pred[[i]] <- ppts$X
    rownames(X.pred[[i]]) <- rawdata[cv.mx[,i],"native_id"]




    cv.pred[cv.mx[,i]] <- predict(model.obj$parms, as.data.frame(ppts$X))

#    cv.var[cv.mx[,i]] <- mark[, 2]
    
    #
    # save the R workspace
 
    print (i)
  } 

	as.list(environment())
}




###note that additional effort has been taken to save components that may be later needed for graphing, etc
###specifically: parms, X, X.pred

#' Function for Running Cross-Validation
#'
#' Running cross-validation.
#' @param rawdata A matrix (?) of the rawdata, which must include the following elements: ???  Defaults to NULL.  This will terminate the function call.  
#' @param desc.vars A list of ?
#' @param pls.comps Integer indicating the number of PLS components.
#' @param UK.varnames A list of ?
#' @param manual.cv Defaults to NULL.  This will?
#' @param factr Defaults to 1e10.  What does it do?
#' @param verbose Boolean.  Defaults to TRUE.
#' @param my.segments Integer; number of segments?
#' @param regional Boolean indicating whether model is regional.  Defaults to TRUE.
#' @param kriging Boolean indicating whether to krig residuals.  Defaults to TRUE.
#' @return 
#' @keywords cross-validation
#' @export
#' @examples 

PLSK.cv.alt <- function(rawdata, desc.vars, pls.comps, UK.varnames=NULL, 
                    manual.cv=NULL, factr=1e9, verbose=FALSE, my.segments=10, regional=TRUE, kriging=TRUE){
  
  # establish CV groups here
  weco.monitors <- poly.int2(weco.polygon, rawdata[, c('longitude','latitude')]) 
  east.monitors <- poly.int2(east.polygon, rawdata[, c('longitude','latitude')])
  west.monitors <- !east.monitors & !weco.monitors

  ###NOTE:one group can be entirely in one region but one region can't be made up of one group
  ###(if one region was made of one group,
  ###there would be no way to predict to that region when that group is left out)

  if(is.null(manual.cv)){
    cv.mx <- matrix(FALSE, nrow=nrow(rawdata), ncol=10)
    for (mons in list(weco.monitors, west.monitors, east.monitors)){
      obs <- sum(mons)
      cv.groups <- trunc(my.segments * 1:obs / (obs+1)) + 1
      cv.groups <- cv.groups[order(runif(obs, 0, 1))]
      for (i in 1:my.segments){
        cv.mx[which(mons)[cv.groups == i], i] = TRUE  #changed from grep(TRUE,mons) to which(mons) for readability
								    #pretty sure this could actually be written as mons & cv.groups==i
	}
    }
  }else{
    my.segments <- length(unique(manual.cv))
    cv.mx <- matrix(FALSE, nrow=nrow(rawdata), ncol=my.segments)
    for(i in 1:my.segments){
      cv.mx[manual.cv==i,i] <- TRUE
    }
    for (mons in list(weco.monitors, west.monitors, east.monitors)){
	stopifnot(sum(as.logical(colSums(cv.mx[mons,])))>1) ##check to see that every region contains multiple groups
    }
  }
  
  stopifnot(all(rowSums(cv.mx)==1))  
  
  cv.pred <- rep(NA, l=nrow(rawdata))
  cv.var <- rep(NA, l=nrow(rawdata))
  names(cv.pred) <- rawdata$native_id
  names(cv.var) <- rawdata$native_id
  
  parms <- list()
  X <- list()
  X.pred <- list()
 

  for (i in 1:my.segments)
  {
    print (i)
    #
    # prune out variables
    cv.rawdata <- rawdata[!cv.mx[,i],] 
    #cv.rawdata is the "90%" sample
    if(verbose){print(paste("dim(cv.rawdata):", dim(cv.rawdata)[1],dim(cv.rawdata)[2]))}

   
    if(kriging){
      if(regional){
        result.temp <- PLSK.full(cv.rawdata, desc.vars, pls.comps, UK.varnames, factr, verbose, regional)
        model.obj$parms <- result.temp$model.obj$parms
	  }
    }else{ #end "if kriging
           # data cleaning needed here?
	    lmdata <- data.frame(y=model.obj$y, 	as.data.frame(model.obj$X))
	    lmRHS <- paste0(paste0("V", 1:ncol(model.obj$X),collapse="+"), "-1")
	    lmformula <- 	as.formula(c("y~", lmRHS ))

	    #model.obj$parms <- lm(model.obj$y~model.obj$X-1)
	    model.obj$parms <- lm(eval(lmformula), data=lmdata)
	    parms[[i]] <- model.obj$parms
    }


    parms[[i]] <- model.obj$parms

    # predict on the leftout group
    ppts <- make_predict_object(rawdata=rawdata[cv.mx[,i],], model.obj, pls.comps)
    b.region.vec    <- v.hash[[b.type]](ppts)
    ppts.region.vec <- v.hash[[v.type]](ppts)
    modl.region.vec <- v.hash[[v.type]](model.obj)
    ppts$X          <- b.hash[[b.type]](ppts, b.region.vec, pls.comps)

    X.pred[[i]] <- ppts$X
    rownames(X.pred[[i]]) <- rawdata[cv.mx[,i],"native_id"]

    if(kriging){
      mark <- c.exp(model.obj$parms,
                  model.obj$X, model.obj$coords, modl.region.vec, model.obj$y,
                  ppts$X, ppts$coords, ppts.region.vec)
      cv.pred[cv.mx[,i]] <- mark[, 1]
      cv.var[cv.mx[,i]] <- mark[, 2]
    }else{
      cv.pred[cv.mx[,i]] <- predict(model.obj$parms, as.data.frame(ppts$X))
    }

    #
    # save the R workspace
 
    print (i)
  } 
  as.list(environment())
}






