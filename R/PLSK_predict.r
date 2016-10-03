
# ================================= #
# Why are there 3 versions of this? #
# ================================= #

#' PLS Prediction
#'
#' A function that 
#' @param model
#' @param rawdata See rawdata under PLSK.full
#' @param desc.vars
#' @param debug.mode
#' @return a list
#' \describe{
#'   \item{mark}{}
#'   \item{modl.region.vec}{}
#'   \item{ppts.region.vec}{}
#'   \item{b.region.vec}{}
#'   \item{ppts}{}
#'   \item{model.obj}{}
#'   \item{result.temp}{}
#'   \item{cv.rawdata}{}
#'   \item{X.pred}{}
#'   \item{X}{}
#'   \item{parms}{}
#'   \item{cv.var}{}
#'   \item{cv.pred}{}
#'   \item{i}{}
#'   \item{cv.groups}{}
#'   \item{obs}{}
#'   \item{mons}{}
#'   \item{cv.mx}{}
#'   \item{west.monitors}{}
#'   \item{east.monitors}{}
#'   \item{weco.monitors}{}
#'   \item{rawdata}{Rawdata as passed into function}
#'   \item{desc.vars}{desc.vars as passed into function}
#'   \item{pls.comps}{number of PLS components as passed into function}
#'   \item{UK.varnames}{Universal kriging variables as passed into function}
#'   \item{manual.cv}{pre-determined cross-validation groups as passed into function}
#'   \item{factr}{?? as passed into function}
#'   \item{verbose}{as passed into function}
#'   \item{my.segments}{}
#'   \item{regional}{logical indicator, 1 = regional, 0 = national}
#'   \item{kriging}{logical indication, 1 = kriging, 0 = LUR}
#'   \item{clean.data}{logical indicator, 1 = do data cleaning, 0 = data cleaning already done}
#' }
#' @keywords 
#' @export
#' @examples 

PLSK.predict <- function(model, rawdata, desc.vars, debug.mode=FALSE){
  if(!(c("location_id") %in% names(rawdata))){ stop("location_id must be in rawdata")}
		
   # Number of PLS components
   model$pls.comps <- as.integer(model$pls.comps) 

   # Time
   data.time       <- system.time({	
     rawdata       <- as.data.table(rawdata)
     setkey(rawdata,location_id)
	
     # Find transformed names and back-transform them
     transformednames.keep <- union(
				setdiff(names(model$rawdata), 
                                names(model$excludevars)),
				c("latitude","longitude","native_id","location_id"))
     originalnames.keep <- backtransform_llnames(backtransform_lognames(transformednames.keep))

     # Drop excluded variables
     remove.vars <- setdiff(names(rawdata), originalnames.keep)
     if(length(remove.vars >0)){
		rawdata[, eval(remove.vars) := NULL, with=TRUE]
     }

     # Remove duplicate locations
     n.dups    <- sum(duplicated(rawdata))
     if(n.dups > 0){
       cat("removing ", n.dups, " duplicated entries (duplicated on location_id)\n")
	 rawdata <- unique(rawdata)		
     }
  
     ## Data processing ##

     rawdata <- log_transform_distances_dt(rawdata, desc.vars=desc.vars)
     rawdata <- combine_a23_ll_dt(rawdata, desc.vars=desc.vars)

     modeled.vars <- rownames(model$model.obj$PLS$coefficients)
     ..missing    <- as.logical(rowSums(is.na( rawdata[,modeled.vars,with=FALSE] )))
     cat("dropping ", sum(..missing), " rows in the rawdata (cohort) data due to missingness in geocovariates\n")
     rawdata      <- rawdata[!..missing,]	

     rawdata      <- as.data.frame(rawdata)
 
     stopifnot(
	 length(setdiff(colnames(model$model.obj$covars), modeled.vars)) +
	 length(setdiff(modeled.vars, colnames(model$model.obj$covars))) ==0
     )  
     ##not sure why they wouldn't be identical...worth lookiing into if this spits out an error message
     stopifnot(all(colnames(model$model.obj$covars) %in% names(rawdata)))	
    
     ## End data processing ##

  })
      
  if(debug.mode) cat("data transformations complete. Total time:\n")
  if(debug.mode) print(data.time)

  # Time for prediction	
  pred.time <- system.time({
    all.regions <- c("east","weco","west")

    #
    # predict 1000 ppts at a time
    maxdata         <- dim(rawdata)[1]
    my.rawdata_pred <- matrix(NA, ncol=4, nrow=maxdata)
    monitors        <- rep(NA,maxdata)
    v.region.vec    <- rep(NA,maxdata)
    LUR.pred        <- rep(NA,maxdata)
    rownames(my.rawdata_pred) <- rawdata$native_id

    pls.list <- list()
    pls.i    <- 1
    for (batch in seq(1, maxdata, 1000)){
      last <- min(maxdata, batch + 999)
      ppts <- make_predict_object(rawdata[batch:last, ], model$model.obj, model$pls.comps)
      pls.list[[pls.i]] <- ppts$PLS$scores
      pls.i <- pls.i + 1
      #
      # set up regional design matrix and regional variogram vector
      b.type <- 'reg'
      v.type <- 'reg'
      b.region.vec <- model$v.hash[[b.type]](ppts)

      if(debug.mode){
        if(!all(all.regions %in% b.region.vec)){
	     cat("encountered batch containing ppts in only ", 
			sum(all.regions %in% b.region.vec), 
			" of ", length(all.regions), " regions\n"
		  )
	  }
      }

      ppts.region.vec  <- model$v.hash[[v.type]](ppts)

      model.region.vec <- model$v.hash[[v.type]](model$model.obj)
    
      ppts$X           <- model$b.hash[[b.type]](ppts, b.region.vec, model$pls.comps)
    
      regions.missing  <- all.regions[!(all.regions %in% unique(b.region.vec))]

      # Conditional expectation - i.e. prediction
      mark <- c.exp(model$model.obj$parms,
                  model$model.obj$X, model$model.obj$coords, model.region.vec, model$model.obj$y,
                  fill_X(ppts$X, regions.missing, pls.comps=model$pls.comps, uk.vars=length(model$UK.varnames)),
                  fill_coords(ppts$coords, regions.missing),
                  fill_regional_vector(ppts.region.vec, regions.missing))
      if(!identical(ppts.region.vec, b.region.vec)){
	  warn("mark defined ppts.region.vec and b.region.vec differently according to v.type versus b.type\n
	        I'm not sure if this was intentional--check into this")
      }
      #
      # take out the last 0,1, or 2 rows if you had
      # to add them in to make the kriging code work
      mark                             <- mark[1:(nrow(mark)-length(regions.missing)),]
      my.rawdata_pred[batch:last, 1:2] <- mark
      my.rawdata_pred[batch:last, 3:4] <- as.matrix(ppts$gps)
	monitors[batch:last]     <- ppts$monitors
      LUR.pred[batch:last]     <- ppts$X %*% model$model.obj$parms$beta
	v.region.vec[batch:last] <- ppts.region.vec
	if(debug.mode){cat("prediction batch ", which(batch== seq(1, maxdata, 1000)),
	 " out of ", length(seq(1, maxdata, 1000)), " complete.\n")}
    }

    my.rawdata_pred <- data.frame(my.rawdata_pred, monitors,
		v.region.vec, LUR.pred)

    names(my.rawdata_pred) <- c("pred_modelingscale","var","longitude","latitude","location_id",
	     "v.region.vec","LUR.pred")

  })
  if(debug.mode) cat("predictions complete. Total time:\n")
  if(debug.mode) print(pred.time)
  if(debug.mode){
    temp <- do.call(rbind.data.frame, pls.list)
    my.rawdata_pred <- cbind.data.frame(my.rawdata_pred, temp)
  }
  my.rawdata_pred
}

PLSK.predict.old <- function(model, rawdata,desc.vars, debug.mode=FALSE){
	if(!(c("location_id") %in% names(rawdata))){ stop("location_id must be in rawdata")}
		
	model$pls.comps <- as.integer(model$pls.comps)
data.time <- system.time({	
	rawdata <- as.data.table(rawdata)
	setkey(rawdata,location_id)
	
	transformednames.keep <- union(
				setdiff(names(model$rawdata), names(model$excludevars)),
				c("latitude","longitude","native_id","location_id")
			)
	
	originalnames.keep <- backtransform_llnames(backtransform_lognames(transformednames.keep))

	remove.vars <- setdiff(names(rawdata), originalnames.keep)
	if(length(remove.vars >0)){
		rawdata[, eval(remove.vars) := NULL, with=TRUE]
	}

	n.dups <- sum(duplicated(rawdata))
	if(n.dups > 0){
		cat("removing ", n.dups, " duplicated entries (duplicated on location_id)\n")
		rawdata <- unique(rawdata)		
	}


	
    rawdata <- log_transform_distances_dt(rawdata, desc.vars=desc.vars)
    rawdata <- combine_a23_ll_dt(rawdata, desc.vars=desc.vars)
    
   
	
    modeled.vars.pls <- rownames(model$model.obj$PLS$coefficients)
    modeled.vars.ukvars <- colnames(model$model.obj$UK.covars)
	if(!is.null(modeled.vars.ukvars)){
	    if(modeled.vars.ukvars %in% modeled.vars.pls) stop("UKcovar %in% pls vars")
	}
    modeled.vars <- c(modeled.vars.pls, modeled.vars.ukvars)

    ..missing <- as.logical(rowSums(is.na( rawdata[,modeled.vars,with=FALSE] )))
    cat("dropping ", sum(..missing), " rows in the rawdata (cohort) data due to missingness in geocovariates\n")
    rawdata <- rawdata[!..missing,]	

    rawdata <- as.data.frame(rawdata)

    stopifnot(
		length(setdiff(colnames(model$model.obj$covars), modeled.vars.pls)) +
		length(setdiff(modeled.vars.pls, colnames(model$model.obj$covars))) ==0
	)  
		##not sure why they wouldn't be identical...worth looking into if this spits out an error message

    stopifnot(all(colnames(model$model.obj$covars) %in% names(rawdata)))
    stopifnot(all( colnames(model$model.obj$UK.covars) %in% names(rawdata)))
	
})
      
	if(debug.mode) cat("data transformations complete. Total time:\n")
	if(debug.mode) print(data.time)
	

	
  pred.time <- system.time({
  all.regions <- c("east","weco","west")


  
  

  #
  # predict 1000 ppts at a time
  maxdata <- dim(rawdata)[1]
  my.rawdata_pred <- matrix(NA, ncol=4, nrow=maxdata)
  monitors <- rep(NA, maxdata)
  v.region.vec <- rep(NA,maxdata)
  LUR.pred <- rep(NA,maxdata)
  rownames(my.rawdata_pred) <- rawdata$native_id
  for (batch in seq(1, maxdata, 1000))
  {
    last <- min(maxdata, batch + 999)
    ppts <- make_predict_object(rawdata[batch:last, ], model$model.obj, model$pls.comps)

    #
    # set up regional design matrix and regional variogram vector
    b.type <- 'reg'
    v.type <- 'reg'
    b.region.vec <- model$v.hash[[b.type]](ppts)

    if(debug.mode){
       if(!all(all.regions %in% b.region.vec)){
	     cat("encountered batch containing ppts in only ", 
			sum(all.regions %in% b.region.vec), 
			" of ", length(all.regions), " regions\n"
		  )
	}
    }

    ppts.region.vec <- model$v.hash[[v.type]](ppts)

    model.region.vec <- model$v.hash[[v.type]](model$model.obj)
    
    ppts$X <- model$b.hash[[b.type]](ppts, b.region.vec, model$pls.comps)
    
    regions.missing <- all.regions[!(all.regions %in% unique(b.region.vec))]

    mark <- c.exp(likfit.obj=model$model.obj$parms,
                  X.mon=model$model.obj$X, coords.mon=model$model.obj$coords, reg.mon=model.region.vec, y=model$model.obj$y,
                  X.pred=fill_X(ppts$X, regions.missing,pls.comps=model$pls.comps,uk.vars=length(model$UK.varnames)),
                  coords.pred=fill_coords(ppts$coords, regions.missing),
                  reg.pred=fill_regional_vector(ppts.region.vec, regions.missing))
    if(!identical(ppts.region.vec, b.region.vec)){
	warn("mark defined ppts.region.vec and b.region.vec differently according to v.type versus b.type\n
	I'm not sure if this was intentional--check into this")
    }
    #
    # take out the last 0,1, or 2 rows if you had
    # to add them in to make the kriging code work
    mark <- mark[1:(nrow(mark)-length(regions.missing)),]
    my.rawdata_pred[batch:last, 1:2] <- mark
    my.rawdata_pred[batch:last, 3:4] <- as.matrix(ppts$gps)
	monitors[batch:last] <- ppts$monitors
    LUR.pred[batch:last] <- ppts$X %*% model$model.obj$parms$beta
	v.region.vec[batch:last] <- ppts.region.vec
	if(debug.mode){cat("prediction batch ", which(batch== seq(1, maxdata, 1000)),
	 " out of ", length(seq(1, maxdata, 1000)), " complete.\n")}
  }

     my.rawdata_pred <- data.frame(my.rawdata_pred, monitors,
		v.region.vec, LUR.pred)

  	names(my.rawdata_pred) <- c("pred_modelingscale","var","longitude","latitude","location_id",
		"v.region.vec","LUR.pred")

})
if(debug.mode) cat("predictions complete. Total time:\n")
if(debug.mode) print(pred.time)
my.rawdata_pred

}


PLSK.predict.old <- function(model, rawdata,desc.vars, debug.mode=FALSE){
	stop("old")
	if(!(c("location_id") %in% names(rawdata))){ stop("location_id must be in rawdata")}
		
	model$pls.comps <- as.integer(model$pls.comps)
data.time <- system.time({	
	rawdata <- as.data.table(rawdata)
	setkey(rawdata,location_id)
	
	transformednames.keep <- union(
				setdiff(names(model$rawdata), names(model$excludevars)),
				c("latitude","longitude","native_id","location_id")
			)
	
	originalnames.keep <- backtransform_llnames(backtransform_lognames(transformednames.keep))

	remove.vars <- setdiff(names(rawdata), originalnames.keep)
	if(length(remove.vars >0)){
		rawdata[, eval(remove.vars) := NULL, with=TRUE]
	}

	n.dups <- sum(duplicated(rawdata))
	if(n.dups > 0){
		cat("removing ", n.dups, " duplicated entries (duplicated on location_id)\n")
		rawdata <- unique(rawdata)		
	}


	
    rawdata <- log_transform_distances_dt(rawdata, desc.vars=desc.vars)
    rawdata <- combine_a23_ll_dt(rawdata, desc.vars=desc.vars)
    
   
	
    modeled.vars <- rownames(model$model.obj$pls$coefficients)
    ..missing <- as.logical(rowSums(is.na( rawdata[,modeled.vars,with=FALSE] )))
    cat("dropping ", sum(..missing), " rows in the rawdata (cohort) data due to missingness in geocovariates\n")
    rawdata <- rawdata[!..missing,]	

    rawdata <- as.data.frame(rawdata)

    stopifnot(
		length(setdiff(colnames(model$model.obj$covars), modeled.vars)) +
		length(setdiff(modeled.vars, colnames(model$model.obj$covars))) ==0
	)  
		##not sure why they wouldn't be identical...worth lookiing into if this spits out an error message

    stopifnot(all(colnames(model$model.obj$covars) %in% names(rawdata)))	
	
})
      
	if(debug.mode) cat("data transformations complete. Total time:\n")
	if(debug.mode) print(data.time)
	

	
  pred.time <- system.time({
  all.regions <- c("east","weco","west")


  
  

  #
  # predict 1000 ppts at a time
  maxdata <- dim(rawdata)[1]
  my.rawdata_pred <- matrix(NA, ncol=4, nrow=maxdata)
  monitors <- rep(NA, maxdata)
  v.region.vec <- rep(NA,maxdata)
  LUR.pred <- rep(NA,maxdata)
  rownames(my.rawdata_pred) <- rawdata$native_id
  for (batch in seq(1, maxdata, 1000))
  {
    last <- min(maxdata, batch + 999)
    ppts <- make_predict_object(rawdata[batch:last, ], model$model.obj, model$pls.comps)

    #
    # set up regional design matrix and regional variogram vector
    b.type <- 'reg'
    v.type <- 'reg'
    b.region.vec <- model$v.hash[[b.type]](ppts)

    if(debug.mode){
       if(!all(all.regions %in% b.region.vec)){
	     cat("encountered batch containing ppts in only ", 
			sum(all.regions %in% b.region.vec), 
			" of ", length(all.regions), " regions\n"
		  )
	}
    }

    ppts.region.vec <- model$v.hash[[v.type]](ppts)

    model.region.vec <- model$v.hash[[v.type]](model$model.obj)
    
    ppts$X <- model$b.hash[[b.type]](ppts, b.region.vec, model$pls.comps)
    
    regions.missing <- all.regions[!(all.regions %in% unique(b.region.vec))]

    mark <- c.exp(model$model.obj$parms,
                  model$model.obj$X, model$model.obj$coords, model.region.vec, model$model.obj$y,
                  fill_X(ppts$X, regions.missing,pls.comps=model$pls.comps),
                  fill_coords(ppts$coords, regions.missing),
                  fill_regional_vector(ppts.region.vec, regions.missing))
    if(!identical(ppts.region.vec, b.region.vec)){
	warn("mark defined ppts.region.vec and b.region.vec differently according to v.type versus b.type\n
	I'm not sure if this was intentional--check into this")
    }
    #
    # take out the last 0,1, or 2 rows if you had
    # to add them in to make the kriging code work
    mark <- mark[1:(nrow(mark)-length(regions.missing)),]
    my.rawdata_pred[batch:last, 1:2] <- mark
    my.rawdata_pred[batch:last, 3:4] <- as.matrix(ppts$gps)
	monitors[batch:last] <- ppts$monitors
    LUR.pred[batch:last] <- ppts$X %*% model$model.obj$parms$beta
	v.region.vec[batch:last] <- ppts.region.vec
	if(debug.mode){cat("prediction batch ", which(batch== seq(1, maxdata, 1000)),
	 " out of ", length(seq(1, maxdata, 1000)), " complete.\n")}
  }

     my.rawdata_pred <- data.frame(my.rawdata_pred, monitors,
		v.region.vec, LUR.pred)

  	names(my.rawdata_pred) <- c("pred_modelingscale","var","longitude","latitude","location_id",
		"v.region.vec","LUR.pred")

})
if(debug.mode) cat("predictions complete. Total time:\n")
if(debug.mode) print(pred.time)
my.rawdata_pred

}