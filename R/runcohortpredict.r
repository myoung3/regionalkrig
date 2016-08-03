
#' Predictf
#'
#' A function that 
#' @param model
#' @param desc.vars
#' @param LUR
#' @param cohort
#' @param load.functions
#' @return NA, produces plot
#' @keywords 
#' @export
#' @examples 

predictf <- function(model, desc.vars, LUR, cohort, load.functions=FALSE){
	if(load.functions){
		require(data.table)
		require(pls)
		library(geoR)
		library(maps)
		library(mvtnorm)
		library(pls)
		library(SpatioTemporal)
		library(corpcor)
		source('H:/exposuremodels/national PLSKrig model/functions_readdata.r')
		source('H:/exposuremodels/national PLSKrig model/PLSK_source.R')
	}

	pred        <- as.data.table(PLSK.predict(model=model, cohort, desc.vars=desc.vars, debug.mode=TRUE))
	predname    <- paste(model$pollutant, "pred",model$year,sep='_')
	LURpredname <- paste(model$pollutant, "LURpred",model$year,sep='_')
	setkey(pred, location_id)
	pred[, eval(predname) := pred_modelingscale^2]

	if(LUR==TRUE){
		pred[, eval(LURpredname) := LUR.pred^2]
	}
	
	pred[, c("pred_modelingscale","LUR.pred", "var","longitude","latitude") := NULL]
	if(is.null(model$`...counter<RESERVED>`)) pred[, v.region.vec := NULL]  #keep region vector only for first model
	pred
}


#' PLS Prediction from Multiple Models
#'
#' A function that 
#' @param allmodels
#' @param cohort
#' @param desc.vars
#' @param LUR
#' @param parallel
#' @param ncores
#' @param onemodel
#' @return NA, produces plot
#' @keywords 
#' @export
#' @examples 


PLSK.batchpredict <- function(allmodels, cohort,desc.vars= c('county','state','state_plane', 'lambert_x','lambert_y'), LUR=FALSE,
	                        parallel=FALSE, ncores=2, onemodel=FALSE){

  if(all(c("rawdata","model.obj") %in% names(allmodels))&!onemodel) stop("it appears that allmodels is only one model. use option onemodel=TRUE")
  if(parallel & onemodel) stop("parallel can't be used with onemodel")

  allmodels[[1]]$`...counter<RESERVED>` <- 1 #keep region vector only for first model

  if(onemodel){
	out  <- predictf(allmodels, desc.vars=desc.vars, LUR=LUR, cohort=cohort)
	setkey(out, location_id)
	return(out)
  } 

  if(parallel){
	require(parallel)
	cl   <- makeCluster(ncores)
	#clusterExport(cl = cl, setdiff(ls(),c("allmodels","cohort")))
	temp <- parLapplyLB(cl=cl,fun=predictf,X=allmodels, desc.vars=desc.vars,
		LUR=LUR,cohort=cohort,load.functions=TRUE)
	stopCluster(cl)
  }else{

   	temp <- lapply(allmodels, predictf, desc.vars=desc.vars, LUR=LUR, cohort=cohort)
  }

  if(length(allmodels) > 1L){
  	out <- temp[[1]][temp[[2]]]
  }else{
	out <- temp[[1]]
  }

  if(length(allmodels) >= 3L){
	for(i in 3:length(temp)){
	  out <- out[temp[[i]]]
	}
  }


  setkey(out, location_id)
  out
}



