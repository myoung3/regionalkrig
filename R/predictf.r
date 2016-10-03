
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


