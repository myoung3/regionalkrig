
###library calls are in this function so it can be run as a cluster function
###... passed to DRfunction includes arguments like LID_drop, vardrop and future arguments for DR functions
runcv <- function(pls, UK.varnames, rawdata=NULL,desc.vars, factr=1e10, manual.cv=NULL,...){

	if(is.null(rawdata)){	stop("DRfunction depreciated.  specify rawdata")}
	PLSK.cv(rawdata=rawdata,desc.vars, pls.comps=pls,UK.varnames=UK.varnames, factr=factr,verbose=TRUE,manual.cv=manual.cv, ...)
}



runfull <- function(pls, UK.varnames, rawdata=NULL, desc.vars, factr=1e10,randomsample=NULL, ...){

	if(is.null(rawdata)){	stop("integrated DRfunction depreciated. specify rawdata")}

	PLSK.full(rawdata=rawdata,desc.vars, pls.comps=pls,UK.varnames=UK.varnames, factr=factr,verbose=TRUE, ...)
}
