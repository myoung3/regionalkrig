
#' Wrapper Function for Running Cross-Validation
#'
#' A wrapper for running cross-validation.
#' @param pls
#' @param UK.varnames A list of ?
#' @param rawdata A matrix (?) of the rawdata, which must include the following elements: ???  Defaults to NULL.  This will terminate the function call.  
#' @param desc.vars A list of ?
#' @param factr Defaults to 1e10.  What does it do?
#' @param manual.cv Defaults to NULL.  This will?
#' @keywords cross-validation
#' @export
#' @examples 

###library calls are in this function so it can be run as a cluster function
###... passed to DRfunction includes arguments like LID_drop, vardrop and future arguments for DR functions
runcv <- function(pls, UK.varnames, rawdata=NULL,desc.vars, factr=1e10, manual.cv=NULL,...){

	if(is.null(rawdata)){	stop("DRfunction deprecated.  specify rawdata")}
	PLSK.cv(rawdata=rawdata,desc.vars, pls.comps=pls,UK.varnames=UK.varnames, factr=factr,verbose=TRUE,manual.cv=manual.cv, ...)
}



#' Wrapper Function for Running a Full Model
#'
#' A wrapper for running a full model.  This function calls PLSK.full; see documentation for PLSK.full for more details.
#' @param pls An integer indicating the number of PLS components
#' @param UK.varnames A list of varnames to include in universal kriging (e.g. satellite data?)
#' @param rawdata A matrix (?) of the rawdata, which must include the following elements: ???  Defaults to NULL.  This will terminate the function call.  
#' @param desc.vars A list of variable names ?
#' @param factr Defaults to 1e10.  What does it do?
#' @param randomsample Defaults to NULL.  This will?
#' @keywords 
#' @export
#' @examples 

runfull <- function(pls, UK.varnames, rawdata=NULL, desc.vars, factr=1e10, randomsample=NULL, ...){

	if(is.null(rawdata)){	stop("integrated DRfunction deprecated. specify rawdata")}
	PLSK.full(rawdata=rawdata,desc.vars, pls.comps=pls,UK.varnames=UK.varnames, factr=factr,verbose=TRUE, ...)
}
