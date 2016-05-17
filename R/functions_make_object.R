

#' Make prediction object
#'
#' Make prediction object
#' @param rawdata See rawdata under PLSK.full
#' @param model.obj
#' @param pls.comps Integer indicating the number of PLS components
#' @details
#' @return A list
#' \describe{
#'   \item{covars.pls}{Stuff}
#'   \item{PLS$scores}{}
#'   \item{UK.covars}{}
#'   \item{monitors}{}
#'   \item{obs}{}
#'   \item{coords}{}
#'   \item{gps}{}
#'   \item{weco.monitors}{}
#'   \item{east.monitors}{}
#'   \item{west.monitors}{}
#' }
#' @keywords 
#' @export
#' @examples 


#unfinished--consider the functions that USE this object (ie, the process of prediction)
##don't forget it's model.obj$PLS not model.obj$pls now
make_predict_object <- function(rawdata, model.obj, pls.comps)
{
  ppts <- c()
  ppts$covars.pls <- scale(rawdata[, colnames(model.obj$covars.pls),drop=FALSE],
                       center=attr(model.obj$covars.pls, 'scaled:center'),
                       scale=attr(model.obj$covars.pls, 'scaled:scale'))
  ppts$PLS$scores <- predict(model.obj$PLS, newdata=ppts, comps=1:pls.comps, type='scores')  
  ppts$UK.covars <- as.matrix(rawdata[, colnames(model.obj$UK.covars),drop=FALSE])
  ppts$monitors <- rawdata$location_id
  ppts$obs <- length(ppts$monitors)
  ppts$coords <- rawdata[, c('lambert_x','lambert_y')]/1000
  ppts$gps <- rawdata[, c('longitude','latitude')]
  ppts$weco.monitors <- ppts$monitors[poly.int2(weco.polygon, ppts$gps)]
  ppts$east.monitors <- ppts$monitors[poly.int2(east.polygon, ppts$gps)]
  if (length(c(ppts$weco.monitors, ppts$east.monitors)) == 0)
    ppts$west.monitors <- ppts$monitors
  else
    if (length(ppts$weco.monitors) == 0)
      ppts$west.monitors <- ppts$monitors[-c(match(ppts$east.monitors, ppts$monitors))]
    else if (length(ppts$east.monitors) == 0)
      ppts$west.monitors <- ppts$monitors[-c(match(ppts$weco.monitors, ppts$monitors))]
    else
      ppts$west.monitors <- ppts$monitors[-match(unlist(list(ppts$east.monitors, ppts$weco.monitors)), ppts$monitors)]
  rownames(ppts$coords) <- ppts$monitors
  rownames(ppts$gps) <- ppts$monitors
  return(ppts)
}




#UK.varnames is a character vector corresponding to the names of variables to be included directly in the UK model (in addition to PLS components)
##new components of model.obj:
##covars.pls are the (scaled) covariates that go into the PLS models--covars except exclude and UK.vars\
##PLS is the pls output, comparable to "pls" in the original version of this funcion (renamed to spit out an error any time it comes across code that hasn't been changed)
##therefore the variables to be used in UK are cbind(1, model$PLS$scores[, 1:ncomp], model.obj$UK.covars)
##where ncomp is the number of pls.components desired


#' Make modeling object
#'
#' Make modeling object.  Executes PLS regression.
#' @param rawdata See rawdata under PLSK.full
#' @param exclude.vars String vector listing variables to exclude from PLS
#' @param UK.varnames Sting vector listing variables that will be included for universal kriging and therefore should
#'   not enter into PLS.
#' @param pls.comps Integer indicating the number of PLS components.  Default is 5.
#' @details
#' @return A list named PLS
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
#' @keywords 
#' @export
#' @examples 


make_modeling_object <- function(rawdata, exclude.vars, UK.varnames, pls.comps=5)
{
  if("data.table" %in% class(rawdata)) rawdata <- as.data.frame(rawdata)

  covars.pls.names <- setdiff(setdiff(colnames(rawdata), exclude.vars), UK.varnames)
  model.obj <- c()
  model.obj$y <- sqrt(rawdata[, 'pollutant_conc'])
  model.obj$UK.covars <- as.matrix(rawdata[, UK.varnames, drop=FALSE])
  model.obj$covars.pls <- scale(as.matrix(rawdata[,covars.pls.names,drop=FALSE]))
  model.obj$PLS <- pls::plsr(y ~ covars.pls, ncomp=pls.comps, validation='none', data=model.obj)
  model.obj$obs <- nrow(rawdata)
  model.obj$coords <- rawdata[, c('lambert_x', 'lambert_y')]/1000
  model.obj$gps <- rawdata[, c('longitude', 'latitude')]
  model.obj$monitors <- rawdata$native_id
  model.obj$weco.monitors <- model.obj$monitors[poly.int2(weco.polygon, model.obj$gps)]
  model.obj$east.monitors <- model.obj$monitors[poly.int2(east.polygon, model.obj$gps)]
  model.obj$west.monitors <- model.obj$monitors[-match(unlist(list( model.obj$east.monitors,  model.obj$weco.monitors)), model.obj$monitors)]
  rownames(model.obj$coords) <- model.obj$monitors
  rownames(model.obj$gps) <- model.obj$monitors
  names(model.obj$y) <- model.obj$monitors
  return (model.obj)
}

