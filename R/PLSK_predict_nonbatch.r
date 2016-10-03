PLSK.predict.nonbatch <- function(model, rawdata,desc.vars, debug.mode=TRUE){
	if(!(c("location_id") %in% names(rawdata))){ stop("variable location_id must be in rawdata")}
	rawdata <- as.data.table(rawdata)
	setkey(rawdata,location_id)
	n.dups <- sum(duplicated(rawdata))
	if(n.dups > 0){
		cat("removing ", n.dups, " duplicated entries (duplicated on location_id)\n")
		rawdata <- unique(rawdata)		
	}

	..missing <- rowSums(is.na( rawdata[,setdiff(names(rawdata),desc.vars),with=FALSE] ))
	cat("dropping ", sum(..missing), " rows in the rawdata (cohort) data due to missingness in geocovariates\n")
	rawdata <- rawdata[..missing ==0,]	
	
	
    rawdata <- log_transform_distances_dt(rawdata,desc.vars=desc.vars)
    rawdata <- combine_a23_ll_dt(rawdata,desc.vars=desc.vars)
    rawdata <- as.data.frame(rawdata)
    stopifnot(all(colnames(model$model.obj$covars) %in% names(rawdata)))	
	
    all.regions <- c("east","weco","west")

    ppts <- make_predict_object(rawdata, model$model.obj, model$pls.comps)

    b.type <- 'reg'
    v.type <- 'reg'
    b.region.vec <- model$v.hash[[b.type]](ppts)
    if(!all(all.regions %in% b.region.vec)){
		 stop("cohortpredictno_batch requires at least one ppt in each region. \n
		 use the batch prediction function\n")
	}
    ppts.region.vec <- model$v.hash[[v.type]](ppts)
    ppts$X <- model$b.hash[[b.type]](ppts, b.region.vec, model$pls.comps)
    model.region.vec <- model$v.hash[[v.type]](model$model.obj)  

    if(!identical(ppts.region.vec, b.region.vec)){
	warn("mark defined ppts.region.vec and b.region.vec differently according to v.type versus b.type\n
	I'm not sure if this was intentional--check into this")
    }


    mark <- c.exp(likfit.obj=model$model.obj$parms,
                   X.mon=model$model.obj$X, coords.mon=model$model.obj$coords,
			reg.mon=model.region.vec, y=model$model.obj$y,
                  X.pred=ppts$X, coords.pred=ppts$coords, reg.pred=ppts.region.vec)
  my.rawdata_pred <- matrix(NA, ncol=4, nrow=nrow(rawdata))
  rownames(my.rawdata_pred) <- rawdata$native_id
  my.rawdata_pred[, 1:2] <- mark
  my.rawdata_pred[, 3:4] <- as.matrix(ppts$gps)
  my.rawdata_pred <- data.frame(my.rawdata_pred, ppts$monitors)
  names(my.rawdata_pred) <- c("pred","var","longitude","latitude","location_id")
  my.rawdata_pred
}