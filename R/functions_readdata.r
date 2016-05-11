

#' Read pollutant request
#'
#' Function for reading pollution data. 
#' @param dataloc Directory containing the data
#' @param reqnum Request number
#' @param poll_list List of pollutants to read
#' @keywords 
#' @examples
#' @return 
#' @export

init.dr <- function(dataloc="data/", reqnum, poll_list){
  stopifnot(substr(dataloc, nchar(dataloc), nchar(dataloc))=="/")
 
  # Read air pollution concentration files
  for (p in poll_list){
    temp <- fread(paste(dataloc, reqnum, "/", reqnum, "_", p, "_annual_avgs.txt",sep=''))
    temp[, native_id:=as.character(native_id)]
    assign(p, temp)
  }

  # Read covariate file
  dr.m.cov <- fread(paste0(dataloc, reqnum, "/", reqnum, "_agency_covars.txt"))
  dr.m.cov[, native_id:=as.character(native_id)]
  dr.m.cov[, location_id:=NULL] #################################remove location id from dr.m.cov because it's in the obs file
  setkey(dr.m.cov, native_id)   #################################merge by native_id
	
	
  function(year, pollutant, randomsample=NULL, LID_drop=NULL,vardrop=NULL, varkeep=NULL, remove.collocated=TRUE){
		stopifnot(sum(!is.null(vardrop), !is.null(varkeep)) <2)  #vardop and varkeep can't be used simultaneously
		if(!is.null(randomsample)){
			if(randomsample == 1){
				 stop("randomsample==1 is untested. To use all data, set randomsample=NULL")
			}
		}

            # ========================================================================================= #
            # I think eveything that follows in this box can be deleted since it was specific to DR0119
		#stopifnot( pollutant %in% c("PM25","NO2"))
		#if(pollutant == "PM25") stop("covariate database missing improve monitors in this DR")
		#if(pollutant=="PM25"){
		#	stopifnot( year %in% 1999:2012 )
		#	dr.m.obs <- PM25
		#}
		#if(pollutant=="NO2"){
		#	stopifnot( year %in% 1990:2012 )
		#	dr.m.obs <- NO2
		#}
            # ========================================================================================= #

            dr.m.obs <- get(pollutant)
		setkey(dr.m.obs,native_id)

            # Modulus operator %%
		pollutant_conc_name <- if( ((year - 1992) %% 4) == 0){
			paste("avg_", pollutant, "_366_",year,"0101",sep='')
		}else{
			paste("avg_", pollutant, "_365_", year, "0101", sep='')
		}

		###subset dr.m.obs to only native_id and year of interest
		monitor <- dr.m.obs[,c("native_id",	"location_id", pollutant_conc_name),with=FALSE]
		setnames(monitor, pollutant_conc_name,"pollutant_conc")
		if(!is.null(varkeep)){
					varkeep <- c(varkeep, "native_id", "latitude", "longitude", "lambert_x",
					"lambert_y", "state_plane", "state", "county")
					dr.m.cov <- dr.m.cov[, varkeep, with=FALSE]		
		}


		out <- dr.m.cov[monitor]
		textvars <- c("native_id","pollutant_conc","location_id",
				"latitude","longitude",	"lambert_x", "lambert_y",
				"state_plane", "state", "county")
		setcolorder(out, c(textvars, setdiff(names(out), textvars)) )
		setkey(out, latitude, longitude)
		#cat("removing ", sum(is.na(out$pollutant_conc)), " rows corresponding to missing pollutants")
		out <- out[!is.na(pollutant_conc),]
		if(!is.null(LID_drop)){
			out <- out[!location_id %in% LID_drop,]
		}
		#stopifnot(is.null(dim(duplicated(out))))
		if(remove.collocated){
			cat("removing ", sum(duplicated(out))," rows corresponding to duplicate locations")
			out <- unique(out)
		}
		if(!is.null(randomsample)){
		 	keepn <- round(nrow(out)*randomsample)
			keep <- sample(1:nrow(out),keepn, repl=FALSE)
			out <- out[keep,] 
		}
		if(out[, sum(is.na(.SD)), .SDcols=setdiff(names(out), textvars)]){ 
			missing_ids <- out[out[, apply(.SD,1, function(.v) any(is.na(.v))), .SDcols=setdiff(names(out), textvars)], location_id]
			stop("missing geocovariates in LIDs: ", paste(missing_ids, collapse=" "))
		}

		out <- as.data.frame(out)
		out <- out[,setdiff(names(out), vardrop)]
		attr(out, "year") <- year
		attr(out, "pollutant") <- pollutant
		out		
  }
}

