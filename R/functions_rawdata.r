
#============================#
# get raw data for model fit #
# from the MESA database     #
#============================#

#' Pulling Data from mysql Database
#'
#' Function for pulling data from mysql.  WARNING: has anyone used this in the last 3 years?!
#' @param pollutant Character string indicating pollutant to pull (options: PM25, O3, NOx, NO2)
#' @param year An integer indicating the year of model, e.g. 2013
#' @keywords 
#' @export
#' @examples 

rawdata_modeling_database <- function(pollutant, year)
{
  channel <- odbcConnect("ozone", uid="markr9", case="tolower")
  tbl = paste(pollutant,'_y',year,sep='')
  a = sqlQuery(channel, "show tables in mark;")
  if (!sum(apply(a, 1, function(x) (return (x==tbl)))))
  {
    procedure.query = paste('call mark.agency_yearly("',pollutant,'", "',year,'-01-01");', sep='')
    sqlQuery(channel, procedure.query)
    wait <- sqlQuery(channel, "show tables;")
  }
  fetch.query = paste('select * from mark.',pollutant,'_y',year,';', sep='')
  pm <- sqlQuery(channel, fetch.query)
  odbcClose(channel)
  return (pm)
}


#=============================#
# get raw data for prediction #
# from the MESA database      #
#=============================#

#' Pulling Data from mysql Database
#'
#' Function for pulling data from mysql.  WARNING: has anyone used this in the last 3 years?!
#' @param cohort Character string indicating locations to pull, e.g. mesa, sisters
#' @param sp A charachter string indicating the state plane, e.g. ALzoneE
#' @keywords 
#' @export
#' @examples 

rawdata_predict_database <- function(cohort, sp)
{
  channel <- odbcConnect("ozone", uid="markr9", case="tolower")
  tbl = paste(cohort,'_',sp,'_covars',sep='')
  a = sqlQuery(channel, "show tables in mark;")
  if (!sum(apply(a, 1, function(x) (return (x==tbl)))))
  {
    procedure.query = paste('call mark.ppt_geocovariates("',cohort,'", "',sp,'");', sep='')
    sqlQuery(channel, procedure.query)
    wait <- sqlQuery(channel, "show tables;")
  }
  fetch.query = paste('select * from mark.',tbl,';', sep='')
  pm <- sqlQuery(channel, fetch.query)
  odbcClose(channel)
  return (pm)
}


#
# remove a table from the database
drop_rawdata_tbl <- function(cohort, sp)
{
  channel <- odbcConnect("ozone", uid="markr9", case="tolower")
  tbl = paste(cohort,'_',sp,'_covars',sep='')
  procedure.query = paste('drop table if exists mark.',tbl,';', sep='')
  retval = sqlQuery(channel, procedure.query)
  wait <- sqlQuery(channel, "show tables;")
}


#==============================#
# get raw data from a flatfile #
#==============================#

rawdata_modeling_flatfile <- function(pollutant, year)
{
  fname <- paste(pollutant,'_',year,'.txt', sep='')
  pm <- read.csv(pmfile)
  return (pm)
}


