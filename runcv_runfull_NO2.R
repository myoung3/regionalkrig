#setwd("C:/Users/micha/Documents/git/regionalkrig")
setwd("C:\\LocalDataAmanda\\git_repo\\regionalkrig")
#library(data.table)

#install.packages(c("SpatioTemporal","corpcor","data.table","geoR","maps","mvtnorm","pls","deldir","ggplot2","mapproj", "devtools","roxygen2"))


#

# regional definitions east.polygon and weco.polygon
attach('R/regional_polygons.RData', pos=2)
assign('weco.polygon', weco.polygon, envir=.GlobalEnv)
assign('east.polygon', east.polygon, envir=.GlobalEnv)


library(devtools)
load_all()

setwd("C:\\LocalDataAmanda\\git_repo")
load_all(pkg = "ruk")

setwd("C:\\LocalDataAmanda\\git_repo\\regionalkrig")



DR0119 <- init.dr(reqnum = "DR0119", poll_list="PM25") #generate a function that reads in DR0119
DR0119 <- init.dr(reqnum = "dr0250", poll_list="NO2", dataloc = "Q:/eac_database/requests/") #generate a function that reads in DR0119
#why not just read in DR0119? This is more efficient:
#it reads in the dataset once, stores the data in the environment of the DR0119 function
#then DR0119 alows you to select subsets of that dataset without rereading it each time
###this seems complicated but it's more efficient when I'm reading in the same data several times (memoization)

#note that r/runmodel_UKvars.R  contains the wrapper functions to actually run the model: runcv and runfull


###establish data parameters
randomsample <- .3 #this is just for my data-read function. it can take a random subset of the dataset. NULL is the full dataset or a value between 0 and 1. .5 is a 50% sample
desc.vars <- c('county','state','state_plane', 'lambert_x','lambert_y')
LID_drop=c("318638" , "318637", "2287")  ###these locations have been problematic (ie, missing covariates) in the past
vardrop <- c("m_to_truck", "m_to_oil", "m_to_6oil", "m_to_main_cityhall", "m_to_local_cityhall", "satellite_NO2","no2_behr_2005","no2_behr_2007")    
#remove these variables from the dataset if you have them. we don't model with them
#m_to_truck has had distance problems in some states, the oil and cityhall variables are not available everythwer, satellite_NO2 is the variable given to us by marshall et al

pollutant <- "NO2" #this is just an input to my dataread function

yearparam <- 2010:2014 #input to my dataread function

rawdata.l <- lapply(yearparam, DR0119, pollutant=pollutant, vardrop=vardrop, LID_drop=LID_drop,randomsample=NULL)  #read in the data using the data-read function DR0119
##the code for this function is in functions_readdata.r

#check for missing rows
stopifnot(!any(unlist(lapply(rawdata.l, apply,1,function(.x){any(is.na(.x))}))))

##identify location_id of missing rows (IF THERE ARE MISSINGS)
lapply(rawdata.l, function(.rawdata){
  missingrow <- apply(.rawdata, 1, function(.x) any(is.na(.x)))
  .rawdata$location_id[missingrow]
})

#you can specify manual cv groups with run.cv(..., manual.cv=<integer vector>)


###########Run CV models

  ###non-parallel version since you're just running one year:
  
  NO2test.cv <- lapply(rawdata.l, function(.rawdata){
    runcv(pls=2, rawdata=.rawdata, manual.cv=NULL, desc.vars=desc.vars, UK.varnames=NULL)
  })
  
#######run full models

NO2test.full <- lapply(rawdata.l, function(.rawdata){
  runfull(pls=2, rawdata=.rawdata, desc.vars=desc.vars, UK.varnames=NULL)
})

###########Run One CV
test1.2010 <- runfull(pls=2, rawdata=rawdata.l[[1]], desc.vars=desc.vars, UK.varnames=NULL)

