#======================================================================================
# LIKELIHOOD FUNCTION FOR REGIONAL VARIOGRAMS
# 
# ini.l.pars: Length 3xR vector of initial (log) starting values for the covariance
# parameters, where "R" is the number of distinct regions; order of starting values 
# should be e.g. c(tau1,sigma1,rho1,tau2,sigma2,rho2,tau3,sigma3,rho3)
# 
# X: Matrix of kriging covariates (in our case, PLS scores)
# coords: Nx2 matrix of coordinates
# reg.ind: N-vector of distinct kriging regions.  For traditional UK, set reg.ind = rep(1,N)
# data: Monitor pollutant exposures
# hess: Logical; whether or not to compute hessian and test whether it is positive-definite.
# A positive-definite hessian is a good indicator that true log-likelihood max has been found
# 
# Returns: List containing maximized log-covariance parameters; kriging regression covariates;
# maximized log-likelihood; indicator of positive definite hessian, if hess=TRUE
#======================================================================================
#======================================================================================
my.likfit <- function(ini.l.pars, X, coords, reg.ind, data, hess=TRUE, trace=0,factr=1e7)
{
  X <- X[order(reg.ind),]
  coords <- coords[order(reg.ind),]
  data <- data[order(reg.ind)]
  reg.ind <- reg.ind[order(reg.ind)]
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  lb <- rep(c(log(0.0001),log(0.0001),log(0.001)),length(n.vec))
  ub <- rep(c(log(10),log(10),log(5000)),length(n.vec))
  opt <- optim(ini.l.pars, prof.lik,coords=coords,
               data=data, X=X, reg.ind= reg.ind, method="L-BFGS-B",
               lower=lb,upper=ub,hessian=hess,control=list(trace=trace,factr=factr))
  # temporary debugging
  print("optim done")
  assign("opt",opt,pos=".GlobalEnv")
  #
  beta <- inf.beta(l.pars = opt$par, X = X, coords = coords, 
                   reg.ind = reg.ind, data  = data)
  if (hess) out <- list('log.cov.pars'=opt$par,'beta'=beta,'max.log.lik'=-opt$value,'hess.pd'=corpcor::is.positive.definite(opt$hessian), 'hessian'=opt$hessian)
  if (!hess) out <- list('log.cov.pars'=opt$par,'beta'=beta,'max.log.lik'=-opt$value)
  out
}


#Fake example:
#us.regions <- c('east','weco','east','east','west','weco','east')
#get.pars <- my.likfit(ini.l.pars = rep(c(0,log(0.5),log(500)),3)
#				X = national$pls$scores[,1:2],
#				coords = national$coords,
#				reg.ind = us.regions)




#===============================
# AUXILIARY FUNCTION FOR c.exp()
#===============================
k.pred <- function(reg.dat,dim.X,beta,log.cov.pars)
{
  #browser()
  ind <- as.vector(by(reg.dat,reg.dat$p.ind,nrow))
  nm <- ind[1]
  nt <- sum(ind)
  y <- reg.dat$y.all[which(reg.dat$p.ind==1)]
 
  mu <- as.matrix(reg.dat[1:nm,3:(2+dim.X)])%*%beta
  vcov <- varcov.eff(reg.dat[,1:2],log.cov.pars)
  vcov12 <- vcov[1:nm,(nm+1):nt]
  vcov11 <- vcov[1:nm,1:nm]
  vcov22 <- vcov[(nm+1):nt,(nm+1):nt]
  pred.var <- vcov22-crossprod(vcov12,solve(vcov11,vcov12))
  reg.dat$y.all[which(reg.dat$p.ind==2)] <- as.matrix(reg.dat[(nm+1):nt,3:(2+dim.X)])%*%beta+t(vcov12)%*%solve(vcov11)%*%(y-mu)
  reg.dat$v.all[which(reg.dat$p.ind==2)] <- diag(pred.var)
  reg.dat[,c('mp.ind','y.all','v.all')]
}





#=======================================================================
# FUNCTION TO COMPUTE CONDITIONAL EXPECTATION (I.E., MAKE PREDICTIONS)
#
# likfit.obj: Object output from my.likfit()
# X.mon: Covariates (i.e., PLS scores) for monitors
# coords.mon: Monitor coordinates
# reg.mon: Vector of regions for monitors
# y: Monitor pollutant concentrations
# X.pred: Covariates (i.e., PLS scores) for prediction locations
# coords.pred: Prediction coordinates
# reg.pred: Vector of regions for prediction sites
#
# Returns: N x 2 matrix, where N is the number of prediction locations;
# 1st column contains predictions, 2nd column contains prediction variances
#=======================================================================
#=======================================================================

c.exp <- function(likfit.obj,
                  X.mon, coords.mon, reg.mon, y,
                  X.pred, coords.pred, reg.pred)
{
  log.cov.pars <- likfit.obj$log.cov.pars
  beta <- likfit.obj$beta
  n1 <- nrow(coords.mon)
  n2 <- nrow(coords.pred)
  p.ind <- c(rep(1,n1),rep(2,n2))
  #X.all <- model.matrix(~rbind(X.mon,X.pred))
  X.all <- rbind(X.mon,X.pred)
  dim.X <- ncol(X.all)
  coords.all <- rbind(coords.mon,coords.pred)
  reg.all <- c(reg.mon,reg.pred)
  y.all <- c(y,rep(NA,n2))
  v.all <- rep(NA,n1+n2)
  pred.id <- 1:n2
  mp.ind <- c(rep(0,n1),pred.id)
  
  dat.all <- cbind(coords.all,X.all,reg.all,mp.ind,p.ind,y.all,v.all)
  c.o <- dat.all[order(reg.all),]
  uq.regs <- unique(dat.all$reg.all)
  n.regs <- length(uq.regs)
   regional.new <- NULL
  #browser()
  for (j in 1:length(uq.regs))
  {
      regional.old <- regional.new
      regional.dat <- c.o[which(c.o$reg.all==uq.regs[j]),]
	 
  ### NOTE: Change made by Casey Olives on 7-11-2013
      #beta.reg <- beta[(3*j-2):(3*j)]
	if(length(unique(regional.dat$p.ind)) > 1){
      	regional.new <- rbind(regional.old,
			k.pred(regional.dat,dim.X,beta,log.cov.pars[(3*j-2):(3*j)]))
	}
  }
  
  predictions <- regional.new[match(pred.id,regional.new$mp.ind),]
  out <- cbind('predictions'=predictions$y.all,'pred.var'=predictions$v.all)
  out

}




