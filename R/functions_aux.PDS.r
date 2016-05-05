#================================
#CREATE EXPONENTIAL COVARIANCE MATRIX
#================================
varcov.eff <- function(coords, l.pars)
{
  eye <- diag(1,nrow(coords))
  d <- as.matrix(dist(coords))
  Sig <- exp(l.pars[1])*eye+exp(l.pars[2])*exp(-d/exp(l.pars[3]))
  Sig
}


#=========================================
#Block diagonal spatial covariance matrix:
#=========================================
block.Sig <- function(l.pars, coords, reg.ind)
{
  coords <- coords[order(reg.ind),]
  reg.ind <- reg.ind[order(reg.ind)]
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  cum <- cumsum(n.vec)
  Sig <- matrix(0,nrow(coords),nrow(coords))
  Sig[1:n.vec[1],1:n.vec[1]]<- varcov.eff(coords[1:n.vec[1],],l.pars[1:3])
  if (length(n.vec)>1) {
    for (j in 2:length(cum)) {
        Sig[(cum[j-1]+1):(cum[j]),(cum[j-1]+1):(cum[j])] <- varcov.eff(coords[(cum[j-1]+1):(cum[j]),],l.pars[(3*j-2):(3*j)])
    }
  }
  Sig.inv <- solve(Sig)
  list('Sig.inv'=Sig.inv,'Sig'=Sig)
}


#====================
#PROFILE LIKELIHOOD
#====================
prof.lik <- function(l.pars, X, coords, reg.ind, data)
{
  m.mat <- model.matrix(~ -1 + as.matrix(X))
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  Sigs <- block.Sig(l.pars,coords,reg.ind)
  Sig <- Sigs$Sig; Sig.inv <- Sigs$Sig.inv
  ###Infer beta from the covariance matrix and evaluate likelihood:
  beta <- solve(crossprod(m.mat,SpatioTemporal::blockMult(Sig.inv,m.mat,n.blocks=length(n.vec),block.sizes=n.vec)),
                crossprod(m.mat,SpatioTemporal::blockMult(Sig.inv,data,n.blocks=length(n.vec),block.sizes=n.vec)))
#                crossprod(m.mat,blockMult(Sig.inv,data,n.x=1,n.blocks=length(n.vec),block.sizes=n.vec)))
  lik <- mvtnorm::dmvnorm(data,m.mat%*%beta,Sig,log=T)
  -lik
}

#==============================
#FUNCTION TO INFER BETA FROM 
#MAXIMIZED COV PARS
#==============================
inf.beta<- function(l.pars, X, coords, reg.ind, data)
{
  # to be sure, reorder by reg.ind as done elsewhere for block operations
  X <- X[order(reg.ind),]
  coords <- coords[order(reg.ind),]
  data <- data[order(reg.ind)]
  reg.ind <- reg.ind[order(reg.ind)]
  #

  #m.mat <- model.matrix(~as.matrix(X))
  m.mat <- X   # X should be a matrix; and the "model.matrix" function
  # would require us to be certain whether to add an intercept.
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  Sig.inv <- block.Sig(l.pars,coords,reg.ind)$Sig.inv
  beta <- solve(crossprod(m.mat,SpatioTemporal::blockMult(Sig.inv,m.mat,n.blocks=length(n.vec),block.sizes=n.vec)),
                crossprod(m.mat,SpatioTemporal::blockMult(Sig.inv,data,n.blocks=length(n.vec),block.sizes=n.vec)))
#                crossprod(m.mat,blockMult(Sig.inv,data,n.x=1,n.blocks=length(n.vec),block.sizes=n.vec)))
  beta
}

