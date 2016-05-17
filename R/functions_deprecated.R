combine_a23_ll <- function(my.data)
{
  a2.vars <- grep("ll_a2", colnames(my.data))
  for (i in a2.vars)
  {
    newcol.index <- 1 + length(colnames(my.data))
    a3.var <- grep(gsub('a2','a3',colnames(my.data)[i]), colnames(my.data))
    my.data[, newcol.index] <- my.data[, i] + my.data[, a3.var]
    colnames(my.data)[newcol.index] <- paste("ll_a23_", strsplit(colnames(my.data)[a3.var], '_')[[1]][3], sep="")
  }
  ll.vars <- grep("ll_a[^1]_s", colnames(my.data))
  my.data <- my.data[, -ll.vars]
  return(my.data)
}

log_transform_distances <- function(my.data)
{
  distance.vars <- grep("^m_to_", colnames(my.data))
  new.varnames <- c()
  for (i in distance.vars)
  {
    newcol.index <- 1 + length(colnames(my.data))
    # throw out records with covariate values of -1
    my.data <- my.data[!my.data[, i] == -1, ]
    my.data[, newcol.index] <- log10(sapply(my.data[, i], function(Col) { max(10, Col)}))
    colnames(my.data)[newcol.index] <- paste('log10_', colnames(my.data)[i], sep='')
  }
  my.data <- my.data[, -distance.vars]
  return (my.data)
}

fill_X_old <- function(X, regions.missing, all.regions=c('east','weco','west'), pls.comps)
{
  if (!length(regions.missing))
    return (X)
  for (i in match(regions.missing, all.regions))
  {
    #X <- rbind(X, rep(0, length(all.regions)*pls.comps))
    X <- rbind(X, rep(0, length(all.regions)*(pls.comps+1)))
    #X[nrow(X), (i*pls.comps-(pls.comps-1)):(i*pls.comps)] = rep(.2, pls.comps)
    X[nrow(X), (i*(pls.comps+1)-pls.comps):(i*(pls.comps+1))] = c(1, rep(.2, pls.comps))
  }
  return (X)
}