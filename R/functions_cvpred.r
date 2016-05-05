
#=======================#
# RETURNS THE ESTIMATES #
#=======================#
cv.pred <- function(national, us.regions, design.matrix, my.pars, fit.monitors, cv.monitors, my.pls.comps)
{
  n.regions <- length(unique(us.regions))
  all.regions <- sort(unique(us.regions))
  fit.design.matrix <- matrix(design.matrix[as.character(fit.monitors), ], ncol=ncol(design.matrix))
  fit.coords <- national$coords[as.character(fit.monitors), ]
  fit.regions <- us.regions[as.character(fit.monitors)]
  fit.y <- national$y[as.character(fit.monitors)]
  pred.design.matrix <- matrix(design.matrix[as.character(cv.monitors), ], ncol=ncol(design.matrix))
  pred.coords <- national$coords[as.character(cv.monitors), ]
  pred.regions <- us.regions[as.character(cv.monitors)]
  mark <- c.exp(my.pars,
                fit.design.matrix, fit.coords, fit.regions, fit.y,
                pred.design.matrix, pred.coords, pred.regions)
  return(mark)
}

