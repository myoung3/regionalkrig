
Output

ini.l.pars Vector of ?
region.vec List of vector of named strings, indicating which monitors are in which regions
b.region.vec List of vector of named strings, indicating which monitors are in which regions
v.type Character indicating ?.  Options include "reg"
b.type Character indicating ?.  Options include "reg"
v.hash List of 2 objects: reg:function and nat:function
b.hash List of 2 objects: reg:function and nat:function
model.obj List of 13 (see below)
exclude.vars List of string vector listing covariates
all.vars  List of string vector listing covariates
miss List of Boolean vector indicating ?
year List of integer indicating year of model
pollutant List of string indicating pollutant
invars List of string vector, listing all variables passed to modeling process
rawdata List of data.frame inluding all original data
desc.vars List of string vector, listing ? 
pls.comps List of integer indicating number of PLS components
UK.varnames List of ?
factr List of numeric ?
verbose List of boolean; indicates?
regional List of boolean; if true then regional, if false then single region

Modeling Object

y Vector of named numbers (concentrations, I assume)
UK.covars ?
covars.pls A matrix of something ?
PLS List of 19 (see below)
obs An integer indicating the number of observations (I think)
coords A data.frame with x and y coordinates (lambert projection)
gps A data.frame with lat and long coordinates
monitors String vector listing monitors
weco.monitors String vector listing monitors in the weco region
east.monitors String vector listing monitors in the east region
west.monitors String vector listing monitors in the west region
X A matrix with rows equal to the number of monitors and columns?
parms List of 5 (see below)

PLS

coefficients List of 3-dimensional matrix (number of... monitors plus 32?, y, and 5 components?)
scores Matrix that is... monitors by 5 components
loadings Matrix that is number of monitors plus 32 by 5
loading.weights a similar matrix?
Yscores matrix that is number of monitors by 5
Yloadings 1 by 5 matrix
projection 366 by 5 matrix
Xmeans numeric vector
Ymeans Number indicating? 
fitted.values 334  by 5 matrix
residuals 334x1x5 matrix
Xvar Vector of 5 named numbers 
Xtotvar Number indicating?
fit.time Named number indicating?
ncomp Integer.... 5.
method "kernelpls"
call  language plsr(formula = y ~ covars.pls, ncomp = pls.comps, data = model.obj, validation = "none")
terms Classes 'terms', 'formula' length 3 y ~ covars.pls
model data frame of y and covars.pls, which is 334 by 366?

Parms

log.cov.pars List of vector of named numbers; includes tau, sigma and rho for each region (variogram parms?)
beta Matrix (instead of vector?) intercept and betas (b1 and b2?) for each region
max.log.lik A list of a number
hess.pd A list of a boolean indicating whether the hessian was positive definite (?)
hessian A matrix (9x9) which is the Hessian.
