##Code for the variance correction paper
### Source codes accompanying the paper ``Revisiting the random shift approach for testing in spatial statistics''
### by Tomáš Mrkvička, Jiří Dvořák, Jonatan Gonzáles and Jorge Mateu,
### submitted to Spatial Statistics.

### This file provides an example of the testing procedures for 
### testing independence of a pair of random fields.

### Version 5.2.2020


### Required external packages ###
##################################

library(spatstat)
library(geoR)


### Simulated data ###
######################

# test.points ... the set of sampling locations, here represented as a point pattern
#                 (the default observation window is the unit square)
# covariateA .... one of the random fields, available in the whole observation window,
#                 here represented as a pixel image
# covariateB .... one of the random fields, can be observed only at sampling points,
#                 here represented as a pixel image but only values at test.points are used

auxA <- attr(rLGCP("exp", mu=0, var=1, scale=0.2, saveLambda=TRUE),"Lambda")
auxB <- attr(rLGCP("exp", mu=0, var=1, scale=0.2, saveLambda=TRUE),"Lambda")
covariateA <- eval.im(log(auxA))
covariateB <- eval.im(log(auxB))
test.points <- runifpoint(100)


### Random shift with torus correction ###
##########################################

# Number of random shifts to be used in the test
N.perm <- 999

observed.covariance <- cov(covariateA[test.points],covariateB[test.points])
values.simulated <- rep(NA, times=N.perm)

for (k in 1:N.perm){
  test.points.shift <- rshift(test.points, edge="torus", radius=0.5)
  values.simulated[k] <- cov(covariateA[test.points.shift],covariateB[test.points])
}

test.rank <- rank(c(observed.covariance,values.simulated))[1]
p.value.torus <- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
p.value.torus


### Random shift with minus correction ###
##########################################

# Number of random shifts to be used in the test
N.perm <- 999

W.central <- owin(xrange=c(1/3,2/3),yrange=c(1/3,2/3))
test.points.central <- test.points[W.central]
valuesA.central <- covariateA[test.points.central]

observed.covariance <- cov(valuesA.central,covariateB[test.points.central])
values.simulated <- rep(NA, times=N.perm)

for (k in 1:N.perm){
  points.shift <- shift.ppp(test.points.central, vec=runif(2, min=-1/3, max=1/3))
  valuesB.shift <- covariateB[points.shift]
  values.simulated[k] <- cov(valuesA.central,valuesB.shift)
}
test.rank <- rank(c(observed.covariance,values.simulated))[1]
p.value.minus <- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
p.value.minus


### Random shift with variance correction (count) ###
#####################################################

# Number of random shifts to be used in the test
N.perm <- 999

values.simulated <- rep(NA, times=N.perm+1)
values.simulated[1] <- cov(covariateA[test.points],covariateB[test.points])
n.simulated <- rep(NA, times=N.perm+1)
n.simulated[1] <- test.points$n

# Random shifts
for (k in 1:N.perm){
  jump <- runifdisc(1, radius = 0.5)
  test.points.shifted <- shift(test.points, c(jump$x,jump$y))
  W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
  test.points.reduced <- test.points[W.reduced]
  test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
  covA <- covariateA[test.points.reduced.backshifted]
  covB <- covariateB[test.points.reduced]
  values.simulated[k+1] <- cov(covA,covB)
  n.simulated[k+1] <- test.points.reduced$n
}

# Variance correction (count)
values.std <- (values.simulated - mean(values.simulated))*sqrt(n.simulated)
test.rank <- rank(values.std)[1]
p.value.count <- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
p.value.count


### Random shift with variance correction (kernel) ###
######################################################

# Number of random shifts to be used in the test
N.perm <- 999

# Value of bandwidth to be used in the kernel regression
bw <- 0.10

# Epanechnikov kernel used for determining the weights in the regression
epanechnikov <- function(x,bw){
  return(pmax(0,(1/bw)*3/4*(1-(x/bw)^2)))
}

values.simulated <- rep(NA, times=N.perm+1)
values.simulated[1] <- cov(covariateA[test.points],covariateB[test.points])
n.simulated <- rep(NA, times=N.perm+1)
n.simulated[1] <- test.points$n
shift.vecs <- matrix(data=NA, nrow=N.perm+1, ncol=2)
shift.vecs[1,] <- c(0,0)

# Random shifts
for (k in 1:N.perm){
  jump <- runifdisc(1, radius = 0.5)
  shift.vecs[k+1,] <- c(jump$x,jump$y)
  test.points.shifted <- shift(test.points, c(jump$x,jump$y))
  W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
  test.points.reduced <- test.points[W.reduced]
  test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
  covA <- covariateA[test.points.reduced.backshifted]
  covB <- covariateB[test.points.reduced]
  values.simulated[k+1] <- cov(covA,covB)
  n.simulated[k+1] <- test.points.reduced$n
}

# Standardization using kernel regression (locally constant regression),
# compute matrix of pairwise distances between the shift vectors
diff_ij <- function(i,j) sqrt(rowSums((shift.vecs[i,]-shift.vecs[j,])^2))
dist.mat <- outer(seq_len(N.perm+1), seq_len(N.perm+1), diff_ij)

# Variance correction (kernel)
values.simulated.std2 <- (values.simulated-mean(values.simulated))^2
mvar <- rep(NA, times=N.perm+1)
for (i in 1:(N.perm+1)){
  wis <- epanechnikov(dist.mat[i,], bw=bw) # weights in the regression
  mvar[i] <- sum(wis*values.simulated.std2)/sum(wis)
}

values.std <- (values.simulated - mean(values.simulated))/sqrt(mvar)
test.rank <- rank(values.std)[1]
p.value.kernel <- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
p.value.kernel



### Random shift with variance correction (exact variance) ###
##############################################################

# Number of random shifts to be used in the test
N.perm <- 999

values.simulated <- rep(NA, times=N.perm+1)
values.simulated[1] <- cov(covariateA[test.points],covariateB[test.points])
var.simulated <- rep(NA, times=N.perm+1)
n <- test.points$n
dist.mat <- pairdist(test.points)

# Estimate variogram model for covariateA
geo.A <- data.frame(x=test.points$x, y=test.points$y, z=covariateA[test.points])
geo.A <- as.geodata(geo.A)
v.A <- variog(geo.A, uvec=seq(0,1,by=0.05), messages=FALSE)
obj.A <- function(par){# par[1] = sigma^2, par[2] = s
  sum( (v.A$v - par[1]*(1-exp(-v.A$u/par[2])) )^2)
}
est.par.A <- optim(par=c(1,1), fn=obj.A)$par
cov.mat.A <- est.par.A[1]*exp(-dist.mat/est.par.A[2])

# Estimate variogram model for covariateB
geo.B <- data.frame(x=test.points$x, y=test.points$y, z=covariateB[test.points])
geo.B <- as.geodata(geo.B)
v.B <- variog(geo.B, uvec=seq(0,1,by=0.05), messages=FALSE)
obj.B <- function(par){
  # par[1] = sigma^2, par[2] = s
  sum( (v.B$v - par[1]*(1-exp(-v.B$u/par[2])) )^2)
}
est.par.B <- optim(par=c(1,1), fn=obj.B)$par
cov.mat.B <- est.par.B[1]*exp(-dist.mat/est.par.B[2])

# Variance of the test statistics for the observed data
r.sums <- rowSums(cov.mat.A) %o% rep(1,times=n)
c.sums <- rep(1,times=n) %o% colSums(cov.mat.A)
rc.sum <- sum(rowSums(cov.mat.A))
A <- cov.mat.A + rc.sum/n^2 - r.sums/n - c.sums/n
r.sums <- rowSums(cov.mat.B) %o% rep(1,times=n)
c.sums <- rep(1,times=n) %o% colSums(cov.mat.B)
rc.sum <- sum(rowSums(cov.mat.B))
B <- cov.mat.B + rc.sum/n^2 - r.sums/n - c.sums/n
var.simulated[1] <- sum(rowSums(A*B))/(n-1)^2

# Random shifts
for (k in 1:N.perm){
  jump <- runifdisc(1, radius = 0.5)
  test.points.shifted <- shift(test.points, c(jump$x,jump$y))
  W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
  test.points.reduced <- test.points[W.reduced]
  test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
  covA <- covariateA[test.points.reduced.backshifted]
  covB <- covariateB[test.points.reduced]
  values.simulated[k+1] <- cov(covA,covB)
  
  # Which points are not discarded by the random shift
  keep <- inside.owin(x=test.points$x, y=test.points$y, w=W.reduced)
  n.reduced <- sum(keep)
  
  # Variance of the test statistics for the shifted data
  cov.mat.A.red <- cov.mat.A[keep,keep]
  cov.mat.B.red <- cov.mat.B[keep,keep]
  r.sums <- rowSums(cov.mat.A.red) %o% rep(1,times=n.reduced)
  c.sums <- rep(1,times=n.reduced) %o% colSums(cov.mat.A.red)
  rc.sum <- sum(rowSums(cov.mat.A.red))
  A <- cov.mat.A.red + rc.sum/n.reduced^2 - r.sums/n.reduced - c.sums/n.reduced
  r.sums <- rowSums(cov.mat.B.red) %o% rep(1,times=n.reduced)
  c.sums <- rep(1,times=n.reduced) %o% colSums(cov.mat.B.red)
  rc.sum <- sum(rowSums(cov.mat.B.red))
  B <- cov.mat.B.red + rc.sum/n.reduced^2 - r.sums/n.reduced - c.sums/n.reduced
  var.simulated[k+1] <- sum(rowSums(A*B))/(n.reduced-1)^2
}

# Variance corrections
values.std <- (values.simulated - mean(values.simulated))/sqrt(var.simulated)
test.rank <- rank(values.std)[1]
p.value.var <- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
p.value.var

### Source codes accompanying the paper ``Revisiting the random shift approach for testing in spatial statistics''
### by Tomáš Mrkvička, Jiří Dvořák, Jonatan Gonzáles and Jorge Mateu,
### submitted to Spatial Statistics.

### This file provides an example of the testing procedures for 
### testing independence in a bivariate point process, based on the cross K-function.

### Version 6.2.2020


### Required external packages ###
##################################

library(spatstat)
library(GET)


### Auxiliary function ###
##########################

K.trans <- function(Z, rr = FALSE){
  # Estimate cross K-function with translation edge-correction
  # Z .... marked point pattern with marks "X" and "Y" (bivariate point pattern)
  # rr ... logical, should also the vector of argument values be returned?
  
  K.hat <- Kcross(Z, "X", "Y", r = seq(0, 0.15, length.out = 50), correction = "translation")
  if (rr) return(list(r = K.hat$r, vals = K.hat$trans)) else return(K.hat$trans)
}


### Simulated data ###
######################

# Z ... bivariate point pattern to be tested, with marks "X" and "Y"

X <- rLGCP("exp", mu = 4.5, var = 1, scale = 0.2, saveLambda = FALSE)
Y <- rLGCP("exp", mu = 4.5, var = 1, scale = 0.2, saveLambda = FALSE)
Z <- superimpose(X = X, Y = Y)


### Random shift with torus correction and cross K-function ###
###############################################################

# Number of random shifts to be used in the test
N.perm <- 999

K.observed <- K.trans(Z, rr = TRUE)
r <- K.observed$r
K.observed <- K.observed$vals

RS.torus <- function() {
  Z.shift <- rshift.ppp(Z, which='X', radius=0.5, edge="torus")
  Kreduced <- K.trans(Z.shift)
  return(Kreduced)
}

K.simulated <- replicate(N.perm, RS.torus())

# Global envelope test
CS <- create_curve_set(list(r = r, obs = K.observed, sim_m = K.simulated))
p.value.torus.K <- attr(rank_envelope(CS, type = "erl"), "p")
p.value.torus.K


### Random shift with minus correction and cross K-function ###
###############################################################

# Number of random shifts to be used in the test
N.perm <- 999

W.central <- owin(xrange=c(1/3,2/3),yrange=c(1/3,2/3))
K.observed <- K.trans(Z[W.central], rr = TRUE)
r <- K.observed$r
K.observed <- K.observed$vals
X <- split(Z)[[1]]
Y <- split(Z)[[2]]

RS.minus <- function() {
  sh.vec <- runif(n=2, min=-1/3, max=1/3)
  X.shift <- as.ppp(ppp(x=X$x+sh.vec[1], y=X$y+sh.vec[2], window=W.central))
  Z.shift.red <- superimpose(X = X.shift, Y = Y[W.central])
  if (sum(Z.shift.red$marks=="X")==0){return(NA)}
  if (sum(Z.shift.red$marks=="Y")==0){return(NA)}
  Kreduced <- K.trans(Z.shift.red)
  return(Kreduced)
}

# Suppress warning messages caused by rejecting shifted points that do not fit into W.central
options(warn = -1)
K.simulated <- replicate(N.perm, RS.minus())
options(warn = 0)

# Global envelope test
CS <- create_curve_set(list(r = r, obs = K.observed, sim_m = K.simulated))
p.value.minus.K <- attr(rank_envelope(CS, type = "erl"), "p")
p.value.minus.K


### Random shift with variance correction (kernel) and cross K-function ###
###########################################################################

# Number of random shifts to be used in the test
N.perm <- 999

# Value of bandwidth to be used in the kernel regression
bw <- 0.10

# Epanechnikov kernel used for determining the weights in the regression
epanechnikov <- function(x,bw){
  return(pmax(0,(1/bw)*3/4*(1-(x/bw)^2)))
}

# Random shift function, alse returning the displacement vector
Rshift <- function(Z, radius = 0.5, which = "X"){
  Y <- split(Z, marks(Z))
  id <- seq_along(Y)
  names(id) <- names(Y)
  iwhich <- id[which]
  nwhich <- id[id != iwhich]
  jump <- runifdisc(1, radius=radius)
  Xn <- shift(Y[[iwhich]], jump)
  Wf <- intersect.owin(Y[[iwhich]]$window, Xn$window)
  Xok <- inside.owin(Xn$x, Xn$y, Wf)
  Yok <- inside.owin(Y[[nwhich]]$x, Y[[nwhich]]$y, Wf)
  XN <- ppp(x=Xn$x[Xok], y=Xn$y[Xok], window=Wf)
  YN <- ppp(x=Y[[nwhich]]$x[Yok], Y[[nwhich]]$y[Yok], window=Wf)
  return(list(Z=superimpose(X=XN, Y=YN), vi=c(jump$x, jump$y)))
}

K.observed <- K.trans(Z, rr = TRUE)
r <- K.observed$r
K.observed <- K.observed$vals

RS.kernel <- function() {
  Zreduc.shift <- Rshift(Z, radius = 0.5)
  Z.shift <- Zreduc.shift$Z
  Kreduced <- K.trans(Z.shift)
  vi <- Zreduc.shift$vi
  return(list(Kr = Kreduced, vi = vi))
}

simu <- replicate(N.perm, RS.kernel())
K.simulated <- sapply(simu[1, ], "[")
vsimu <- sapply(simu[2, ], "[")

K.simulated <- cbind(K.simulated, K.observed) # Observed K-functions
Kmean <- apply(K.simulated, 1 , mean)         # Sampling mean
Si <- sweep(K.simulated, 1, Kmean, FUN = "-") # Raw differences

# Kernel smoothing
vsimu <- cbind(vsimu, c(0,0))
aij <- pairdist.default(t(vsimu))
KK <- epanechnikov(aij,bw=bw) # weights for the kernel regression
KK <- matrix(KK,nrow=N.perm+1)
wij <- KK/sum(KK)
si2 <- (Si ^ 2) %*% wij
S <- (si2 ^ (-1/2)) * (Si)
S[is.nan(S)] <- 0

# Global envelope test
CS <- create_curve_set(list(r = r, obs = S[ , N.perm + 1], sim_m = S[ , 1:N.perm]))
p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
p.value.kernel.K


