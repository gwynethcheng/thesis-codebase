## -------------------------------------------
## Fitting the model using AECM algorithm
## -------------------------------------------

## An example of applying SS-FPCA to Taruo Lake Chlorophyll data

source("/Users/gwyneth/Desktop/Masters_Thesis/Codebase/SSFPCA/Functions_SS_FPCA.R")

load("/Users/gwyneth/Desktop/Masters_Thesis/Codebase/SSFPCA/Taro_Chl_trimmed.RData")

## Step 0: Preparations
# Center the data
Taro.mean <- apply(Taro.trim[,,-1], MARGIN=3, FUN=mean, na.rm=T)
Taro.data <- apply(Taro.trim[,,-1], MARGIN=3, FUN=as.vector)
Taro.data <- Taro.data - matrix(rep(Taro.mean, each=nrow(Taro.data)), nrow=nrow(Taro.data), ncol=100)

# create the basis matrix
k1 <- 4   # number of interior knots along longitude
k2 <- 2   # number of interior knots along latitude

knot1 <- quantile(lon.trim, prob=(1:k1)/(k1+1))
knot2 <- quantile(lat.trim, prob=(1:k2)/(k2+1))
bs.lon <- bs(lon.trim, knots=knot1, degree=3, intercept=T)  # longitude
bs.lat <- bs(lat.trim, knots=knot2, degree=3, intercept=T)  # latitude
bs.grid <- kronecker(X=as.matrix(bs.lon), Y=as.matrix(bs.lat)) 

temp <- t(bs.lon)%*%bs.lon
R.lon <- t(chol(temp)) 
B.lon <- t(solve(R.lon, t(bs.lon)))
temp <- t(bs.lat)%*%bs.lat
R.lat <- t(chol(temp)) 
B.lat <- t(solve(R.lat, t(bs.lat)))
Base <- kronecker(X=B.lat, Y=B.lon)    # The tensor spline basis
df <- ncol(Base)   # the basis dimension

# some other quantities
P <- 6   # the number of functional PCs
Pmiss <- 0.9   # the Kalman filtering threshold

n <- nrow(Taro.data)
N <- ncol(Taro.data)
n.obs <- sum(!is.na(Taro.data))
NA.matrix <- is.na(Taro.data)
index.lk <- mask.mat == 1   # indices of lake pixels (mask.mat gives the lake mask)
index.lk[is.na(index.lk)] <- FALSE

## Step 1: Initialization
# the state space componnt
beta0 <- mu0 <- rep(0, times=df)
P0 <- Sigma0 <-  100*diag(df)
Mt <- 1*diag(df)
sigma.h <- sum(apply(Taro.data, MARGIN=1, FUN=var, na.rm=T), na.rm=T)/n
HHt <- sigma.h*diag(df)
cHt <- chol(HHt)

# the FPCA component
par.zero <- fpca.ini(data=Taro.data, basis=Base, K=P, take.mean=F, pert=0.001)
Theta <- par.zero$Theta
xi <- Base %*% Theta
Lambda <- par.zero$D
sigma <- par.zero$sigma

## Step 2: The AECM iterations
loglike0 <- 0
loglike1 <- 1
it <- 0

while ((abs((loglike1 - loglike0) / loglike0) > 0.0005)&&(loglike1 - loglike0 >= 0)&&(it <= 10)) {

  begin <- proc.time()

  it <- it + 1
  loglike0 <- loglike1

  ## CYCLE 1 : the state space model
  ## Using a set of notation consistent with the KFilter functions

  ## Apply the filter and smoother (E-step)
  KS <- sparse.Psmooth.G(y=Taro.data, A=Base, Phi=Mt, mu0=mu0, Sigma0=Sigma0, cHt=cHt,
                         Theta=Theta, Lambda=Lambda, sigma=sigma,
                         Pmiss=Pmiss, NA.mat=NA.matrix)
  a.f <- KS$xf
  P.f <- KS$Pf
  a.s <- KS$xs
  P.s <- KS$Ps
  a0 <- KS$x0n
  P0 <- KS$P0n
  Jmat <- KS$J
  J0 <- KS$J0
  Kgain <- KS$K
  Plag.s <- sparse.Plag1(y=Taro.data, A=Base, Phi=Mt, Pf=P.f, J=Jmat, J0=J0, K=Kgain)

  ## Compute the likelihood and get the MLEs (M-step)
  state.lik <- AECM.KFa.lik(Phi=Mt, KS=KS, Plag.s=Plag.s, mu0=mu0, Sigma0=Sigma0, HHt=HHt)
  a0.sum <- state.lik$a0.sum
  at.sum <- state.lik$at.sum
  S00 <- state.lik$S00
  S01 <- state.lik$S01
  S11 <- state.lik$S11

  ## The MLEs (including the scoring method)
  Mt.temp <- diag(S01) / diag(S00)
  Mt.mle <- diag(Mt.temp)

  Htemp <- 1/N * (S11 - 2*S01%*%t(Mt.mle) + Mt.mle%*%S00%*%t(Mt.mle))
  H.mle <- (t(Htemp) + Htemp)/2

  ## Update the smoothed state with the new MLEs
  Mt <- Mt.mle
  HHt <- H.mle
  cHt <- chol(HHt)

  KS.cycle1 <- sparse.Psmooth.G(y=Taro.data, A=Base, Phi=Mt, mu0=mu0, Sigma0=Sigma0, cHt=cHt,
                                Theta=Theta, Lambda=Lambda, sigma=sigma,
                                Pmiss=Pmiss, NA.mat=NA.matrix)
  a.f <- KS.cycle1$xf
  P.f <- KS.cycle1$Pf
  a.s <- KS.cycle1$xs
  P.s <- KS.cycle1$Ps
  a0 <- KS.cycle1$x0n
  P0 <- KS.cycle1$P0n
  a.p <- KS.cycle1$xp
  P.p <- KS.cycle1$Pp
  Jmat <- KS.cycle1$J
  J0 <- KS.cycle1$J0
  Kgain <- KS.cycle1$K
  Plag.cycle1 <- sparse.Plag1(y=Taro.data, A=Base, Phi=Mt, Pf=P.f, J=Jmat, J0=J0, K=Kgain)

  ## CYCLE2 : the mixed model FPCA
  meanfun <- Base%*%a.s

  ## The E-step
  par.E <- AECM.fpca.Estep(data=Taro.data, basis=Base, P=P, meanfun=meanfun, a.s=a.s, P.s=P.s,
                           Theta=Theta, sigma=sigma, D=Lambda)
  alpha <- par.E$alpha
  alpha2 <- par.E$alpha2
  alphabeta <- par.E$alphabeta

  ## M-step
  par.M <- AECM.fpca.Mstep(data=Taro.data, basis=Base, P=P, meanfun=meanfun, a.s=a.s, P.s=P.s,
                           Theta=Theta, alpha=alpha, alpha2=alpha2, alphabeta=alphabeta)
  Lambda.mle <- par.M$D
  Theta.mle <- par.M$Theta
  sigma.mle <- par.M$sigma
  rss <- par.M$rss

  # update the parameters and related quantities
  Theta <- Theta.mle
  Lambda <- Lambda.mle
  sigma <- sigma.mle

  ## Update the joint log-likelihood
  lik.each <- 1:N
  nobs <- 1:N
  for (i in 1:N){
    index <- which(!is.na(Taro.data[,i]))
    nobs[i] <- length(index)
    Z <- Taro.data[index,i]
    B <- Base[index,]

    lik.obs <- - 0.5*nobs[i]*log(sigma)
    lik.alpha <- - 0.5*sum(log(diag(Lambda))) - 0.5*t(alpha[,i])%*%solve(Lambda)%*%alpha[,i]
    lik.each[i] <-  lik.obs + lik.alpha
  }

  beta.lik <- AECM.KFa.lik(Phi=Mt, KS=KS.cycle1, Plag.s=Plag.cycle1, mu0=mu0, Sigma0=Sigma0, HHt=HHt)
  lik.beta <- - 0.5*(beta.lik$at.sum + beta.lik$a0.sum)
  lik.joint <- - 0.5/sigma*rss + sum(lik.each) + lik.beta   # - 0.5*(sum(nobs) + df + P)*log(2*pi)
  loglike1 <- lik.joint

  end <- proc.time() - begin
  print(paste('Iteration', it, '; Timer', end[1], '; loglike = ', loglike1))

}

# expected results: larger jump in log likelihood from iteration 1 to 2, followed by a smaller jump from 2 to 3, and so on
# rapid improvement early on as parameters move from their initialization toward the optimum, then progressively smaller gains as the algorithm converges
# the while loop's convergence criterion will stop when that relative change drops below 0.05%

# actual results:
# [1] "Iteration 1 ; Timer 1372.844 ; loglike =  269031.711211574"
# [1] "Iteration 2 ; Timer 1295.106 ; loglike =  274408.24708853"
# [1] "Iteration 3 ; Timer 1259.248 ; loglike =  275256.579083301"
# [1] "Iteration 4 ; Timer 1300.835 ; loglike =  275566.59969419"
# [1] "Iteration 5 ; Timer 1294.189 ; loglike =  275756.518880994"
# [1] "Iteration 6 ; Timer 1250.63 ; loglike =  275896.30731276"
# [1] "Iteration 7 ; Timer 1252.221 ; loglike =  276002.862922499"
# total time: about 2.5 hours
# under this criterion, the paper indicates that the AECM algorithm converged after 6 iterations

## Step 3: Save the results
AECM.sigma <- sigma
AECM.Theta <- Theta

AECM.Mt <- Mt
AECM.HHt <- HHt
AECM.beta <- a.s
AECM.vbeta <- P.s
AECM.D <- Lambda

AECM.tmean <- Base %*% AECM.beta
final <- AECM.fpca.orth(data=Taro.data, basis=Base, P=P, meanfun=AECM.tmean, a.s=AECM.beta,
                        P.s=AECM.vbeta, Theta=AECM.Theta, sigma=AECM.sigma, D=AECM.D)
AECM.newTheta <- final$newTheta
AECM.Lambda <- final$Lambda
AECM.alpha <- final$score

dim(AECM.alpha)   # results: 6 x 100
# AECM.alpha contains the PC scores
# 6 functional PCs (spatial patterns) x 100 monthly images
# each of the 6 rows corresponds to one functional PC, and each of the 100 columns corresponds to one time point (month)
# so entry (p,t) tells you how strongly the p-th spatial pattern (eigenfunction) was expressed at time t
# these are what get plotted as the "PC scores" time series in Figure 6 of the paper

diag(AECM.Lambda) # the 6 eigenvalues - 364.53935 259.58523 187.49505 90.73545 80.84569 43.33438
# PC1 (364.5) captures the most spatial variation and PC6 (43.3) the least
# the proportion of variation explained by each PC within the FPCA component is computed from these
# e.g. PC1 explains 364.5 / (364.5 + 259.6 + ... + 43.3) of the FPCA component's variance
# author's eigenvalues: 356.18 254.21 and 170.95

diag(AECM.Mt)     # the 48 diagonal VAR(1) coefficients
# [1]  0.99071988  0.25676050  0.91664143  0.80556774  0.42166228  0.61007698  0.65916999
# [8]  0.40343003  0.05680648  0.13047491  0.07195282  0.33603752  0.11706907  0.03838722
# [15]  0.39293071  0.34362607  0.03227352  0.15557162 -0.13955674  0.44942758  0.26087274
# [22]  0.22913813 -0.13051509  0.26323127  0.22119567  0.20722971  0.13972063  0.28259556
# [29]  0.11686617  0.27138875  0.28628643  0.08558775  0.38137290  0.11028403  0.20422169
# [36]  0.70556868  0.21891737  0.22276455  0.10804944  0.28067328  0.93134967  0.83816534
# [43]  0.67878649  0.38540827  0.28943592  0.54206124  0.83412551  0.78128318

# the model first compresses each monthly image into a much smaller set of numbers - 48 "basis coefficients"
# think of these 48 basis coefficients as a mathematical summary of the spatial shape of that month's image
# each value asks "how much does this particular basis coefficient this month look like it did last month?"
# each coefficient evolves as its own process - this month's value = m_k x last month's value + small random nudge
# values closer to 1 suggest stable or persistent spatial features from month to month
# e.g. basis coefficients 1, 41, 47, 48 (values 0.83–0.99) are highly persistent, suggesting that whatever spatial variation those basis functions encode tends to evolve slowly and smoothly over time 
# values closer to 0 suggest spatial features that are essentially noise with no temporal memory

AECM.sigma        # scalar error variance - 0.090862
# this is the leftover noise
# the part of each pixel's log-Chl value that the model couldn't explain with either the spatial patterns or the temporal dynamics
# since the data is log-transformed chlorophyll, 0.091 is quite small, suggesting the model is fitting the data reasonably well and most of the meaningful variation has been captured

##########

# plotting Figure 6 of the paper

library(fields)

col_pal <- colorRampPalette(c("darkgreen", "#57C200", "#E7E604", "#EEB58F", "#EFDBD9"))(200)

par(mfrow=c(2,3), mar=c(4,4,3,5))

# eigenfunctions
# final$eigenfun is a matrix with 4635 rows (one per lake pixel) and 6 columns (one per PC)
# column 1 is the first eigenfunction — a single value per pixel describing the spatial pattern of PC1
# this line reshapes that flat vector of 4635 values into a 103×45 grid (longitude × latitude) so it can be plotted as a 2D map
# same is for ef2 and ef3
ef1_mat <- matrix(final$eigenfun[,1], nrow=nlon, ncol=nlat)
ef2_mat <- matrix(final$eigenfun[,2], nrow=nlon, ncol=nlat)
ef3_mat <- matrix(final$eigenfun[,3], nrow=nlon, ncol=nlat)

# multiplying an eigenfunction by −1 gives an equally valid solution mathematically
# this flips the sign so that the spatial pattern matches the orientation shown in the paper
# the corresponding PC scores are also flipped (-AECM.alpha[1,]) to keep the two consistent with each other
ef1_mat <- -ef1_mat

# mask.mat is a 103×45 binary matrix where 1 = lake pixel, NA = land/border
# reshape mask.mat into the same grid format as the eigenfunctions
lake_mask <- matrix(as.vector(mask.mat), nrow=nlon, ncol=nlat)
# set any pixel that isn't a lake pixel to NA - appears blank in the plot
ef1_mat[lake_mask != 1 | is.na(lake_mask)] <- NA
ef2_mat[lake_mask != 1 | is.na(lake_mask)] <- NA
ef3_mat[lake_mask != 1 | is.na(lake_mask)] <- NA

# drawing the 2D colour maps
image.plot(lon_unique, lat_unique, ef1_mat,
           col=col_pal,
           xlab="Longitude", ylab="Latitude", main="Eigenfunction 1")

image.plot(lon_unique, lat_unique, ef2_mat,
           col=col_pal,
           xlab="Longitude", ylab="Latitude", main="Eigenfunction 2")

image.plot(lon_unique, lat_unique, ef3_mat,
           col=col_pal,
           xlab="Longitude", ylab="Latitude", main="Eigenfunction 3")

# line plot for PC scores
# AECM.alpha is a 6×100 matrix of PC scores — one score per PC per time point
# row 1 gives the scores for PC1 across all 100 months
plot(years_used, -AECM.alpha[1,], type="l", xlab="", ylab="PC scores 1")
abline(h=0, lty=2, col="grey60")
plot(years_used, AECM.alpha[2,], type="l", xlab="", ylab="PC scores 2")
abline(h=0, lty=2, col="grey60")
plot(years_used, AECM.alpha[3,], type="l", xlab="", ylab="PC scores 3")
abline(h=0, lty=2, col="grey60")
