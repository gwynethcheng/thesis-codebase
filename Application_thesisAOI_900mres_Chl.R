## -------------------------------------------
## Fitting the model using AECM algorithm
## -------------------------------------------

## Applying SS-FPCA to McMurdo Region Chl data
## Functions_SS_FPCA_mcmurdo.R is a slightly edited version
## to optimise for large n (valid pixels)

source("/Users/gwyneth/Desktop/Masters_Thesis/thesis-codebase/Functions_SS_FPCA.R") # using author's original code
library(splines)

## -----------------------------------------------------------
## Reads the CSV outputs and saves each variable as a
## self-contained .RData file ready for the SS-FPCA analysis
## -----------------------------------------------------------

## configuring paths
INPUT_DIR <- '/Users/gwyneth/Desktop/oct2022_mar2023_ssfpca_900m'
OUTPUT_DIR <- '/Users/gwyneth/Desktop/oct2022_mar2023_ssfpca_900m'   # can be the same folder or different
STEM <- 'ssfpca_data'
VARIABLES <- c('chl_nn', 'chl_oc4me')

for (varname in VARIABLES) {
  
  cat(sprintf('\nLoading variable: %s\n', varname))
  
  # file paths
  data_path <- file.path(INPUT_DIR, sprintf('%s_%s_data_matrix.csv', STEM, varname))
  coords_path <- file.path(INPUT_DIR, sprintf('%s_%s_pixel_coords.csv', STEM, varname))
  meta_path <- file.path(INPUT_DIR, sprintf('%s_%s_metadata.csv', STEM, varname))
  
  # read data matrix (n x N), column headers are date strings
  data_raw <- read.csv(data_path, check.names = FALSE)
  dates <- as.Date(colnames(data_raw)) # N date objects
  data.mat <- as.matrix(data_raw) # n x N numeric matrix
  storage.mode(data.mat) <- 'double'
  rownames(data.mat) <- NULL
  
  cat(sprintf('Data matrix: %d pixels x %d days\n', nrow(data.mat), ncol(data.mat)))
  cat(sprintf(' Missing: %.1f%%\n',
              100 * mean(is.na(data.mat))))
  cat(sprintf('Value range : [%.4f, %.4f]\n',
              min(data.mat, na.rm = TRUE),
              max(data.mat, na.rm = TRUE)))
  
  # read pixel coordinates
  coords <- read.csv(coords_path)
  px_x <- coords$px_x       # projected x of each pixel centre (MSLC2000)
  px_y <- coords$px_y       # projected y of each pixel centre (MSLC2000)
  row_idx <- coords$row_idx    # row index in the full H x W canvas
  col_idx <- coords$col_idx    # col index in the full H x W canvas
  
  # read grid metadata
  meta <- read.csv(meta_path)
  grid_height <- meta$grid_height
  grid_width <- meta$grid_width
  pixel_size <- meta$pixel_size_m
  epsg <- meta$epsg
  
  cat(sprintf('Canvas size : %d rows x %d cols  (%.0f m pixels, EPSG:%d)\n',
              grid_height, grid_width, pixel_size, epsg))
  
  # ocean mask: logical H x W matrix, TRUE = valid pixel
  # reconstructed from row_idx / col_idx so we can map n back to 2D
  ocean_mask <- matrix(FALSE, nrow = grid_height, ncol = grid_width)
  ocean_mask[cbind(row_idx + 1L, col_idx + 1L)] <- TRUE # +1: R is 1-indexed
  
  # save as .RData
  
  # all objects needed for the SS-FPCA analysis are bundled together.
  # data.mat   — n x N matrix  (equivalent to mcmurdo.data in the author's script)
  # dates      — length-N Date vector
  # px_x       — length-n projected x coordinates
  # px_y       — length-n projected y coordinates
  # ocean_mask — H x W logical matrix for reconstructing 2D images
  # grid_height, grid_width, pixel_size, epsg — grid metadata
  
  out_path <- file.path(OUTPUT_DIR, sprintf('%s_%s_oct2022_mar2023.RData', STEM, varname))
  save(data.mat, dates,
       px_x, px_y, row_idx, col_idx,
       ocean_mask,
       grid_height, grid_width, pixel_size, epsg,
       file = out_path)
  
  cat(sprintf('Saved : %s\n', out_path))
}

## -----------------------------------------------------------
## Run SS-FPCA on chl_oc4me, repeat for chl_nn
## -----------------------------------------------------------

load("/Users/gwyneth/Desktop/oct2022_mar2023_ssfpca_900m/ssfpca_data_chl_oc4me_oct2022_mar2023.RData")
dim(data.mat)

cat("px_x range:", range(px_x), "\n")  # expect 7,000,000+ = easting
cat("px_y range:", range(px_y), "\n")  # expect 5,000,000+ = northing

## Step 0: Preparations
# center the data by subtracting per-time-point mean before running AECM

col.means <- colMeans(data.mat, na.rm = TRUE)
mcmurdo.data <- data.mat - matrix(rep(col.means, each = nrow(data.mat)), nrow = nrow(data.mat))

# create the basis matrix using projected coordinates (px_x, px_y, in metres)
# the basis functions (B-splines) are piecewise polynomials that tile across the spatial domain
# interior knots (k1, k2) are the "join points" where one piece ends and the next begins
# more knots = more flexible, wiggly surface; fewer knots = smoother, simpler surface
# for my data, changed to 2,2 since range (px_x) and range (px_y) shows almost 1:1 aspect, with very sparse time points
# also to keep df low since fpca.ini works on columns with no. of valid values above df
k1 <- 2   # number of interior knots along longitude
k2 <- 2   # number of interior knots along latitude

knot1 <- quantile(px_x, prob=(1:k1)/(k1+1))
knot2 <- quantile(px_y, prob=(1:k2)/(k2+1))
bs.lon <- bs(px_x, knots=knot1, degree=3, intercept=T)  # longitude
bs.lat <- bs(px_y, knots=knot2, degree=3, intercept=T)  # latitude
#bs.grid <- kronecker(X=as.matrix(bs.lon), Y=as.matrix(bs.lat)) 

temp <- t(bs.lon)%*%bs.lon
R.lon <- t(chol(temp)) 
B.lon <- t(solve(R.lon, t(bs.lon)))
temp <- t(bs.lat)%*%bs.lat
R.lat <- t(chol(temp)) 
B.lat <- t(solve(R.lat, t(bs.lat)))
# Base <- kronecker(X=B.lat, Y=B.lon)    # The tensor spline basis
# edited here:
# author's Kronecker product worked because her bs.lon and bs.lat were evaluated on a grid
# lon.trim had 103 unique values and lat.trim had 45 unique values, so the Kronecker product gave (103*45) x (7*7) = 4635 x 49
# my px_x and px_y each have 5325 values (one per pixel), so the Kronecker product explodes

# fix:
# build the tensor product basis row-by-row
# for pixel i: basis row = B.lon[i,] %x% B.lat[i,]  (1x7 kron 1x7 = 1x49)
# vectorised across all n pixels
# Row-wise tensor product — explicit and shape-safe
n.px <- nrow(B.lon)   # 5325
df.lon <- ncol(B.lon) # 7
df.lat <- ncol(B.lat) # 7
df.total <- df.lon * df.lat  # 25

Base.raw <- matrix(NA, nrow=n.px, ncol=df.total)
for (i in seq_len(n.px)) {
  Base.raw[i, ] <- kronecker(B.lat[i, ], B.lon[i, ])
}

## IMPORTANT: orthonormalise the full tensor product basis.
## Individual orthonormality of B.lon and B.lat does NOT carry over
## to their row-wise Kronecker product — Base.raw must be
## orthonormalised as a whole, otherwise Lambda eigenvalues inflate
## by ~300x and H.mle diverges.
temp <- t(Base.raw) %*% Base.raw
cat("Condition number of Base.raw Gram matrix:", round(kappa(temp), 1), "\n")
R    <- t(chol(temp))
Base <- t(solve(R, t(Base.raw)))
df   <- ncol(Base)

## Verify orthonormality
BTB <- t(Base) %*% Base
cat("Diagonal range of t(Base)%*%Base (should be ~1):",
    round(range(diag(BTB)), 8), "\n")
cat("Max off-diagonal (should be ~0):",
    max(abs(BTB - diag(df))), "\n")

# some other quantities
P <- 6   # the number of functional PCs
Pmiss <- 0.9   # the Kalman filtering threshold
# only skip filtering if more than 90% of pixels are missing in a time point

n <- nrow(mcmurdo.data)
N <- ncol(mcmurdo.data)
n.obs <- sum(!is.na(mcmurdo.data))
NA.matrix <- is.na(mcmurdo.data)
# edited slightly here:
# index.lk and mask.mat no needed because my data already contains only valid pixels
# so no masking is necessary
#index.lk <- mask.mat == 1   # indices of lake pixels (mask.mat gives the lake mask)
#index.lk[is.na(index.lk)] <- FALSE

# checking time points with no valid pixels
# due to resampling into geotiffs, which were used as the base for the ss-fpca matrix
# the date columns: "what dates had at least one GeoTIFF file land in my input folder?"
# but some of those geotiffs might be empty due to resampling (900m pixels won't exist in a raster with one or two sparse 300m pixels)
n.obs.per.time <- colSums(!is.na(mcmurdo.data))
cat("Time points with zero observations:\n")
print(which(n.obs.per.time == 0))
cat("Total:", sum(n.obs.per.time == 0), "\n")

# identify and remove fully-empty time points
min.obs <- 500 # keep time points with at least roughly 10% of valid pixels observed
valid.times <- which(n.obs.per.time > min.obs)
cat("Keeping", length(valid.times), "days\n")
# RECORD: for oct 2022-mar2023, 23 days were kept

mcmurdo.data <- mcmurdo.data[, valid.times]
dates        <- dates[valid.times]        # keep dates in sync
NA.matrix    <- is.na(mcmurdo.data)

# update N
N <- ncol(mcmurdo.data)
n.obs <- sum(!is.na(mcmurdo.data))

# check what time points are left
cat("Missing rate:", round(100*mean(is.na(mcmurdo.data)), 1), "%\n")
# for oct2022-mar2023, missingness is now about 81%
cat("Obs per time point:\n")
print(summary(colSums(!is.na(mcmurdo.data))))

## Step 1: Initialization
# the state space componnt
beta0 <- mu0 <- rep(0, times=df)
P0 <- Sigma0 <-  100*diag(df)
Mt <- 1*diag(df)
sigma.h <- sum(apply(mcmurdo.data, MARGIN=1, FUN=var, na.rm=T), na.rm=T)/n
HHt <- sigma.h*diag(df)
cHt <- chol(HHt)

# the FPCA component
par.zero <- fpca.ini(data=mcmurdo.data, basis=Base, K=P, take.mean=F, pert=0.001)
Theta <- par.zero$Theta
xi <- Base %*% Theta
Lambda <- par.zero$D
sigma <- par.zero$sigma

# check Lambda is now on a sensible scale
cat("diag(Lambda):", round(diag(Lambda), 4), "\n")
# expect values roughly on the order of hundreds, not 100000

## Step 2: The AECM iterations
loglike0 <- 0
loglike1 <- 1
it <- 0

while ((abs((loglike1 - loglike0) / loglike0) > 0.0005)&&(loglike1 - loglike0 >= 0)&&(it <= 10)) {

  begin <- proc.time()

  it <- it + 1
  loglike0 <- loglike1
  
  ## added block
  # regularise H.mle to ensure positive definiteness again
  # numerical errors in the MLE can cause near-zero or negative eigenvalues
  eig.H <- eigen(HHt, symmetric=TRUE)
  if (any(eig.H$values < 1e-8)) {
    cat(sprintf("  Iteration %d: regularising HHt at top of loop (min eigenvalue = %.2e)\n",
                it, min(eig.H$values)))
    eig.H$values <- pmax(eig.H$values, 1e-8)
    HHt <- eig.H$vectors %*% diag(eig.H$values) %*% t(eig.H$vectors)
    HHt <- (t(HHt) + HHt) / 2
    cHt <- chol(HHt)
  }

  ## CYCLE 1 : the state space model
  ## Using a set of notation consistent with the KFilter functions

  ## Apply the filter and smoother (E-step)
  KS <- sparse.Psmooth.G(y=mcmurdo.data, A=Base, Phi=Mt, mu0=mu0, Sigma0=Sigma0, cHt=cHt,
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
  Plag.s <- sparse.Plag1(y=mcmurdo.data, A=Base, Phi=Mt, Pf=P.f, J=Jmat, J0=J0, K=Kgain)
  
  ## added block
  ## checkpoint in case last iteration is a decrease in loglike, which also stops the while loop
  ## returns back to previous iteration with higher loglike value
  Theta.prev  <- Theta
  Lambda.prev <- Lambda
  sigma.prev  <- sigma
  a.s.prev    <- a.s
  Mt.prev     <- Mt
  HHt.prev    <- HHt

  ## Compute the likelihood and get the MLEs (M-step)
  state.lik <- AECM.KFa.lik(Phi=Mt, KS=KS, Plag.s=Plag.s, mu0=mu0, Sigma0=Sigma0, HHt=HHt)
  a0.sum <- state.lik$a0.sum
  at.sum <- state.lik$at.sum
  S00 <- state.lik$S00
  S01 <- state.lik$S01
  S11 <- state.lik$S11

  ## The MLEs (including the scoring method)
  # edited here to indicate M = I
  Mt.mle <- diag(df)
  #Mt.temp <- diag(S01) / diag(S00)
  #Mt.mle <- diag(Mt.temp)

  Htemp <- 1/N * (S11 - 2*S01 + S00)
  #Htemp <- 1/N * (S11 - S01%*%t(Mt.mle) - Mt.mle%*%t(S01) + Mt.mle%*%S00%*%t(Mt.mle))
  #Htemp <- 1/N * (S11 - 2*S01%*%t(Mt.mle) + Mt.mle%*%S00%*%t(Mt.mle))
  H.mle <- (t(Htemp) + Htemp)/2
  
  ## added block
  # regularise H.mle to ensure positive definiteness again
  # numerical errors in the MLE can cause near-zero or negative eigenvalues
  eig.H <- eigen(H.mle, symmetric=TRUE)
  if (any(eig.H$values < 1e-8)) {
    cat(sprintf("  Iteration %d: regularising H.mle (min eigenvalue = %.2e)\n",
                it, min(eig.H$values)))
    eig.H$values <- pmax(eig.H$values, 1e-8)
    H.mle <- eig.H$vectors %*% diag(eig.H$values) %*% t(eig.H$vectors)
    H.mle <- (t(H.mle) + H.mle) / 2   # enforce symmetry
  }

  ## Update the smoothed state with the new MLEs
  Mt <- Mt.mle
  HHt <- H.mle
  cHt <- chol(HHt)

  KS.cycle1 <- sparse.Psmooth.G(y=mcmurdo.data, A=Base, Phi=Mt, mu0=mu0, Sigma0=Sigma0, cHt=cHt,
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
  Plag.cycle1 <- sparse.Plag1(y=mcmurdo.data, A=Base, Phi=Mt, Pf=P.f, J=Jmat, J0=J0, K=Kgain)

  ## CYCLE2 : the mixed model FPCA
  meanfun <- Base%*%a.s

  ## The E-step
  par.E <- AECM.fpca.Estep(data=mcmurdo.data, basis=Base, P=P, meanfun=meanfun, a.s=a.s, P.s=P.s,
                           Theta=Theta, sigma=sigma, D=Lambda)
  alpha <- par.E$alpha
  alpha2 <- par.E$alpha2
  alphabeta <- par.E$alphabeta

  ## M-step
  par.M <- AECM.fpca.Mstep(data=mcmurdo.data, basis=Base, P=P, meanfun=meanfun, a.s=a.s, P.s=P.s,
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
    index <- which(!is.na(mcmurdo.data[,i]))
    nobs[i] <- length(index)
    Z <- mcmurdo.data[index,i]
    B <- Base[index,]

    lik.obs <- - 0.5*nobs[i]*log(sigma)
    lik.alpha <- - 0.5*sum(log(diag(Lambda))) - 0.5*t(alpha[,i])%*%solve(Lambda)%*%alpha[,i]
    lik.each[i] <-  lik.obs + lik.alpha
  }

  beta.lik <- AECM.KFa.lik(Phi=Mt, KS=KS.cycle1, Plag.s=Plag.cycle1, mu0=mu0, Sigma0=Sigma0, HHt=HHt)
  lik.beta <- - 0.5*(beta.lik$at.sum + beta.lik$a0.sum)
  lik.joint <- - 0.5/sigma*rss + sum(lik.each) + lik.beta   # - 0.5*(sum(nobs) + df + P)*log(2*pi)
  loglike1 <- lik.joint
  
  ## roll back to previous iteration if likelihood decreased
  if (loglike1 < loglike0) {
    cat("  Likelihood decreased — rolling back to previous iteration\n")
    Theta  <- Theta.prev
    Lambda <- Lambda.prev
    sigma  <- sigma.prev
    a.s    <- a.s.prev
    Mt     <- Mt.prev
    HHt    <- HHt.prev
    cHt    <- chol(HHt.prev)
  }

  end <- proc.time() - begin
  print(paste('Iteration', it, '; Timer', end[1], '; loglike = ', loglike1))

}

# expected results: larger jump in log likelihood from iteration 1 to 2, followed by a smaller jump from 2 to 3, and so on
# rapid improvement early on as parameters move from their initialization toward the optimum, then progressively smaller gains as the algorithm converges
# the while loop's convergence criterion will stop when that relative change drops below 0.05%

# exporting values
# save entire workspace after AECM convergence
save.image("/Users/gwyneth/Desktop/oct2022_mar2023_ssfpca_900m/AECMresults_oct2022_mar2023_chl_oc4me.RData")

# if loading into a new session, use code below
# restores everything — mcmurdo.data, Base, df, P, N, n, dates,
# all AECM. objects, alpha, sigma, Theta, Lambda, a.s, Mt, HHt etc.

#source("/Users/gwyneth/Desktop/Masters_Thesis/thesis-codebase/Functions_SS_FPCA.R")
#library(splines)
#library(fields)
#load("/Users/gwyneth/Desktop/oct2022_mar2023_ssfpca_900m/AECM_results_chl_oc4me.RData")

## Step 3: Save the results
AECM.sigma <- sigma
AECM.Theta <- Theta

AECM.Mt <- Mt
AECM.HHt <- HHt
AECM.beta <- a.s
AECM.vbeta <- P.s
AECM.D <- Lambda

AECM.tmean <- Base %*% AECM.beta
final <- AECM.fpca.orth(data=mcmurdo.data, basis=Base, P=P, meanfun=AECM.tmean, a.s=AECM.beta,
                        P.s=AECM.vbeta, Theta=AECM.Theta, sigma=AECM.sigma, D=AECM.D)
AECM.newTheta <- final$newTheta
AECM.Lambda <- final$Lambda
AECM.alpha <- final$score

## checking the results

## orthogonality of eigenfunctions
ef      <- final$eigenfun          # n x P matrix
ef.gram <- t(ef) %*% ef
cat("t(eigenfun)%*%eigenfun (should be ~identity):\n")
print(round(ef.gram, 4))
cat("Max off-diagonal:", max(abs(ef.gram - diag(P))), "\n")
cat("Diagonal range (should be ~1):", round(range(diag(ef.gram)), 6), "\n")

## residual sum of squares
cat("\nResidual diagnostics:\n")
fitted <- AECM.tmean + Base %*% (AECM.newTheta %*% AECM.alpha)
resid  <- mcmurdo.data - fitted

cat("RSS:", round(sum(resid^2, na.rm=TRUE), 4), "\n")
cat("Mean squared residual:", round(mean(resid^2, na.rm=TRUE), 6), "\n")
cat("AECM.sigma (should be close to MSR):", round(AECM.sigma, 6), "\n")

## summary of results
dim(AECM.alpha) # should be P x N (6 x no. of days)
diag(AECM.Lambda) # the P eigenvalues in decreasing order
round(diag(AECM.Lambda) / sum(diag(AECM.Lambda)) * 100, 1) # variance proportions
diag(AECM.Mt)     # the df diagonal VAR(1) coefficients
AECM.sigma        # scalar error variance # should be small and positive

##########

# plotting Figure 6 of the paper

# just to check that results weren't edited by visualisation code
snapshot <- list(
  alpha  = AECM.alpha,
  Lambda = AECM.Lambda,
  sigma  = AECM.sigma,
  beta   = AECM.beta
)

library(fields)

col_pal <- colorRampPalette(c("#e84855", "#fffd82", "#5fad56"))(100)

# different from author's code: reconstruct 2D eigenfunction grids from irregular pixels
# final$eigenfun is n x P (one row per valid pixel, one col per PC)
# place each pixel's value back onto the full H x W canvas using row_idx and col_idx

make_ef_grid <- function(ef_vec, row_idx, col_idx, grid_height, grid_width) {
  grid <- matrix(NA, nrow=grid_height, ncol=grid_width)
  grid[cbind(row_idx + 1L, col_idx + 1L)] <- ef_vec   # +1: R is 1-indexed
  return(grid)
}

ef1_grid <- make_ef_grid(final$eigenfun[,1], row_idx, col_idx, grid_height, grid_width)
ef2_grid <- make_ef_grid(final$eigenfun[,2], row_idx, col_idx, grid_height, grid_width)
ef3_grid <- make_ef_grid(final$eigenfun[,3], row_idx, col_idx, grid_height, grid_width)

# axis coordinates
# x_coords and y_coords are pixel centre coordinates in MSLC2000 metres

# if px_x/px_y are not on a regular grid (which they may not be), 
# use col_idx and row_idx to build axis vectors instead:
x_axis <- seq(min(px_x), max(px_x), length.out=grid_width) # lon
y_axis <- seq(min(px_y), max(px_y), length.out=grid_height)  # decreasing: top=north # lat

# sign flip if needed (same logic as author)
# inspect the plots first before flipping
# ef1_grid <- -ef1_grid   # uncomment if needed to match expected orientation

# dates for x axis of score plots
# dates is your length-N Date vector (already filtered to valid.times)
date_labels <- as.numeric(format(dates, "%Y")) +
  (as.numeric(format(dates, "%j")) - 1) / 365

# output png
png("/Users/gwyneth/Desktop/oct2022_mar2023_ssfpca_900m/oct2022_mar2023_chl_oc4me_Rplot.png",
    width=16, height=9, units="in", res=300)

# plot
par(mfrow=c(2,3), mar=c(4,6,3,6))

# plot function — grid is stored row1=north, so we flip rows before transposing
# after flipping, row1=south, which matches y_axis increasing from south to north
plot_ef <- function(grid, x_axis, y_axis, title) {
  grid_flipped <- grid[nrow(grid):1, ]   # flip rows so south is row 1
  image.plot(x_axis, y_axis, t(grid_flipped),
             col=col_pal,
             xlab="Easting (m)", ylab="Northing (m)",
             main=title
             )
}

plot_ef(ef1_grid, x_axis, y_axis, "Eigenfunction 1")
plot_ef(ef2_grid, x_axis, y_axis, "Eigenfunction 2")
plot_ef(ef3_grid, x_axis, y_axis, "Eigenfunction 3")

# PC scores
plot(dates, AECM.alpha[1,], type="l", xlab="", ylab="PC scores 1",
     xaxt="n", las=1)
axis.Date(1, at=pretty(dates), format="%b %Y")
abline(h=0, lty=2, col="grey60")

plot(dates, AECM.alpha[2,], type="l", xlab="", ylab="PC scores 2",
     xaxt="n", las=1)
axis.Date(1, at=pretty(dates), format="%b %Y")
abline(h=0, lty=2, col="grey60")

plot(dates, AECM.alpha[3,], type="l", xlab="", ylab="PC scores 3",
     xaxt="n", las=1)
axis.Date(1, at=pretty(dates), format="%b %Y")
abline(h=0, lty=2, col="grey60")

dev.off()

# re-checking results
cat("alpha unchanged:", identical(snapshot$alpha,  AECM.alpha),  "\n")
cat("Lambda unchanged:", identical(snapshot$Lambda, AECM.Lambda), "\n")
cat("sigma unchanged:", identical(snapshot$sigma,  AECM.sigma),  "\n")
cat("beta unchanged:", identical(snapshot$beta,   AECM.beta),   "\n")
