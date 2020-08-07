#-----------------------------------------------------------------------------------------

as.owin.ext <- function(rasExt, mask=NULL, maskValue=NA, ...) {

  # Convert an Extent object from raster package to an owin object from spatstat package.
  # Can have a mask as a raster layer (where extent(mask) == rasExt).
  #
  # USES ...
  # spatstat library version 1.45-0
  # raster library version 2.3-24
  #
  # ARGUMENTS ...
  # rasExt:    an Extent object of a raster.
  # mask:      can be a raster layer or other mask objects as specified in function owin.
  # maskValue: the value that is used to mask other data in the argument 'mask'.  That is,
  #            any cell in mask, whose value equals maskValue, will not be included in the
  #            observation window.

  # Check extent argument is valid class.
  if ( !inherits(rasExt, "Extent") ) {
    stop("Invalid class for argument 'rasExt' in function 'as.owin.ext'.")
  }

  # Check mask argument, if given.
  myMask <- NULL
  if ( !is.null(mask) ) {
    if ( inherits(mask, "RasterLayer") ) {
      # Is a raster layer, so extract the relevant values into right form for owin function.
      myMask <- as.data.frame(rasterToPoints(mask))
      if ( is.na(maskValue) ) {
        myMask[ ,3] <- !is.na(myMask[ ,3])
      } else {
        myMask[ ,3] <- myMask[ ,3] != maskValue   # Could be dodgy, test!
      }
      #myMask[,3] <- as.logical(myMask[,3])  # Converts third column back to TRUE or FALSE.
    } else {
      # Leave as is and let owin function deal with it.
      myMask <- mask
    }
  }

  # Get the observation window object.
  xrange <- c(rasExt@xmin, rasExt@xmax)
  yrange <- c(rasExt@ymin, rasExt@ymax)
  myOwin <- owin(xrange, yrange, mask=myMask, ...)

  # Return value.
  return(myOwin)
}

#-----------------------------------------------------------------------------------------

lambda.cell <- function (myFormula, coeffs, covars, lnLambda=FALSE) {

  # Intensity function.  For inhomogeneous Poisson point distribution and model of this as
  #
  #        ln(lambda(cell)) = formula with coeffs and covars(cell)
  #
  # Returns a vector containing the value of lambda(cell) for all cells in covars.
  #
  # ARGUMENTS ...
  # myFormula: a formula object that contains how to combine the information for an intensity function.
  # coeffs:    a vector containing the coefficients for the intensity formula.
  # covars:    a data.frame containing columns of covariate data with the values in each
  #            cell forming the rows.  Only contains information for cells in the domain.
  #            NB: covariate columns need to be named the same as in the formula!
  # lnLambda:  return ln(lambda(cell)) for TRUE or lambda(cell) for FALSE

  # Prepare the formula by removing any response variable.
  myTerms <- terms(myFormula, keep.order = TRUE)
  myTerms <- delete.response(myTerms)

  # Get the relevant data.
  mat <- as.matrix(model.frame(myTerms, covars))

  # Is there an intercept term?
  if ( attr(myTerms, "intercept") ) {
    ones <- rep(1, length.out=nrow(covars))
    mat <- cbind(ones, mat)
  }

  # Evaluate the function.
  lnRes <- mat %*% coeffs
  if ( lnLambda ) {
    return(drop(lnRes))
  } else {
    return(drop(exp(lnRes)))
  }

}

#-----------------------------------------------------------------------------------------

lambda.xy <- function (x, y, myFormula, coeffs, covars, lnLambda=FALSE) {

  # Intensity function.  For inhomogeneous Poisson point distribution and model of this as
  #
  #        ln(lambda(x,y)) = formula with coeffs and covars(x,y)
  #
  # Returns a vector containing the value of lambda(x,y) for all points (x,y).
  # Assumes (x,y) will have valid values in covars (i.e. only use (x,y) that do!!!)
  #
  # ARGUMENTS ...
  # x:         a vector of longitude values of the points
  # y:         a vector of latitude values of the points
  # myFormula: a formula object that contains how to combine the information for an intensity function.
  # coeffs:    a vector containing the coefficients for the intensity formula.
  # covars:    a raster stack containing numEnvirs layers and values for all x and y.
  # lnLambda:  return ln(lambda(x,y)) for TRUE or lambda(x,y) for FALSE

  # Prepare the formula by removing any response variable.
  myTerms <- terms(myFormula, keep.order = TRUE)
  myTerms <- delete.response(myTerms)

  # Get the relevant data.
  xy <- cbind(x,y)
  vals <- davesExtract.v3(covars, xy)   # NB: one column each for each raster layer in stack
  vals <- as.data.frame(cbind(xy, vals), stringsAsFactors=FALSE)
  mat <- as.matrix(model.frame(myTerms, vals))

  # Is there an intercept term?
  if ( attr(myTerms, "intercept") ) {
    ones <- rep(1, length.out=length(x))
    mat <- cbind(ones, mat)
  }

  # Evaluate the function.
  lnRes <- mat %*% coeffs
  if ( lnLambda ) {
    return(drop(lnRes))
  } else {
    return(drop(exp(lnRes)))
  }

}

#-----------------------------------------------------------------------------------------

davesExtract.v3 <- function(myRaster, xy){

  # Extract the values from myRaster at the points in xy.  Faster than raster extract function!

  # Get one layer of the raster
  if ( ! inherits(myRaster, c("RasterStack", "RasterLayer")) ) {
    stop("Unrecognised class for first argument.  Expected RasterStack or RasterLayer.")
  }

  # Get the extent of the raster (range for both x and y)
  myExt <- extent(myRaster)

  # Get the cell size for the raster
  cellSize <- res(myRaster)
  numCols <- ncol(myRaster)
  numRows <- nrow(myRaster)

  # Get the cell indices (one for row and one for column) for each xy point.
  # Indices start at top left (as for a matrix). Cells contain top border and left border.
  # Right most border becomes part of cells to the left, and, bottom most border becomes
  # part of cells above.
  icol <- floor((xy[,1] - myExt@xmin) / cellSize[1]) + 1
  icol[icol > numCols] <- numCols
  irow <- floor((myExt@ymax - xy[,2]) / cellSize[2]) + 1
  irow[irow > numRows] <- numRows
  ixy <- cbind(irow,icol)

  # Extract values from these cells.
  numLayers <- nlayers(myRaster)
  vals <- NULL
  for ( i in 1:numLayers ) {
    myMat <- as.matrix(myRaster[[i]])
    vals <- cbind(vals, myMat[ixy])
  }

  # Add column names to the return values.
  colnames(vals) <- names(myRaster)

  # Return values.
  return(vals)
}

#-----------------------------------------------------------------------------------------

fithian <- function(isP) {

  ## Family object
  # isP : vector that is TRUE for Poisson family rows and FALSE for binomial family rows.

  N <- length(isP)
  isB <- !isP
  P <- poisson(link=log)
  B <- binomial(link=cloglog)

  linkfun <- function(mu) {
    R <- double(N)
    R[isP] <- P$linkfun(mu[isP])
    R[isB] <- B$linkfun(mu[isB])
    R
  }
  linkinv <- function(eta) {
    R <- double(N)
    R[isP] <- P$linkinv(eta[isP])
    R[isB] <- B$linkinv(eta[isB])
    R
  }
  mu.eta <- function(eta) {
    R <- double(N)
    R[isP] <- P$mu.eta(eta[isP])
    R[isB] <- B$mu.eta(eta[isB])
    R
  }
  valideta <- function(eta) TRUE
  variance <- function(mu) {
    R <- double(N)
    R[isP] <- P$variance(mu[isP])
    R[isB] <- B$variance(mu[isB])
    R
  }
  validmu <- function(mu) all(is.finite(mu)) && all(mu>0 & (isP | mu<1))
  dev.resids <- function(y,mu,wt) {
    R <- double(N)
    R[isP] <- P$dev.resids(y[isP],mu[isP],wt[isP])
    R[isB] <- B$dev.resids(y[isB],mu[isB],wt[isB])
    R
  }

  aic <- function(y,n,mu,wt,dev) {
    P$aic(y[isP],n[isB],mu[isP],wt[isP],dev)+B$aic(y[isB],n[isB],mu[isB],wt[isB],dev)
  }

  ## Hard to get this right
  initialize <- substitute({
    n <- rep.int(1, nobs)
    mustart <- ifelse(isP,y+0.1,(y+0.5)/2)
  },list(isP=isP))

  structure(list(family = "fithian",
                 link = "fithian",
                 linkfun = linkfun,
                 linkinv = linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 aic = aic,
                 mu.eta = mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = valideta),
            class = "family")
}


# popa <- function(isP) {
#
#   ## Family object
#   # isP : vector that is TRUE for Poisson family rows and FALSE for binomial family rows.
#
#   N <- length(isP)
#   isB <- !isP
#   P <- poisson(link=log)
#   B <- binomial(link=cloglog)
#
#   linkfun <- function(mu) {
#     R <- double(N)
#     R[isP] <- P$linkfun(mu[isP])
#     R[isB] <- B$linkfun(mu[isB])
#     R
#   }
#   linkinv <- function(eta) {
#     R <- double(N)
#     R[isP] <- P$linkinv(eta[isP])
#     R[isB] <- B$linkinv(eta[isB])
#     R
#   }
#   mu.eta <- function(eta) {
#     R <- double(N)
#     R[isP] <- P$mu.eta(eta[isP])
#     R[isB] <- B$mu.eta(eta[isB])
#     R
#   }
#   valideta <- function(eta) TRUE
#   variance <- function(mu) {
#     R <- double(N)
#     R[isP] <- P$variance(mu[isP])
#     R[isB] <- B$variance(mu[isB])
#     R
#   }
#   validmu <- function(mu) all(is.finite(mu)) && all(mu>0 & (isP | mu<1))
#   dev.resids <- function(y,mu,wt) {
#     R <- double(N)
#     R[isP] <- P$dev.resids(y[isP],mu[isP],wt[isP])
#     R[isB] <- B$dev.resids(y[isB],mu[isB],wt[isB])
#     R
#   }
#
#   aic <- function(y,n,mu,wt,dev) {
#     P$aic(y[isP],n[isB],mu[isP],wt[isP],dev)+B$aic(y[isB],n[isB],mu[isB],wt[isB],dev)
#   }
#   ## Hard to get this right
#   initialize <- substitute({
#     n <- rep.int(1, nobs)
#     mustart <- ifelse(isP,y+0.1,(y+0.5)/2)
#   },list(isP=isP))
#
#
#   ## Extra stuff for gam
#   Pd2link <- function(mu) -1/mu^2
#   Bd2link <- function(mu) {
#     l1m <- log1p(-mu)
#     -1/((1 - mu)^2 * l1m) * (1 + 1/l1m)
#   }
#   d2link <- function(mu) {
#     R <- double(N)
#     R[isP] <- Pd2link(mu[isP])
#     R[isB] <- Bd2link(mu[isB])
#     R
#   }
#
#
#   Pd3link <- function(mu) 2/mu^3
#   Bd3link <- function(mu) {
#     l1m <- log1p(-mu)
#     mu3 <- (1 - mu)^3
#     (-2 - 3 * l1m - 2 * l1m^2)/mu3/l1m^3
#   }
#   d3link <- function(mu) {
#     R <- double(N)
#     R[isP] <- Pd3link(mu[isP])
#     R[isB] <- Bd3link(mu[isB])
#     R
#   }
#
#
#   Pd4link <- function(mu) -6/mu^4
#   Bd4link <- function(mu) {
#     l1m <- log1p(-mu)
#     mu4 <- (1 - mu)^4
#     (-12 - 11 * l1m - 6 * l1m^2 - 6/l1m)/mu4/l1m^3
#   }
#   d4link <- function(mu) {
#     R <- double(N)
#     R[isP] <- Pd4link(mu[isP])
#     R[isB] <- Bd4link(mu[isB])
#     R
#   }
#
#   Pdvar <- function(mu) rep.int(1, length(mu))
#   Bdvar <- function(mu) 1 - 2 * mu
#   dvar <- function(mu) {
#     R <- double(N)
#     R[isP] <- Pdvar(mu[isP])
#     R[isB] <- Bdvar(mu[isB])
#     R
#   }
#
#   Pd2var <- function(mu) rep.int(0, length(mu))
#   Bd2var <- function(mu) rep.int(-2, length(mu))
#   d2var <- function(mu) {
#     R <- double(N)
#     R[isP] <- Pd2var(mu[isP])
#     R[isB] <- Bd2var(mu[isB])
#     R
#   }
#
#   Pd3var <- function(mu) rep.int(0, length(mu))
#   Bd3var <- function(mu) rep.int(0, length(mu))
#   d3var <- function(mu) {
#     R <- double(N)
#     R[isP] <- Pd3var(mu[isP])
#     R[isB] <- Bd3var(mu[isB])
#     R
#   }
#
#   ls <- function(y, w, n, scale) {
#     c(-aic(y,n,y,w,0)/2,0,0)
#   }
#
#   structure(list(family = "popa",
#                  link = "popa",
#                  canonical = "popa",
#                  linkfun = linkfun,
#                  linkinv = linkinv,
#                  variance = variance,
#                  dev.resids = dev.resids,
#                  aic = aic,
#                  mu.eta = mu.eta,
#                  initialize = initialize,
#                  validmu = validmu,
#                  valideta = valideta,
#                  d2link=d2link,
#                  d3link=d3link,
#                  d4link=d4link,
#                  dvar=dvar,
#                  d2var=d2var,
#                  d3var=d3var,
#                  ls=ls),
#             class = "family")
# }

#-----------------------------------------------------------------------------------------

bglm <- function(A, B, Y, weights=1, offset=0, family, control=glm.control()) {

  # Arguments ...
  # A:       species specific partition of design matrix (nrow = num samples,
  #          ncol = 1 + num species specific coefficients)
  # B:       non species specific partition of design matrix (nrow = num samples,
  #          ncol = num non species specific coefficients)
  # Y:       response matrix (nrow = num samples, ncol = num species)
  # weights: scalar or matrix (nrow = num Samples, ncol = num species)
  # offset:  offset vector (same for all species, length = num samples)

  invperm <- function(p) {
    p[p] <- seq_along(p)
    p
  }

  # Dimensions ...
  nobs <- nrow(Y)
  mY <- ncol(Y)
  mA <- ncol(A)
  mB <- ncol(B)
  mbeta <- mY*mA+mB

  if ( nrow(A) != nobs )
    stop("Number of rows in A not equal to number of observations in Y.")
  if ( nrow(B) != nobs )
    stop("Number of rows in B not equal to number of observations in Y.")
  if ( length(offset) > 1 && length(offset) != nobs)
    stop("Length of offset vector not equal to number of observations in Y.")

  yNA <- is.na(Y)


  ## Initialize
  beta <- double(mbeta)
  js <- (mY*mA)+seq_len(mB)

  ## IRLS iterations
  for(iter in seq_len(control$maxit)) {
    if ( control$trace ) message(paste0("Iteration = ", iter))

    ## Construct the (weighted) normal equation XTWX = XTWZ
    XTWX <- matrix(0,mbeta,mbeta)
    XTWZ <- double(mbeta)
    for(k in seq_len(mY)) {
      is <- ((k-1)*mA)+seq_len(mA)

      ## IRLS weights and adjusted response
      y <- Y[,k]
      if(iter>1) {
        eta <-  drop(A%*%beta[is] + B%*%beta[js] + offset)
      } else {
        eta <- local({
          weights <- if(is.matrix(weights)) weights[,k] else weights
          weights <- rep(weights,length.out=nobs)
          weights[yNA[,k]] <- 0
          eval(family$initialize)
          family$linkfun(mustart)
        })
      }
      mu <-  family$linkinv(eta)
      mu.eta <- family$mu.eta(eta)
      if ( any(mu.eta == 0) )
        stop(paste("Divisor mu.eta = 0 for iter =", iter,"and k =",k))
      varmu <- family$variance(mu)
      if ( any(varmu == 0) )
        stop(paste("Divisor varmu = 0 for iter =", iter,"and k =",k))
      z <- (eta-offset) + (y-mu)/mu.eta
      w <- if(is.matrix(weights)) weights[,k] else weights
      w <- w*(mu.eta^2)/varmu
      w[yNA[,k]] <- 0

      ## Contribution to XTWX
      XTWX[is,is] <- crossprod(A,w*A)
      XTWX[js,is] <- t(XTWX[is,js] <- crossprod(A,w*B))
      XTWX[js,js] <- XTWX[js,js]+crossprod(B,w*B)

      ## Contribution to XTWZ
      wz <- w*z
      wz[yNA[,k]] <- 0
      XTWZ[is] <- crossprod(A,wz)
      XTWZ[js] <- XTWZ[js]+crossprod(B,wz)
    }
    beta.old <- beta

    ## Solve by Cholesky decomposition
    U <- suppressWarnings(chol(XTWX,pivot=TRUE))
    rank <- attr(U,"rank")
    p <- attr(U,"pivot")
    if(rank < mbeta) {
      U <- U[seq_len(rank),seq_len(rank)]
      p <- p[seq_len(rank)]
    }
    beta0 <- backsolve(U,backsolve(U,XTWZ[p],transpose=TRUE))
    beta <- double(mbeta)
    beta[p] <- beta0
    converged <- sqrt(crossprod(beta-beta.old)) < control$epsilon
    if (converged) break
  }

  ## Compute the deviance
  deviance <- 0
  Mu <- Y
  for(k in seq_len(mY)) {
    is <- ((k-1)*mA)+seq_len(mA)
    eta <-  A%*%beta[is] + B%*%beta[js] + offset
    Mu[,k] <-  family$linkinv(eta)
    wts <- if(is.matrix(weights)) weights[,k] else rep(weights,length.out=length(y))
    devres <- family$dev.resid(Y[,k],Mu[,k],wts)
    deviance <- deviance + sum(devres[!yNA[,k]])
  }

  beta <- rep(NA,mbeta)
  beta[p] <- beta0
  V <- matrix(NA,mbeta,mbeta)
  V[p,p] <- chol2inv(U)

  list(coefficients=beta,
       cov=V,
       deviance=deviance,
       fitted.values=Mu,
       prior.weights=weights,
       converged=converged,
       iter=iter)
}

#-----------------------------------------------------------------------------------------

spatial.covariate <- function(nx,ny,theta) {
  ## Define covariance structure
  cov <- Exp.image.cov(grid=list(x=seq.int(nx),y=seq.int(ny)),theta=theta,setup=TRUE)
  ## Simulate covariate over grid
  sim.rf(cov)
}

#-----------------------------------------------------------------------------------------

sim.covariates <- function(m,nx,ny,theta,V=diag(1,m,m),covarPrefix="x") {
  d <- matrix(0,nx*ny,m)
  for(j in seq_len(m)) {
    simcov <- spatial.covariate(nx,ny,sample(theta,1))
    d[,j] <-  scale(as.vector(t(simcov)),center=TRUE,scale=TRUE)  # Sam added transpose.
  }
  d <- d %*% chol(V)
  colnames(d) <- paste0(covarPrefix,seq.int(m))
  as.data.frame(d)
}

#-----------------------------------------------------------------------------------------

covariates.asRasterStack <- function(nRows, nCols, data, xyGiven=FALSE) {

  # Convert data.frame data into raster stack with one layer per covariate.

  # Number of columns in the data.
  numDataCols <- dim(data)[2]

  # Use x and y to specify raster?
  if ( xyGiven ) {
    indXYCols <- which(c("x","y") %in% colnames(data))
    covarNames <- colnames(data)[-indXYCols]
    rsCovars <- stack()
    for ( i in 1:length(covarNames) ) {
      rlCovar <- rasterFromXYZ(data[ ,c("x","y",covarNames[i])])
      rsCovars <- addLayer(rsCovars, rlCovar)
    }
  } else {
    nCovars <- numDataCols
    rlCovar <- raster(nrows=nRows, ncols=nCols, xmn = 0.5, xmx=nCols+0.5, ymn=0.5, ymx=nRows+0.5)
    rsCovars <- stack()
    for ( i in 1:nCovars ) {
      covar <- matrix(data[ ,i], nrow=nRows, ncol=nCols)
      rlCovar <- setValues(rlCovar, apply(covar,2,rev))
      rsCovars <- addLayer(rsCovars, rlCovar)
    }
    names(rsCovars) <- colnames(data)
  }

  # Return raster stack.
  return(rsCovars)
}

#-----------------------------------------------------------------------------------------

makeCovariates <- function(nRows, nCols, nCovars, theta, covarNames=paste0("x",1:nCovars),
                           checkPlot=FALSE){

  # Simulate the specified number of covariates for the given domain (nRows + nCols).

  covars <- sim.covariates(nCovars,nCols,nRows,theta)
  colnames(covars) <- covarNames
  dfData <- data.frame(x=rep(1:nCols, each=nRows), y=rep(nRows:1, times=nCols), covars)
  rsCovars <- covariates.asRasterStack(nRows, nCols, dfData, TRUE)

  if ( checkPlot ) {
    opar <- par(mfrow=c(2,2))
    plot(rsCovars, asp=1, xlab=names(rsCovars))
    par(opar)
  }

  return(list(data=dfData,raster=rsCovars))

}

#-----------------------------------------------------------------------------------------

cellFromMask <- function(x, mask, maskValue=NA, index=TRUE) {

  # Get the cell numbers of x that are defined to be included by the mask.  Note that x
  # and mask can have different raster structures as long as mask extent is >= x extent.
  #
  # Arguments ...
  # x:         a raster from which the cell numbers are determined.
  # mask:      raster that defines included and excluded parts of the extent.
  # maskValue: value that is used in the cells for the excluded areas of the extent.
  # index:     TRUE returns the cell number of the included cells, FALSE returns a logical
  #            vector with a TRUE for included cells and a FALSE for excluded cells.

  # Get the centre points of each cell of x.
  centres <- data.frame(xyFromCell(x, 1:ncell(x)), stringsAsFactors = FALSE)

  # Get the cells of x that need to be included.
  if ( !is.null(mask) ) {
    # Use mask to determine centre points to be included or excluded.
    maskValues <- drop(davesExtract.v3(mask, centres))
    if ( is.na(maskValue) ) {
      indIncludeCells <- !is.na(maskValues)
    } else {
      indIncludeCells <- maskValues != maskValue
    }
  } else {
    # Include all points.
    indIncludeCells <- rep(TRUE, times=dim(centres)[1])
  }

  # Set return value.
  if ( index ) {
    # Index vector (i.e. vector of cell numbers that are to be included.)
    return(which(indIncludeCells))
  } else {
    # Logical vector of whether or not the cell is included (in raster cell order).
    return(indIncludeCells)
  }

}

#-----------------------------------------------------------------------------------------

diffCoeffs <- function(true, est, abs=TRUE, scale=TRUE) {

  # Calculates the difference between the true and estiamted coefficients.
  #
  # true:  a vector containing the true value of the coefficients (numCoeffs)
  # est:   a vector containing the estimated value of the coefficients (numCoeffs)
  # abs:   if true, return the absoulute value.
  # scale: if true, scale the result by the true coefficient values.

  # Check arguments are of the same length.
  numCoeffs <- length(true)
  if ( length(est) != numCoeffs ) stop("Coefficient vectors are not the same length.")

  # Difference.
  retVal <- true - est

  # Scale result?
  if ( scale ) retVal <- retVal/true

  # Absolute result?
  if ( abs ) retVal <- abs(retVal)

  # Return value.
  return(retVal)

}

#-----------------------------------------------------------------------------------------

proj4str.longlat <- function() {

  # Returns the string that gives the long-lat projection string.

  return("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
}

#-----------------------------------------------------------------------------------------

proj4str.laea <- function(units="m", lonOrigin="0", latOrigin="-90") {

  # Returns the string that gives the laea projection string.

  projStr <- paste("+proj=laea +lat_0=", latOrigin, " +lon_0=", lonOrigin,
                   " +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0",
                   " +units=", units, sep="")
  return(projStr)
}

#-----------------------------------------------------------------------------------------

makeFileName <- function(fnames, dir, ext, dirSep="/", extSep=".") {
  
  # Make a full file name/s from the given components.  If fnames is a vector, then a full
  # file name will be made for each element of the vector and a vector returned.
  
  fileNames <- fnames
  
  if ( ! missing(dir) ) {
    # Add a directory path to the file name.
    fileNames <- paste0(dir, dirSep, fileNames)
  }
  
  if ( ! missing(ext) ) {
    # Add an extension to the file name.
    fileNames <- paste0(fileNames, extSep, ext)
  }
  
  # Return file name.
  return(fileNames)
}

#-----------------------------------------------------------------------------------------

colourBlindRGB <- function(col=c("black", "yellow", "blue", "red", "green", "skyBlue",
                                 "orange", "purple", "white"), maxColorValue=255) {
  
  # Returns the RGB for these special colour blind friendly versions of the colours.
  # col can be a single value or a vector.  It can be the character strings above or integer
  # values 1:8 corresponding to the strings above.  Anything else will result in an error.
  
  # Character string versions.
  namesValidColours <- c("black", "yellow", "blue", "red", "green", "skyBlue", "orange",
                         "purple", "white")
  numValidColours <- length(namesValidColours)
  
  # RGB values
  rgbValidColours <- c(rgb(0, 0, 0, maxColorValue = maxColorValue),       # black
                       rgb(240, 228, 66, maxColorValue = maxColorValue),  # yellow (goldy yellow)
                       rgb(0, 114, 178, maxColorValue = maxColorValue),   # blue (marine blue)
                       rgb(213, 94, 0, maxColorValue = maxColorValue),    # red (tomato red)
                       rgb(0, 158, 115, maxColorValue = maxColorValue),   # green (blue green)
                       rgb(86, 180, 233, maxColorValue = maxColorValue),  # sky blue
                       rgb(230, 159, 0, maxColorValue = maxColorValue),   # orange
                       rgb(204, 121, 167, maxColorValue = maxColorValue), # purple (reddish purple)
                       rgb(255, 255, 255, maxColorValue = maxColorValue)) # white
  names(rgbValidColours) <- namesValidColours
  
  # Check requested colours are valid.
  if ( inherits(col, "integer") ) {
    if ( max(col) > numValidColours || min(col) < 1 ) stop("Invalid colour number in argument 'col'.")
  } else if ( inherits(col, "character") ) {
    for ( colour in col ) match.arg(colour, namesValidColours, several.ok=FALSE)   
  } else {
    stop("Invalid class for argument 'col'.")
  }
  
  # Return RGB values.
  return(rgbValidColours[col])
  
}