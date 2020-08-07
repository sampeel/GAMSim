#-----------------------------------------------------------------------------------------

setupSim <- function(nRows, nCols, nEnvirs, nBiases, nSpecies, seed=10, plotCheck = FALSE){
          
  # Setup for all experiments.  Should only need to be run once.
  # Returns list containing a domain object and a cells object.
  
  # Create domain.
  simExt <- extent(0.5,nCols+0.5, 0.5,nRows+0.5)
  domainObj <- makeDomain(simExt, nrows = nRows, ncols = nCols, proj = proj4str.longlat())
  
  # Simulate covariates.
  cellsObj <- simCh3Covars(domainObj, nEnvirs, nBiases, seed[1], plotCheck)
  
  # Set species.
  cellsObj <- setSpecies(cellsObj, paste0("sp", 1:nSpecies))
  
  # Return values.
  return(list(domain = domainObj, cells = cellsObj))
  
}

#-----------------------------------------------------------------------------------------

runExperiment <- function(cellsObj, domainObj, nSamps, expType = c("mixed","quad"), 
                          corBias = NULL, useMethods = c("MGCV","JAGS","TMB"), 
                          useDataMgcv = c("PA", "PO","AB","POPA"),  centreData = FALSE, 
                          sim=0,  alpha = rep(0.0, cellsObj$numSpecies), beta = NULL, 
                          gamma = NULL, delta = NULL, saveDir=getwd(),
                          plotCheck = FALSE, plotResults = TRUE, quiet = TRUE) {
  
  
  # Run the GAM experiment on the given domain and with the given covariates (in cellsObj).
  # The number of AB/PA samples can be altered.  The number of PO observations is not varied.
  # Plots are only for checking purposes.  Results are saved so that paper quality plots
  # can be made externally.
  # Returns a list of fits for each method used.  Also saves results to a directory created
  # within this function.
  #
  # Arguments ...
  # cellsObj:    a cells object that already contains cells and covariate values.
  # domainObj:   a domain object (needed for creating sample locations)
  # nSamps:      the number of survey sample locations to simulate (random location only!)
  # expType:     species component curves are either a mixture of curve types or all quadratic.
  #              (FYI: all bias component curves are quadratic and do not vary with this setting)
  # corBias:     an integer vector of length 2 that sets the corBias[2] sampling bias 
  #              covariate to be the same values as the corBias[1] species covariate.  This 
  #              simulates the sampling bias being correlated (i.e. the same covariate) with 
  #              the species intensity.  Set this equal to NULL or c(0,0) if no 
  #              correlation is required.
  # useMethods:  the methods to use to produce results (can be one or more of these value).
  # useDataMgcv: the type of data to be used in the mgcv method (i.e. function 'gam'). 
  #              Only one type at a time is valid.  Only used if "MGCV" included in useMethods.
  # centreData:  centre the covariate data based on each covariate's mean value in data.
  # sim:         simulation number. A quick way to generate a few different by changing 
  #              the random seed by this amount simulations (FYI: sim=0 for thesis sim).
  # alpha:       see function 'simCh3Intensity'.  Only used when expType == "quad".
  # beta:        see function 'simCh3Intensity'.  Only used when expType == "quad".
  # gamma:       see function 'simCh3SampBias.
  # delta:       see function 'simCh3SampBias.
  # saveDir:     the experiment directory where it will save the results (*.RData files)
  #              set to NULL for no save.
  # plotCheck:   plot check plots for simulated data.
  # plotResults: plot check plots for results.
  # quiet:       if FALSE, will print sim results at different stages of the experiment.


  # Check objects have been setup.
  if ( is.null(cellsObj) || is.null(cellsObj$namesSpecies) || is.null(domainObj) ) 
    stop("Please run function 'setupSim' first.")

  # Check multi-option arguments.
  expType <- match.arg(expType)
  useMethods <- match.arg(useMethods, several.ok = TRUE)
  useDataMgcv <- match.arg(useDataMgcv, several.ok = TRUE)
  nMethods <- length(useMethods)
  if ( length(useDataMgcv) == 1 ) {
    useDataMgcv <- rep(useDataMgcv, times=nMethods)
  } else if ( length(useDataMgcv) != nMethods) {
    stop("The length of 'useDataMgcv' should be one or the same as 'useMethods'.")
  }
  
  # Experiment set up.
  if ( is.null(corBias) || ( corBias[1] == 0 && corBias[2] == 0 )) {
    corBias <- NULL
  } else if ( corBias[1] != 0 && corBias[2] != 0 ) {
    # Set up duplicated sampling bias covariate.
    nameCorBias <- cellsObj$namesBiases[corBias[2]]
    nameCorEnvir <- cellsObj$namesCovars[corBias[1]]
    cellsObjBiasesOrig <- cellsObj$biases[ ,nameCorBias]
    cellsObj$biases[ ,nameCorBias] <- cellsObj$covars[ ,nameCorEnvir]
  } else {
    stop("Correlation of bias indicator 'corBias' not set up properly.  Please check.")
  }

  # Create save directory, if necessary.
  if ( ! is.null(saveDir) ) {
    if ( ! dir.exists(saveDir) ) dir.create(saveDir)
    if ( ! quiet ) message("Save to directory: ", saveDir)
  }
  
  # Simulate the true intensity and the number of individuals per cell for each species.
  if ( expType == "quad") {
    if ( is.null(beta) ) stop("For experiment type 'quad', beta needs to be given.")
    cellsObj <- simCh3Intensity(cellsObj, seed = c(1,9+sim), alpha = alpha, beta = beta, 
                                plotCheck = plotCheck)
  } else {
    cellsObj <- simCh3Intensity(cellsObj, seed = c(1,9+sim), plotCheck = plotCheck)
  }
  if ( ! quiet ) {
    msg <- paste0("Number of individuals in species ", 1:cellsObj$numSpecies, ": ", 
                  apply(cellsObj$N, 2, sum), collapse = "\n")
    message(msg)
    message(paste0("Using these function types: ", paste0(cellsObj$compFuncs, collapse=" ")))
  }
  
  # Simulate sample locations and value per species.
  surveysObj <- simCh3Samples(domainObj, cellsObj, nSamps, seed=3+sim, plotCheck = plotCheck)
  if ( ! quiet ) {
    msg <- paste0("Number of individuals sampled for species ", 1:cellsObj$numSpecies, ": ", 
                  apply(surveysObj$N, 2, sum), collapse = "\n")
    message(msg)
    msg <- paste0("Number of sample presences for species ", 1:cellsObj$numSpecies, ": ", 
                  apply(surveysObj$Y, 2, sum), collapse = "\n")
    message(msg)
  }
  
  # Simulate probability of observation per cell per species and simulate observations.
  cellsObj <- simCh3SampBias(cellsObj, 2+sim, gamma, delta, plotCheck=plotCheck)
  if ( ! quiet ) {
    msg <- paste0("Number of observed species ", 1:cellsObj$numSpecies, ": ", 
                  apply(cellsObj$P, 2, sum), collapse = "\n")
    message(msg)
  }
  
  ### Data formation ...   
  
  # Names and numbers of stuff.
  namesSpecies <- cellsObj$namesSpecies
  namesEnvirs <- cellsObj$namesCovars
  namesBiases <- cellsObj$namesBiases
  nBiases <- cellsObj$numBiases

  # Setup data
  d.po <- cbind.data.frame(cellsObj$P[ ,namesSpecies],
                           cellsObj$covars[ ,namesEnvirs,drop=FALSE],
                           cellsObj$biases[ ,namesBiases,drop=FALSE],
                           area = cellsObj$areaCell,
                           isP = TRUE)
  d.pa <- cbind.data.frame(surveysObj$Y[ ,namesSpecies],
                           cellsObj$covars[surveysObj$rowsInCells,namesEnvirs,drop=FALSE],
                           matrix(0, nrow=nSamps, ncol=nBiases, dimnames = list(NULL, namesBiases)),
                           area = surveysObj$areas,
                           isP = FALSE)
  d.ab <- cbind.data.frame(surveysObj$N[ ,namesSpecies],
                           cellsObj$covars[surveysObj$rowsInCells,namesEnvirs,drop=FALSE],
                           matrix(0, nrow=surveysObj$numSamples, ncol=nBiases, 
                                  dimnames = list(NULL, namesBiases)),
                           area = surveysObj$areas, isP = FALSE)
  d.popa <- rbind(d.pa,d.po)
  
  # Centre data, if requested
  if ( centreData ) {
    for ( covar  in namesEnvirs ) {
      d.popa[ ,covar] <- scale(d.popa[ ,covar], scale = FALSE, center = TRUE)
      d.ab[ ,covar] <- scale(d.ab[ ,covar], scale = FALSE, center = TRUE)
      d.pa[ ,covar] <- scale(d.pa[ ,covar], scale = FALSE, center = TRUE)
    }
    for ( bias in namesBiases ) {
      d.popa[d.popa$isP,bias] <- scale(d.popa[d.popa$isP,bias],scale = FALSE,center = TRUE)
      d.po[ ,bias] <- scale(d.po[ ,bias], scale = FALSE, center = TRUE)
    }
  }
  if ( ! quiet ) message("Input data for methods has been made.")
  
  
  ### Run methods.
  retLst <- list(expDir=saveDir, expType=expType, corBias=corBias, useMethods=useMethods, 
                 useDataMgcv=useDataMgcv, cellsObj=cellsObj, surveysObj)
  for ( m in 1:nMethods ) {
    # Initialise return value.
    fit <- NULL
    
    if ( useMethods[m] == "MGCV" ) {
      # mgcv::gam ...
      if ( ! quiet ) message("Starting MGCV method with ",useDataMgcv[m]," data.")
      start.time <- Sys.time()
      namesBiasesMGCV <- namesBiases
      if ( useDataMgcv[m] == "AB") {
        data.MGCV <- d.ab
      } else if ( useDataMgcv[m] == "PA" ) {
        data.MGCV <- d.pa
      } else if ( useDataMgcv[m] == "PO" ) {
        data.MGCV <- d.po
        if ( ! is.null(corBias) ) {
          namesBiasesMGCV <- namesBiases[-corBias[2]]
          colCorBias <- which(names(data.MGCV) == namesBiases[corBias[2]])
          data.MGCV <- data.MGCV[ ,-colCorBias]
        }
      }
      
      tryOut <- try(runMGCV(data.MGCV, namesSpecies, namesEnvirs, namesBiasesMGCV, 
                                      useDataMgcv[m], plotResults = plotResults))
      if ( inherits(tryOut, "try-error")) { 
        lstFit.MGCV <- NULL 
      } else {
        lstFit.MGCV <- tryOut
      }
      end.time <- Sys.time()
      total.time.MGCV <- end.time - start.time
      
      # Save results.
      if ( ! is.null(saveDir) ) {
        mgcvFile <- paste0("mgcvResults-n", useDataMgcv[m], NROW(data.MGCV),"-S",simNum)
        saveResults <- makeFileName(mgcvFile, saveDir, "RData")
        useData <- useDataMgcv[m]
        save(lstFit.MGCV, total.time.MGCV, useData, data.MGCV, cellsObj, surveysObj, 
             domainObj, corBias, file = saveResults)
      }
      retLst[[paste0("fit.MGCV.",useDataMgcv[m])]] <- lstFit.MGCV
      if ( ! quiet ) message("Finished MGCV method with ",useDataMgcv[m]," data.")
    
    } else if ( useMethods[m] == "TMB" ) {
      # popaGAM::popaTMB
      if ( ! quiet ) message("Starting TMB method.")
      start.time <- Sys.time()
      data.TMB <- d.popa
      
      tryOut <- try(runTMB(data.TMB, namesSpecies, namesEnvirs, namesBiases, 
                           plotResults = plotResults))
      if ( inherits(tryOut, "try-error")) { 
        fit.TMB <- NULL 
      } else {
        fit.TMB <- tryOut
      }
      end.time <- Sys.time()
      total.time.TMB <- end.time - start.time
      
      # Save results
      if ( ! is.null(saveDir) ) {
        tmbFile <- paste0("tmbResults-nPA",nSamps,"-S",simNum)
        saveResults <- makeFileName(tmbFile, saveDir, "RData")
        save(fit.TMB, total.time.TMB, data.TMB, cellsObj, surveysObj, domainObj, file = saveResults)
      }
      retLst$fit.TMB <- fit.TMB
      if ( ! quiet ) message("Finished TMB method.")
    
    } else if ( useMethods[m] == "JAGS" ) {
      # popaGAM::popaJAGS ...
      if ( ! quiet ) message("Starting JAGS method.")
      start.time <- Sys.time()
      data.JAGS <- d.popa
      
      tryOut <- try(runJAGS(data.JAGS, namesSpecies, namesEnvirs, namesBiases, 
                            plotResults = plotResults))
      if ( inherits(tryOut, "try-error")) { 
        fit.JAGS <- NULL 
      } else {
        fit.JAGS <- tryOut
      }
      end.time <- Sys.time()
      total.time.JAGS <- end.time - start.time
      
      # Save results
      if ( ! is.null(saveDir) ) {
        jagsFile <- paste0("jagsResults-nPA", nSamps, "-S",simNum)
        saveResults <- makeFileName(jagsFile, saveDir, "RData")
        save(fit.JAGS, total.time.JAGS, data.JAGS, cellsObj, surveysObj, domainObj, file = saveResults)
      }
      retLst$fit.JAGS <- fit.JAGS
      if ( ! quiet ) message("Finished JAGS method.")

    } else {
      warning("Unrecognised method ",useMethods[m]," in argument 'useMethods'.")
    }
  }
  
  # Reset original cell bias covariate values, if necessary.
  if ( ! is.null(corBias) ) cellsObj$biases[ ,nameCorBias] <- cellsObjBiasesOrig

  # Return value ...
  return(retLst)
  
}
  
#-----------------------------------------------------------------------------------------

runMGCV <- function(data, namesSpecies, namesEnvirs, namesBiases=NULL, 
                    dataType=c("AB", "PA", "PO"), offset = log(data$area), plotResults = TRUE) {
  
  # Fit GAM using R package mgcv with specified type of data (and approriate family and 
  # link function).  
  # Returns a list of gam fits, one fit per species.
  
  # Check arguments.
  if ( is.null(data) ) stop("The data has not been given.")
  namesDataCols <- names(data)
  if ( is.null(namesSpecies) ) 
    stop("The names of the species have not been given.")
  if ( ! all(namesSpecies %in% namesDataCols) ) 
    stop("One or more species names do not match data column names.")
  if ( is.null(namesEnvirs) ) 
    stop("The species covariate names have not been given.")
  if ( ! all(namesEnvirs %in% namesDataCols) ) 
    stop("One or more environment covariate names do not match data column names.")
  dataType <- match.arg(dataType)
  if ( dataType == "PO" ) {
    if ( is.null(namesBiases) ) 
      stop("The sampling bias covariate names must be given when dataType = 'PO'.")
    if ( ! all(namesBiases %in% namesDataCols) ) 
      stop("One or more sampling bias covariate names do not match data column names.")
  }

  # What are the gam terms.
  namesCovars <- namesEnvirs
  if (  dataType == "PO" ) namesCovars <- c(namesEnvirs, namesBiases)
  gam.terms <- paste0("s(", namesCovars, ",bs='cr')", collapse = " + ")

  # # Set up result plotting, if necessary.
  # nCovars <- length(namesCovars)
    nSpecies <- length(namesSpecies)
  # if ( plotResults ) {
  #   oma.new <- c(0,0,4,0)
  #   mar.new <- c(5.1,4.1,0.1,0.5)
  #   opar <- par(mfcol=c(nCovars,nSpecies), oma=oma.new, mar=mar.new)
  # }
  
  # Set up link family for given data type.
  if ( dataType == "PO" || dataType == "AB" ) {
    gam.family <- poisson(link="log")
  } else if ( dataType == "PA") {
    gam.family <- binomial(link="cloglog")
  }
  
  # Fit GAM for each species.
  gamLst <- vector("list",nSpecies)
  for ( k in 1:nSpecies ) {
    # GAM
    species <- namesSpecies[k]
    gam.formula <- as.formula(paste0(species, "~", gam.terms))
    gamLst[[k]] <- gam(gam.formula, gam.family, data, offset = offset, select = TRUE)
    
    # Plot gam.
    if ( plotResults ) {
      plot.gam(gamLst[[k]], pages=1, scale=0, 
               main=paste0("MGCV estimates (",namesSpecies[k],")"))
    }
  }

  # Return value.
  return(gamLst)
  
}

#-----------------------------------------------------------------------------------------

runJAGS <- function(data, namesSpecies, namesEnvirs, namesBiases=NULL, 
                    area = data$area, diagonalize = TRUE, plotResults = TRUE) {
  
  # Fit popaJAGS using R package popaGAM. Assumes cubic regression splines as the basis 
  # for the smoothing (see mgcv::smooth.terms for more info).
  # Returns a popaGAM fit object.

  # Check arguments.
  if ( is.null(data) ) stop("The data has not been given.")
  namesDataCols <- names(data)
  if ( is.null(namesSpecies) ) 
    stop("The names of the species have not been given.")
  if ( ! all(namesSpecies %in% namesDataCols) ) 
    stop("One or more species names do not match data column names.")
  if ( is.null(namesEnvirs) ) 
    stop("The species covariate names have not been given.")
  if ( ! all(namesEnvirs %in% namesDataCols) ) 
    stop("One or more environment covariate names do not match data column names.")
  if ( is.null(namesBiases) ) 
    stop("The sampling bias covariate names have not been given.")
  if ( ! all(namesBiases %in% namesDataCols) ) 
    stop("One or more sampling bias covariate names do not match data column names.")
  
  # Set up formulae for TMB with smooth terms.
  lambda.terms <- paste0("s(", namesEnvirs, ",bs='cr')", collapse = " + ")
  lambda.response <- paste0("cbind(",paste0(namesSpecies, collapse=","),")")
  lambda.formula <- as.formula(paste0(lambda.response, " ~ ", lambda.terms))
  bias.terms <- paste0("s(", namesBiases, ",bs='cr')", collapse = " + ")
  bias.formula <- as.formula(paste0("isP ~ ", bias.terms))
  
  # Set up control for rjags.
  popa.control <- gam.control()
  popa.control$n.chains <- 4
  popa.control$n.adapt <- 500
  
  # Set up GAM popa model for use in rjags. (equivalent popaGAM:::popaJAGS + jagsModel functions)
  resLst <- popaGAM:::popaJAGS(lambda.formula, bias.formula, data, area, diagonalize)
  cat(resLst$jags$model,file=paste0(getwd(),"/tmpJmod.txt"), sep="")
  # resLst <- popa.model(lambda.formula, bias.formula, data = data, area = area, 
  #                      control=popa.control, diagonalize = diagonalize)
  # 

  # Run rjags to estimate parameters of GAM popa model.
  samples <- jags.samples(resLst$jags, resLst$jags$vars, n.iter = 25000, thin = 5) # n.iter = 50000, thin = 10)
  save(data, resLst, samples, file = paste0(getwd(),"/jagsSamples.RData"))

  # Plot checks (FYI: part of method so don't just do when plotCheck == TRUE!!!)
  if ( plotResults ) {
    plot.samples(samples$betaA, 1:length(namesSpecies), namesSpecies,"betaA")
    plot.samples(samples$betaA, titleResponse = "Bias")
  }

  # Form popa object ready for plotting from JAGS results
  # names(samples)[1] <- "betaA"
  # names(samples)[3] <- "rhoA"
  # names(samples)[4] <- "betaB"
  # names(samples)[6] <- "rhoB"
  # names(resLst)[1] <- "p.gam"
  fit.JAGS <- popaGAM:::jags2popaGAM(resLst, samples, burnin = 2000)
  #fit.JAGS <- jags2popaGAM(resLst, samples, burnin = 2000)
  
  # Plot results ...
  if ( plotResults ) {
    oma.new <- par("oma")
    oma.new[3] <- 2.1
    mar.new <- par("mar")
    mar.new[3] <- 1.1
    opar <- par(oma=oma.new, mar=mar.new)
    for ( k in 0:length(namesSpecies) ) {
      popaGAM:::plot.popaGAM(fit.JAGS,species=k,pages=1)
      if ( k == 0 ) {
        title("JAGS estimates (sampling bias)", outer=TRUE)
      } else {
        title(paste0("JAGS estimates (", namesSpecies[k], ")"), outer=TRUE)
      }
    }
    par(opar)
  }
  
  # Return value.
  return(fit.JAGS)
  
}

#-----------------------------------------------------------------------------------------

runTMB <- function(data, namesSpecies, namesEnvirs, namesBiases, 
                   area = data$area, plotResults = TRUE) {
  
  # Fit popaTMB using R package popaGAM. Assumes cubic regression splines as the basis 
  # for the smoothing (see mgcv::smooth.terms for more info).
  # Returns a popaGAM fit object.
  
  # Check arguments.
  if ( is.null(data) ) stop("The data has not been given.")
  namesDataCols <- names(data)
  if ( is.null(namesSpecies) ) 
    stop("The names of the species have not been given.")
  if ( ! all(namesSpecies %in% namesDataCols) ) 
    stop("One or more species names do not match data column names.")
  if ( is.null(namesEnvirs) ) 
    stop("The species covariate names have not been given.")
  if ( ! all(namesEnvirs %in% namesDataCols) ) 
    stop("One or more environment covariate names do not match data column names.")
  if ( is.null(namesBiases) ) 
    stop("The sampling bias covariate names have not been given.")
  if ( ! all(namesBiases %in% namesDataCols) ) 
    stop("One or more sampling bias covariate names do not match data column names.")

  # Set up formulae for TMB with smooth terms.
  lambda.terms <- paste0("s(", namesEnvirs, ",bs='cr')", collapse = " + ")
  lambda.response <- paste0("cbind(",paste0(namesSpecies, collapse=","),")")
  lambda.formula <- as.formula(paste0(lambda.response, " ~ ", lambda.terms))
  bias.terms <- paste0("s(", namesBiases, ",bs='cr')", collapse = " + ")
  bias.formula <- as.formula(paste0("isP ~ ", bias.terms))

  # Run TMB version of popaGAM.
  fit.TMB <- popaGAM:::popaTMB(lambda.formula, bias.formula, data, area)
  
  # Plot results ...
  if ( plotResults ) {
    oma.new <- par("oma")
    oma.new[3] <- 2.1
    mar.new <- par("mar")
    mar.new[3] <- 1.1
    opar <- par(oma=oma.new, mar=mar.new)
    for ( k in 0:length(namesSpecies) ) {
      popaGAM:::plot.popaGAM(fit.TMB,species=k,pages=1)
      if ( k == 0 ) {
        title("TMB estimates (sampling bias)", outer=TRUE)
      } else {
        title(paste0("TMB estimates (", namesSpecies[k], ")"), outer=TRUE)
      }
    }
    par(opar)
  }

  # Return value.
  return(fit.TMB)
  
}

#-----------------------------------------------------------------------------------------

simCh3Covars <- function(domainObj, nEnvirs=3, nBiases=2, seed=10, plotCheck=FALSE) {
  
  # Simulate domain (res = c(1,1)) and covariates for chapter 3 work (GAM chapter).
  # Names of covariates are set within the function (x for species covariates, z for bias).
  # Return cells object with simulated covariates.

  # Initialise stuff ...
  nCovars <- nEnvirs + nBiases
  namesEnvir <- paste0("x",1:nEnvirs)
  namesBias <- paste0("z",1:nBiases)
  namesCovar <- c(namesEnvir, namesBias)
  
  # Make cells in the domain ...
  cellsObj <- makeCells(domainObj$mask, domainObj$maskValue)

  # Make the covariate data
  set.seed(seed)    # Important to run this line each time same covariates are required!
  theta <- sample(10:15,nCovars,replace=TRUE)
  retLst <- makeCovariates(domainObj$nrows, domainObj$ncols, nCovars, theta, namesCovar, 
                           plotCheck)

  # # Make indCommon[r,2] bias covariate the same as the indCommon[r,1] envir covariate
  # # (i.e. simulate use of same covariate in both species distribution and sampling bias)
  # if ( ! is.null(indCommon) ) {
  #   for ( r in 1:NROW(indCommon)) {
  #     retLst$data[ ,namesBias[indCommon[r,2]]][] <- retLst$data[ ,namesEnvir[indCommon[r,1]]][]
  #     retLst$raster[[namesBias[indCommon[r,2]]]] <- retLst$raster[[namesEnvir[indCommon[r,1]]]]
  #   }
  # }
  # 

  # Store covariates in cells object ...
  cellsObj <- makeCovarData(cellsObj,
                            dropLayer(retLst$raster, namesBias),
                            dropLayer(retLst$raster, namesEnvir))
  if ( plotCheck ) {
    opar <- par(mfrow=c(2,3))
    plotCellVals(cellsObj, "covars", titles = namesEnvir)
    plotCellVals(cellsObj, "biases", titles = namesBias)
    par(mfrow=c(1,1))
  }
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

simCh3Intensity <- function(cellsObj, alpha = rep(0.0, cellsObj$numSpecies), beta = NULL,
                            nIndivids = rep(50*cellsObj$numCells, cellsObj$numSpecies),
                            seed=c(1,9), plotCheck=FALSE) {
  
  # Simulate the number of individuals per cell per species.  To do this, first create the
  # true intensity values (and component function values) from the randomly selected functions 
  # (as available in function 'makeGAMTrue').  Then scale the intercept so that approximately
  # the given number of individuals exist in the domain, for each species.  Finally, 
  # simulate the number of individuals in each cell using true intensity values.
  # Returns cells object with "trueLambda", "trueComp", "compFuncs" and "N" set.
  
  # Initialise stuff.
  nSpecies <- cellsObj$numSpecies
  namesSpecies <- cellsObj$namesSpecies
  namesEnvir <- cellsObj$namescovars
  nEnvirs <- cellsObj$numCovars
  
  # Set up true intensity and component functions per species.
  if ( is.null(beta) ) {
    set.seed(seed[1])
    whichFuncs <- matrix(sample(1:6, 6, replace=FALSE), nrow = nEnvirs, ncol=nSpecies,
                         dimnames = list(namesEnvir, namesSpecies))
  } else {
    dimBeta <- dim(beta)
    if ( length(dimBeta) == 3 ) {
      # Use quadratic form for all functions.
      whichFuncs <- matrix(7, nrow = nEnvirs, ncol=nSpecies,
                           dimnames = list(namesEnvir, namesSpecies)) 
    } else {
      # FYI: full checking of beta is done within makeIntensity.GAM function.
      stop("Argument 'beta' should have three dimensions, please check.")
    }
  }
  cellsObj <- makeIntensity.GAM(cellsObj, whichFuncs, alpha, beta, namesSpecies)
  
  # Scaling of alpha and true intensity to get approx nIndivids across domain, per species.
  for ( k in 1:nSpecies ) {
    species <- namesSpecies[k]
    scaleLambda <- nIndivids[k] / sum(cellsObj$trueLambda[ ,species])
    cellsObj$trueLambda[ ,species] <- cellsObj$trueLambda[ ,species] * scaleLambda 
    alpha[k] <- alpha[k] + log(scaleLambda)
    cellsObj$trueComp[ ,species,1] <- alpha[k]
  }
  if ( plotCheck ) 
    plotCellVals(cellsObj, "trueLambda", 
                                titles = paste0("True lambda (", namesSpecies,")"))
  
  # Plot component functions.
  if ( plotCheck ) {
    # Plot true component curves (all on one page).
    opar <- par(mfrow=c(nEnvirs,nSpecies), oma=c(0,0,0,0), mar=c(4,3,2,1))
    for ( ev in 1:nEnvirs ) {
      # Covariate values ...
      x <- cellsObj$covars[ ,ev]
      
      # Plot component for each species
      for ( k in 1:nSpecies ) {
        y <- cellsObj$trueComp[ ,k,ev+1]
        plot(x, y, xlab=namesEnvir[ev], 
             main=paste0("True component functions (",namesSpecies[k], ")"))
      }
    }
    par(opar)
  }
  
  
  # Simulate numbers of each species in each cell.
  set.seed(seed[2])
  cellsObj <- makeNumIndivids(cellsObj)
  if ( plotCheck ) plotCellVals(cellsObj, "N", 
                                titles = paste0("Simulated num individuals (", namesSpecies,")"))

  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

simCh3Samples <- function(domainObj, cellsObj, nSamps=200, sampAreaRange = c(0.01,0.02),
                          seed = 3, plotCheck = FALSE){
  
  # Simulate samples (assumed to be random across whole domain with 100% detection).
  # Returns newly created surveys object.
  
  # Make sample locations
  minSampleArea <- sampAreaRange[1]
  maxSampleArea <- sampAreaRange[2]
  set.seed(seed)
  surveysObj <- makeSampleLocations(domainObj$mask, cellsObj, nSamps, minSampleArea,
                                    maxSampleArea)
  
  # Assign gear types to samples (makeNSampled function needs gear info, so fake it).
  nGears <- 1
  namesGears <- paste0("gt", 1:nGears)
  surveysObj <- assignGearUsedSamples(surveysObj, nGears, namesGears, "rand")
  
  # Make samples (both AB and PA data)
  surveysObj <- makeNSampled(surveysObj, cellsObj$N, cellArea = cellsObj$areaCell)
  
  # Check simulated data
  if ( plotCheck ) {
    namesSpecies <- cellsObj$namesSpecies
    for ( species in namesSpecies ) {
      # Abundance
      plotCellVals(cellsObj, "N", species, title="")
      thisSpeciesN <- surveysObj$N[ ,species]
      points(surveysObj$xy, pch=as.character(thisSpeciesN), col="black")
      title(paste0("Simulated number sampled (", species,")"))
      
      # Presence-absence
      plotCellVals(cellsObj, "N", species, title="")
      thisSpeciesY <- surveysObj$Y[ ,species] + 1
      pchChar <- c("A","P")
      points(surveysObj$xy, pch=pchChar[thisSpeciesY], col="black")
      title(paste0("Simulated presence-absence (", species,")"))
    }
  }

  # Return value.
  return(surveysObj)
  
}

#-----------------------------------------------------------------------------------------

simCh3SampBias <- function(cellsObj, seed = 2, gamma = c(-4.9,-5.5), 
                           delta = c(0,-0.1,-0.4,0.2), plotCheck = FALSE){

  # Set up probability of observation for sampling bias (function b_k(s) in Fithian, et al 2015).
  # Assumes b_k(s) is a quadratic function for species k and that only the intercepts differ 
  # across species. Returns cells object with "trueB" and "trueBComp" (for each covariate),
  # and the simulated data "P" and "PO".
  #
  # Arguments ...
  # cellsObj:  a cells object that already contains cells and covariate values.
  # seed:      the PO data is generated randomly using this random seed.
  # gamma:     an nSpecies length vector of intercepts for the function b_k(s) 
  # delta:     an nBiases x 2 length vector.  E.g. the sample bias covariate functions are
  #                h_1(z_1(s)) = delta[1] z_1(s) + delta[nBiases+1] z_1(s)^2
  #                h_2(z_2(s)) = delta[2] z_2(s) + delta[nBiases+2] z_2(s)^2
  #            etc. for the number of sampling bias covariates specified in cells object.
  # plotCheck: plot check plots for simulated data.

  
  # Initialise stuff ...
  namesBias <- names(cellsObj$biases)
  nBiases <- length(namesBias)
  namesSpecies <- cellsObj$namesSpecies
  nSpecies <- cellsObj$numSpecies
  namesGears <- "gt1"
  nGears <- length(namesGears)
  nCells <- cellsObj$numCells
  
  # Check length of gamma and delta.
  if ( is.null(gamma) ) stop("Argument 'gamma' is NULL.")
  if ( length(gamma) != nSpecies )  stop("Length of argument 'gamma' is incorrect.")
  if ( is.null(delta) ) stop("Argument 'delta' is NULL.")
  if ( length(delta) != (nBiases * 2) ) stop("Length of argument 'delta' is incorrect.")
  
  # Set up formula for sampling bias (this is hardwired to be quadratic)  
  bias.formula <- as.formula(paste("~",paste(namesBias, collapse="+"),
                                   "+", paste("I(",namesBias, "^2)", collapse="+")))
  namesBiasTerms <- attr(terms(bias.formula),"term.labels")
  nBiasTerms <- length(namesBiasTerms)
  namesBiasCoeffs <- c("(bias)", namesBiasTerms)
  
  # Set up coefficients for nSpecies sample bias functions (gamma[k] + h_1(z_1(s)) + ... )
  coeffs <- matrix(nrow=nBiasTerms+1, ncol=nSpecies,
                  dimnames=list(namesBiasCoeffs, namesSpecies))
  coeffs[seq.int(2,nBiasTerms+1),1] <- delta  
  coeffs[-1, ] <- coeffs[-1,1]              # Covariate coefficients SAME for all species
  coeffs[1, ] <- gamma

  # Simulate probability of observation true values.
  cellsObj <- makeProbObs(cellsObj, bias.formula, coeffs["(bias)", ], coeffs[-1,1])
  if ( plotCheck )
    plotCellVals(cellsObj, vals="trueB", 
                 titles=paste0("Prob of observation (",namesSpecies,")"))
  
  # Get components of sampling bias.
  namesGamma <- paste0("(bias.", namesSpecies, ")")
  cellsObj$trueBComp <- matrix(0, nrow = nCells, ncol = nSpecies + nBiases,  
                              dimnames = list(NULL, c(namesGamma,paste0("fb",1:nBiases))))
  for ( k in 1:nSpecies ) {
    # Species specific intercepts.
    cellsObj$trueBComp[ ,k] <- coeffs[1,k]
  }
  for ( bs in 1:nBiases ) {
    # Component for each bias covairate.
    x <- cellsObj$biases[ ,bs]
    y <- coeffs[bs+1,1] * x + coeffs[bs+1+nBiases,1] * x * x
    cellsObj$trueBComp[ ,nSpecies+bs] <- y
  }
  
  # Set up probability of detection (NB: this is just to reuse old function).
  namesDetCoeffs <- paste0("z", 1:nGears)
  zeta <- matrix(0, nGears, nSpecies, dimnames = list(namesDetCoeffs, namesSpecies))  # PrDet = 1
  cellsObj <- setProbDet(cellsObj, zeta, namesGears)
  
  # Simulate presence-only observations.
  set.seed(seed)
  cellsObj <- makePOObservations(cellsObj, 0, nGears, namesGears)
  if ( plotCheck ) {
    plotCellVals(cellsObj, vals="P", 
                 titles = paste0("Simulated num presence-only per cell (", namesSpecies, ")"))
  }
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

plotCompareMethods <- function(expDir, expDataFiles, plotMethods = c("MGCV","JAGS","TMB"),
                               includeBiasCovars = FALSE, includeTrueFunctions = TRUE,
                               plotDevice = "RStudioGD", plotDir = expDir, plotUnits = "cm",
                               plotWidth = 15.8, plotHeight = 20, plotRes = 600, 
                               plotFileName = "estCompFuncs", doPlotIdentifiers=TRUE,
                               posIdentifiers = "topleft", includeRug = TRUE) {
  
  # This function loads results so it can be performed independently of any data simulation 
  # and result estimation.  It uses the data specified by the files in 'expDataFiles' 
  # (assumes these have '.RData' extensions).  Plots the component covariate functions 
  # (as mgcv::plot.gam function does) as estimated by the requested methods.  Can plot
  # estimated bias component curves when they are available (not for MGCV with PA or AB data).
  # Can plot true values for the component covariate functions (at the given data locations).
  # Set plotDevice = "png" (or other valid graphics device) to plot to file.
  
  
  # Check there are the same number of data files as there are methods to plot.
  nMethods <- length(plotMethods)
  if ( length(expDataFiles) != nMethods ) 
    stop("The number of files does not match the number of methods.")

  # Load data 
  for ( i in 1:nMethods ) {
    # Load this method's data.
    loadFile <- makeFileName(expDataFiles[i], expDir, "RData")
    load(loadFile)
  }
  
  # Names and numbers of things.
  namesSpecies <- cellsObj$namesSpecies
  nSpecies <- cellsObj$numSpecies
  namesEnvir <- names(cellsObj$covars)
  nEnvirs <- length(namesEnvir)
  if ( includeBiasCovars ) {
    namesBiases <- names(cellsObj$biases)
    nBiases <- length(namesBiases) 
    covars <- cbind(cellsObj$covars, cellsObj$biases)
  } else {
    namesBiases <- NULL
    nBiases <- 0
    covars <- cellsObj$covars
  }
  nCovars <- nEnvirs + nBiases
  namesCovars <- c(namesEnvir, namesBiases)
  
  # Plot settings
  oma.old <- par("oma")
  mar.old <- par("mar")
  nPlotCols <- nCovars
  if ( plotDevice != "RStudioGD" ) {
    # No need for title space.
    plotToFile <- TRUE
    oma.new <- c(2,2,0,0)
    nPlotRows <- nMethods * nSpecies
    nPlotPages <- 1
  } else {
    plotToFile <- FALSE
    oma.new <- c(2,2,4,0)
    nPlotRows <- nMethods
    nPlotPages <- nSpecies
  }
  mar.new <- c(1.1,2.1,0.1,0.1)
  trueCol <- colourBlindRGB("orange")
  trueSym <- "|"  # 20
  
  # Get values to plot
  gamVals <- vector(mode = "list", length = nSpecies)
  jagsVals <- gamVals
  tmbVals <- gamVals
  for ( k in 1:nSpecies ) {
    # Create plotable objects by writting to an emtpy file.
    # TO DO: work out how to do this for bias covars.
    if ( "MGCV" %in% plotMethods ) {
      pdf(file=NULL)
      gamVals[[k]] <- plot(lstFit.MGCV[[k]], scale=0)
      dev.off()
    }
    if ( "JAGS" %in% plotMethods ) {
      pdf(file=NULL)
      jagsVals[[k]] <- popaGAM:::plot.popaGAM(fit.JAGS, species=k)
      jagsVals[[k]] <- c(jagsVals[[k]], popaGAM:::plot.popaGAM(fit.JAGS, species = 0))
      dev.off()
    }
    if ( "TMB" %in% plotMethods ) {
      pdf(file=NULL)
      tmbVals[[k]] <- popaGAM:::plot.popaGAM(fit.TMB, species = k)
      tmbVals[[k]] <- c(tmbVals[[k]], popaGAM:::plot.popaGAM(fit.TMB, species = 0))
      dev.off()
    }
  }
  
  # Set up plotting to file, if necessary (will print to RStudio plot window otherwise)  
  if ( plotToFile ) {
    plotFile <- makeFileName(plotFileName, expDir, plotDevice)
    argList <- list(filename = plotFile, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
  }
  opar <- par(mfrow=c(nPlotRows,nPlotCols), oma=oma.new, mar=mar.new, cex=0.50)
  plotIdentifier <- paste0("(",letters[1:(nPlotRows*nPlotCols)],")")
  
  iPlotCount <- 0
  # Plot file for each species.
  for ( k in 1:nSpecies ) {
    if ( ! plotToFile ) iPlotCount <- 0
    
    for ( pm in 1:nMethods ) {
      # Plot method for this row of the figure.
      method <- plotMethods[pm]
      
      # Prepare true component functions values for this species.
      if ( includeTrueFunctions ) {
        # Just this species' functions, not its intercept (i.e. alpha_k).
        trueComp <- cellsObj$trueComp[ ,k,2:(nEnvirs+1)]
        
        # Add sampling bias functions (without species specific intercepts), if required.
        if ( includeBiasCovars ) 
          trueComp <- cbind(trueComp, cellsObj$trueBComp[ ,(nSpecies+1):(nSpecies+nBiases)])
      }
      
      # Plot these columns for the given species x method combination
      for ( cv in 1:nCovars ) {
        xlimRange <- range(covars[ ,cv])
        nameCovar <- namesCovars[cv]
        skipPlot <- FALSE
        
        # Which method has been requested?
        xRug <- NULL
        if ( method == "MGCV" ) {
          if ( useDataMgcv %in% c("PA","AB") ) {
            if ( cv <= nEnvirs ) {
              # If want rug plot, get covariate values for sample locations only!
              xRug <- data.MGCV[ ,nameCovar]
              if ( includeTrueFunctions ) {
                rows <- surveysObj$rowsInCells
                yConstraint <- mean(trueComp[rows,cv])
              }
              
              # GAM on survey data.
              gamTmp <- gamVals[[k]]
              x <- gamTmp[[cv]]$x
              y <- gamTmp[[cv]]$fit
              yplusSE <- y+gamTmp[[cv]]$se
              yminusSE <- y-gamTmp[[cv]]$se
            } else {
              skipPlot <- TRUE
            }
              
          } else if ( useDataMgcv == "PO" ) {
            # GAM on PO data.
            gamTmp <- gamVals[[k]]

            # If want rug plot, get covariate values for all cells?
            xRug <- data.MGCV[ ,nameCovar]
            if ( includeTrueFunctions ) {
              yConstraint <- mean(trueComp[ ,cv])
            }
            
            if ( is.null(corBias) ) {
              # No correlation on sampling bias
              x <- gamTmp[[cv]]$x
              y <- gamTmp[[cv]]$fit
              yplusSE <- y+gamTmp[[cv]]$se
              yminusSE <- y-gamTmp[[cv]]$se
              
            } else {
              # Correlation on sampling bais
              posCorBias <- nEnvirs + corBias[2]
              if ( cv < posCorBias ) {
                # No change in position.
                x <- gamTmp[[cv]]$x
                y <- gamTmp[[cv]]$fit
                yplusSE <- y+gamTmp[[cv]]$se
                yminusSE <- y-gamTmp[[cv]]$se
                
              } else if ( cv == posCorBias ) {
                skipPlot <- TRUE  
                
              } else {
                # Change in position.
                x <- gamTmp[[cv-1]]$x
                y <- gamTmp[[cv-1]]$fit
                yplusSE <- y+gamTmp[[cv-1]]$se
                yminusSE <- y-gamTmp[[cv-1]]$se
              }
            }
          } else {
            stop("Unrecognised type of data used in MGCV method.")
          }
          
        } else if ( method == "JAGS" ) {
          # TO DO: work out how to do this for bias covars.
          # If want rug plot, get covariate values for sample locations only!
          xRug <- covars[ ,cv]
          if ( includeTrueFunctions ) {
            rows <- surveysObj$rowsInCells
            yConstraint <- mean(c(trueComp[rows,cv],trueComp[ ,cv]))
          }
          
          # popaJags on PAPO data.
          # TO DO: work out how to do this for bias covars.
          jagsTmp <- jagsVals[[k]]
          x <- jagsTmp[[cv]]$x
          y <- jagsTmp[[cv]]$fit
          yplusSE <- y+jagsTmp[[cv]]$se
          yminusSE <- y-jagsTmp[[cv]]$se
          
        } else if ( method == "TMB" ) {
          # If want rug plot, get covariate values for sample locations only!
          xRug <- covars[ ,cv]
          if ( includeTrueFunctions ) {
            rows <- surveysObj$rowsInCells
            yConstraint <- mean(c(trueComp[rows,cv],trueComp[ ,cv]))
          }
          
          # popaTMB on PAPO data.
          tmbTmp <- tmbVals[[k]]
          x <- tmbTmp[[cv]]$x
          y <- tmbTmp[[cv]]$fit
          yplusSE <- y+tmbTmp[[cv]]$se
          yminusSE <- y-tmbTmp[[cv]]$se
        } else {
          stop(paste0("Method '", method, "' does not have plotting code."))
        }
        
        # Plot ...
        if ( ! skipPlot ) {
          ytrue <- trueComp[ ,cv] - yConstraint
          ylimRange <- range(c(ytrue,y,yplusSE,yminusSE))
          plot(x, y, type="l", ann=FALSE, xlim=xlimRange, ylim=ylimRange, xaxt="n")
          axis(1, at=NULL, labels=FALSE, tick=TRUE)    # tick marks only!
          
          # Plot identifier labels in top right corner of each plot, if requested.
          if ( doPlotIdentifiers ) {
            iPlotCount <- iPlotCount + 1
            if ( posIdentifiers == "topleft" ) {
              posId <- c(xlimRange[1], ylimRange[2])
              adjId <- c(0,1)
            } else if ( posIdentifiers == "topright") {
              posId <- c(xlimRange[2], ylimRange[2])
              adjId <- c(1,1)
            } else if ( posIdentifiers == "bottomleft" ) {
              posId <- c(xlimRange[1], ylimRange[1])
              adjId <- c(0,0)
            } else if ( posIdentifiers == "bottomright" ) {
              posId <- c(xlimRange[2], ylimRange[1])
              adjId <- c(1,0)
            } else {
              stop("This value of posIdentifers is not valid.")
            }
            text(x=posId[1], y=posId[2], label=plotIdentifier[iPlotCount], adj=adjId, cex=0.8)
          }
          
          # Method label on first column y-axis.
          if ( cv == 1 ) {
            midyRange <- ((ylimRange[2] - ylimRange[1])/2.0) + ylimRange[1]
            ylabel=substitute(hat(f)[i*sp] * " from " * m, list(sp=k, m=method))
            axis(2, at=midyRange, labels=ylabel, tick=FALSE, padj=-1.5)
          }
          
          # Add confidence intervals to plot
          lines(x, yplusSE, lty="dashed")
          lines(x, yminusSE, lty="dashed")
          
          # Add rug of sample or observation covariate values used to form estimate.
          if ( includeRug ) rug(xRug, ticksize=0.03)
          
          # Add true values for comparison to plot
          if ( includeTrueFunctions ) {
            xtrue <- covars[ ,cv]
            indOrdered <- sort(xtrue,index.return=TRUE)$ix
            lines(xtrue[indOrdered], ytrue[indOrdered], lty="solid", col=trueCol)
          }

          # Add x-axis labels to last row only.
          if ( k == nSpecies && pm == nMethods ) {
            axis(1, at=NULL, labels=TRUE, tick=FALSE)   # tick mark labels
            midxRange <- ((xlimRange[2] - xlimRange[1])/2.0) + xlimRange[1]
            xlabel <- namesCovars[cv]
            axis(1, at=midxRange, labels=xlabel, tick=FALSE, padj=2.0) # x-axis labels
          }
        } else {
          plot.new()
        }      
      } # for ( cv in 1:nCovars )   ...  i.e. columns of plots
    }  # for ( pm in 1:nMethods ) ...  i.e. inner rows of plots if plotting to file
    
    # Title for page.
    if ( ! plotToFile) title(main="Estimated component functions for ",namesSpecies[k], outer=TRUE)
    
  }  # for ( k in 1:nSpecies )  ...  i.e. outer rows of plots if plotting to file
  
  # Turn off the plotting device, if it is a file.
  par(mfrow = c(1,1), oma = oma.old, mar = mar.old)
  if ( plotToFile ) dev.off()
  
}

#-----------------------------------------------------------------------------------------

plotSimilarity <-function(expDir, 
                          expDataFiles,
                          numSpecies,
                          plotColTitles = paste("Exp",1:length(expDataFiles)),
                          plotDevice = "RStudioGD",
                          plotDir = expDir,
                          plotUnits = "cm",
                          plotWidth = 12,
                          plotHeight = 10,
                          plotRes = 600,
                          plotFileName = "similarityBoxPlot",
                          doPlotIdentifiers = TRUE,
                          posIdentifiers = "bottomright", 
                          quiet = TRUE) {
    
  # Function to box plot the measure of similarity in the different experiments contained
  # in the given data files.  Assumes the same number of simulations have been performed 
  # in each experiment.
  
  # Check there are the same number of data files as there are methods to plot.
  nExps <- length(expDataFiles)
  if ( length(plotColTitles) != nExps ) 
    stop("The number of plot column titles does not match the number of experiments.")
  
  # Plotting initialisation ...  
  if ( plotDevice != "RStudioGD" ) {
    plotFile <- makeFileName(plotFileName, expDir, plotDevice)
    argList <- list(filename = plotFile, width = plotWidth, height = plotHeight, 
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
  }
  opar <- par(mfcol=c(nSpecies,nExps), oma=c(2,2,2,0.1), mar=c(1.5,1.5,0,0), cex=0.66)
  ylimRange <- c(0.0,1.0)   # similarity measure is correlation so know this.

  # Run loops to produce Box plot sub-plots.
  plotIdentifiers <- matrix(letters[1:(nExps*nSpecies)], nrow=nSpecies, ncol=nExps, byrow = TRUE)
  for ( i in 1:nExps ) {
    # Load similarity data
    loadFile <- makeFileName(expDataFiles[i], expDir, "RData")  
    load(loadFile)     # load in corLambda and corEta  (FYI: eta = ln(lambda)).
    
    # Initialise bits and pieces
    nSims <- dim(corLambda)[1]
    useMethods <- dimnames(corLambda)[[2]]
    nMethods <- dim(corLambda)[2]
    nSpecies <- dim(corLambda)[3] 
    if ( nSpecies != numSpecies ) 
      warning("The number of species in the data does not match number specified for plotting.")
    xlimRange <- c(0.5,nMethods+0.5)
    
    # Plot for each species (species are rows of figure)
    for ( k in 1:nSpecies ) {
      dfLambda <- as.data.frame(corLambda[ , ,k])
      if ( ! quiet ) {
        # Sims with any non-converged methods.
        for ( m in 1:nMethods) {
          nSimsNA <- sum( is.na(corLambda[ ,m,k]) )
          if ( nSimsNA > 0 )  
            message( paste0("Species ",k," has ", nSimsNA, 
                            " simulations that didn't converge for method ", useMethods[m],
                            " in experiment ", i, ".") )
        }
      }
      
      # Box plot similarity measures for each method.
      boxplot(dfLambda, ylim=ylimRange, xaxt="n", yaxt="n", xlab=NULL, ylab=NULL)
      abline(h=1.0,lty="dotted")
      axis(1, at=c(1:nMethods), labels=FALSE, tick=TRUE)    # tick marks only!
      axis(2, at=NULL, labels=FALSE, tick=TRUE)    # tick marks only!
      
      # Title and x-axis labels
      if ( k == 1 ) {
        # First row of plot.
        colTitle <- plotColTitles[[i]]
        midxRange <- ((nMethods - 1.0)/2.0) + 1.0
        axis(3, at=midxRange, labels=colTitle, tick=FALSE, padj=0.5) # axis label
      } else if ( k == nSpecies ) {
        # Last row of plot.
        axis(1, at=c(1:nMethods), labels = useMethods, tick=FALSE) # tick mark labels
      } else {
        # nothing at this stage.
      }
      
      # Y-axis labels
      if ( i == 1 ) {
        # First column of plot.
        midyRange <- ((ylimRange[2] - ylimRange[1])/2.0) + ylimRange[1]
        ylabel=substitute("similarity of " * hat(mu)[c*sp] * " to " * mu[c*sp], list(sp=k))
        #ylabel=paste("similarity of estimated intensities, k =", k)
        axis(2, at=midyRange, labels=ylabel, tick=FALSE, padj=-1.1) # axis label
        axis(2, labels=TRUE, tick=FALSE) # tick mark labels
      }
      
      # Plot identifier labels in top right corner of each plot, if requested.
      if ( doPlotIdentifiers ) {
        if ( posIdentifiers == "topleft" ) {
          posId <- c(xlimRange[1], ylimRange[2])
          adjId <- c(0,1)
        } else if ( posIdentifiers == "topright") {
          posId <- c(xlimRange[2], ylimRange[2])
          adjId <- c(1,1)
        } else if ( posIdentifiers == "bottomleft" ) {
          posId <- c(xlimRange[1], ylimRange[1])
          adjId <- c(0,0)
        } else if ( posIdentifiers == "bottomright" ) {
          posId <- c(xlimRange[2], ylimRange[1])
          adjId <- c(1,0)
        } else {
          stop("This value of posIdentifers is not valid.")
        }
        text(x=posId[1], y=posId[2], label=paste0("(", plotIdentifiers[k,i], ")"), 
             adj=adjId, cex=0.8)
      }
      
    } 
  }

  # Reset plotting device.
  par(opar)
  if ( plotDevice != "RStudioGD" ) dev.off()
  
}

#-----------------------------------------------------------------------------------------
