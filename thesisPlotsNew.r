#-----------------------------------------------------------------------------------------
#
# Test example with various curves for environment and sampling bias covariates in
# Fithian model.  Solve using a GAM like popa function (either popaTMB or popaJAGS).
# Compare results with GAM solution using data (i.e. samples recording counts).
#
#-----------------------------------------------------------------------------------------

library(popaGAM)
library(fields,quietly = TRUE)
library(raster)
library(spatstat)
library(mgcv)
library(rjags)
library(TMB)

mySourceDir <- "/perm_storage/home/sampeel/chap3/sim2/samR/"
# HOME: "C:/Users/sampe/OneDrive - University of Tasmania/Thesis/Paper 3/Spoon/popaGAM2/samR/"
# UNI:  "C:/Users/slpeel/OneDrive - University of Tasmania/Thesis/Paper 3/Spoon/popaGAM2/samR/"
# VM:   "/perm_storage/home/sampeel/chap3/sim2/samR/"
source(paste0(mySourceDir,"cellsObj.r"))
source(paste0(mySourceDir,"domainObj.r"))
source(paste0(mySourceDir,"surveysObj.r"))
source(paste0(mySourceDir,"utils.r"))
source(paste0(mySourceDir,"popaJagam3.r"))
source(paste0(mySourceDir,"simFuncs.r"))

# Experiment set up.
lstObjs <- setupSim(nRows = 60, nCols = 80, nEnvirs = 3, nBiases = 2, nSpecies = 2)

# Sample bias coefficients
gamma <- c(-4.9,-5.5)                 # Sample bias intercepts 
delta <- c(0,-0.1,-0.4,0.2)           # Coefficient order: z_1, z_2, z_1^2, z_2^2, ...

#-----------------------------------------------------------------------------------------
#
# Quadratic component function example
#
#-----------------------------------------------------------------------------------------

# Species distribution coefficients
cellsObj <- lstObjs$cells
alpha <- c(4.3, 4.1)                  # overall intercept per species
beta <- array(dim = c(3, cellsObj$numCovars, cellsObj$numSpecies), 
              dimnames=list(c("b0","b1","b2"), cellsObj$namesCovars, cellsObj$namesSpecies))
beta[ ,1,1] <- c(0.0,-0.2,-0.4)       # covariate component per species ("b0" is component intercept)
beta[ ,2,1] <- c(0.0,-0.1, 0.1)
beta[ ,3,1] <- c(0.0, 0.4,-0.2)
beta[ ,1,2] <- c(0.0,-0.1,-0.4)
beta[ ,2,2] <- c(0.0,-0.5, 0.0)
beta[ ,3,2] <- c(0.0,-0.4,-0.2)


# Figure 4.1 - Example 1 for high numbers of PA data ...
nSamps <- 2000
expType <- "quad"
useDataMgcv <- "PA"
corBias <- NULL
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
useMethods <- c("MGCV","TMB")
for (simNum in 1:1) { 
  # Experiment run.
  lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                          useMethods = useMethods, useDataMgcv = useDataMgcv, 
                          alpha = alpha, beta = beta, sim = simNum-1, saveDir = expDir,
                          quiet = TRUE, plotCheck = FALSE, plotResults = FALSE)
  
  
  # Plot results
  expDataFiles <- useMethods
  for ( m in 1:length(useMethods)) {
    expDataFiles[m] <- paste0(tolower(useMethods[m]),"Results-n",useDataMgcv,nSamps,"-S",simNum)
  }
  plotFileName <- paste0("estCompFuncs-n",useDataMgcv,nSamps, "-nocor-S",simNum)
  plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", plotHeight=15,
                     plotWidth=12, plotFileName = plotFileName, 
                     includeBiasCovars = FALSE, posIdentifiers = "topright")
}                   

# Figure 4.3 - Example 1 for low numbers of PA data ...
nSamps <- 200
expType <- "quad"
useDataMgcv <- "PA"
corBias <- NULL
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
useMethods <- c("MGCV","TMB")
for (simNum in 1:1) {
  # Experiment run.
  lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                          useMethods = useMethods, useDataMgcv = useDataMgcv, 
                          alpha = alpha, beta = beta, sim = simNum-1, saveDir = expDir,
                          quiet = TRUE, plotResults = FALSE)
  
  # Quadratic curves, no correlation, nPA = 200.
  expDataFiles <- useMethods
  for ( m in 1:length(useMethods)) {
    expDataFiles[m] <- paste0(tolower(useMethods[m]),"Results-n",useDataMgcv,nSamps,"-S",simNum)
  }
  plotFileName <- paste0("estCompFuncs-n",useDataMgcv,nSamps, "-nocor-S",simNum)
  plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", plotHeight=15,
                     plotWidth=12, plotFileName = plotFileName, 
                     includeBiasCovars = FALSE)
  
}


# Figure 4.? - Simulations to look at generality of result.
start.time <- Sys.time()
nSamps <- 200  # 2000
expType <- "quad"
useDataMgcv <- "PA"
corBias <- NULL
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
useMethods <- c("MGCV","TMB")
nMethods <- length(useMethods)
nSims <- 100
nSpecies <- lstObjs$cells$numSpecies
corEta <- array(dim=c(nSims, nMethods, nSpecies),
                dimnames = list(1:nSims, useMethods, lstObjs$cells$namesSpecies))
corLambda <- corEta
doCompareSingleSimPlot <- FALSE
for (simNum in 1:nSims) {
  # Experiment run: Mixed curves, no correlation, nPA = 200.
  start.time.sim <- Sys.time()
  message("Simulation number: ", simNum)
  lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                          useMethods, alpha = alpha, beta = beta, gamma=gamma, delta=delta, 
                          useDataMgcv = useDataMgcv, sim = simNum-1, 
                          saveDir = NULL, quiet = TRUE, plotCheck =  FALSE, plotResults = FALSE)
  sim.time <- Sys.time() - start.time.sim
  message("Time taken: ", sim.time)
  
  # Calculate correlation to true ln(intensity).
  for ( k in 1:nSpecies ) {
    trueLambda <- lstRes$cellsObj$trueLambda[ ,k]
    trueEta <- log(trueLambda)
    for ( m in 1:nMethods ) {
      noMethod <- FALSE
      if (  useMethods[m] == "MGCV" ) {
        if ( ! is.null(lstRes$lstFit.MGVC) ) 
          estEta <- mgcv::predict.gam(lstRes$lstFit.MGVC[[k]], lstRes$cellsObj$covars, 
                                      type="link")
      } else if ( useMethods[m] == "TMB" ) {
        if ( ! is.null(lstRes$fit.TMB) )
          estEta <- popaGAM:::predict.popaGAM(lstRes$fit.TMB, lstRes$cellsObj$covars, 
                                              "intensity", k, type="link")
      } else if ( useMethods[m] == "JAGS" ) {
        if ( ! is.null(lstRes$fit.JAGS) )
          estEta <- popaGAM:::predict.popaGAM(lstRes$fit.JAGS, lstRes$cellsObj$covars, 
                                              "intensity", k, type="link")
      } else {
        warning(paste("Method", useMethods[m], "has not been predicted."))
        noMethod <- TRUE
      }
      if ( ! noMethod ) {
        corEta[simNum, useMethods[m], k] <- cor(trueEta,estEta)
        estLambda <- exp(estEta)
        corLambda[simNum, useMethods[m], k] <- cor(trueLambda, estLambda)
      }
    }
  }
  
  # plot comparison between methods.
  if ( doCompareSingleSimPlot ) {
    expDataFiles <- useMethods
    for ( m in 1:length(useMethods)) {
      expDataFiles[m] <- paste0(tolower(useMethods[m]),"Results-n",useDataMgcv,nSamps,"-S",simNum)
    }
    plotFileName <- paste0("estCompFuncs-n",useDataMgcv,nSamps, "-nocor-S",simNum)
    plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", includeBiasCovars = FALSE,
                       plotWidth=12, plotHeight = 15, plotFileName = plotFileName)
  }
}
saveFile <- makeFileName(paste0("corRes-n",useDataMgcv,nSamps,"-nSim",nSims),expDir,"RData")
save(corEta, corLambda, file=saveFile)
Total.time <- difftime(Sys.time(), start.time) 
Total.time # message("Total time taken: ", Total.time)

# Box plots of similarity measures (i.e. correlations).
expType <- "quad"
useDataMgcv <- "PA"
nSims <- 100
nSpecies <- 2
nSampsInExps <- c(2000,200)
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
nExps <- length(nSampsInExps)
expDataFiles <- vector("character",nExps)
plotColTitles <- list()
for ( i in 1:nExps  ) {
  nSamp <- nSampsInExps[i]
  expDataFiles[i] <- (paste0("corRes-n",useDataMgcv,nSamp,"-nSim",nSims))
  plotColTitles[[i]] <- substitute(n[samp] * " = " * ns, list(ns=nSamp))
}

plotFile <- paste0("similarityPlot-n",useDataMgcv,"-nSim",nSims)
plotSimilarity(expDir, expDataFiles, nSpecies, plotColTitles, plotDevice = "png",
               plotWidth = 12,  plotHeight = 10, plotFileName = plotFile, quiet = FALSE)



#-----------------------------------------------------------------------------------------
#
# Mixed component function example - no correlation
#
#-----------------------------------------------------------------------------------------

# Figure 4.2 & 4.4 - Example 2 for high or lower numbers of PA data ...
nSamps <- 2000  # 200
expType <- "mixed"
useDataMgcv <- "PA"
corBias <- NULL
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
useMethods <- c("MGCV","TMB")
for (simNum in 1:1) {
  # Experiment run: Mixed curves, no correlation, nPA = 2000.
  lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                          useMethods, useDataMgcv = useDataMgcv, sim = simNum-1, 
                          saveDir = expDir, quiet = TRUE, plotResults = FALSE)

  # plot comparison between methods.
  expDataFiles <- useMethods
  for ( m in 1:length(useMethods)) {
    expDataFiles[m] <- paste0(tolower(useMethods[m]),"Results-n",useDataMgcv,nSamps,"-S",simNum)
  }
  plotFileName <- paste0("estCompFuncs-n",useDataMgcv,nSamps, "-nocor-S",simNum)
  plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", includeBiasCovars = FALSE,
                     plotWidth=12, plotHeight = 15, plotFileName = plotFileName)
  
  # Experiment run: Mixed curves, no correlation, nPA = 2000, mgcv uses PO.
  # (FYI: TMB and JAGS will produce same results as above experiment.)
  # lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias,
  #                         useMethods = c("MGCV"), useDataMgcv = "PO", sim = simNum-1, 
  #                         saveDir = expDir, quiet = TRUE, plotResults = FALSE)
  # 
  # expDataFiles <- useMethods
  # useDataMgcvTmp <- c("PO",repuseDataMgcv)
  # nData <- c(lstObjs$cells$numCells,nSamps)
  # for ( m in 1:length(useMethods)) {
  #   expDataFiles[m] <- paste0(tolower(useMethods[m]), "Results-n", useDataMgcvTmp[m],
  #                             nData[m],"-S", simNum)
  # }
  # plotFileName <- paste0("estCompFuncs-n",useDataMgcv, nSamps, "nPO-nocor-S",simNum)
  # plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", includeBiasCovars = TRUE,
  #                    plotWidth=12, plotHeight = 15, plotFileName = plotFileName )
}


# Figure 4.5 - Simulations to look at generality of result.
start.time <- Sys.time()
nSamps <- 2000  # 200
expType <- "mixed"
useDataMgcv <- "PA"
corBias <- NULL
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
useMethods <- c("MGCV","TMB")
nMethods <- length(useMethods)
nSims <- 100
nSpecies <- lstObjs$cells$numSpecies
corEta <- array(dim=c(nSims, nMethods, nSpecies),
                   dimnames = list(1:nSims, useMethods, lstObjs$cells$namesSpecies))
corLambda <- corEta
doCompareSingleSimPlot <- FALSE
for (simNum in 1:nSims) {
  # Experiment run: Mixed curves, no correlation, nPA = 200.
  start.time.sim <- Sys.time()
  message("Simulation number: ", simNum)
  lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                          useMethods, gamma=gamma, delta=delta, useDataMgcv = useDataMgcv, sim = simNum-1, 
                          saveDir = NULL, quiet = TRUE, plotCheck =  FALSE, plotResults = FALSE)
  sim.time <- Sys.time() - start.time.sim
  message("Time taken: ", sim.time)
  
  # Calculate correlation to true ln(intensity).
  for ( k in 1:nSpecies ) {
    trueLambda <- lstRes$cellsObj$trueLambda[ ,k]
    trueEta <- log(trueLambda)
    for ( m in 1:nMethods ) {
      noMethod <- FALSE
      if (  useMethods[m] == "MGCV" ) {
        if ( ! is.null(lstRes$lstFit.MGVC) ) 
          estEta <- mgcv::predict.gam(lstRes$lstFit.MGVC[[k]], lstRes$cellsObj$covars, 
                                    type="link")
      } else if ( useMethods[m] == "TMB" ) {
        if ( ! is.null(lstRes$fit.TMB) )
          estEta <- popaGAM:::predict.popaGAM(lstRes$fit.TMB, lstRes$cellsObj$covars, 
                                            "intensity", k, type="link")
      } else if ( useMethods[m] == "JAGS" ) {
        if ( ! is.null(lstRes$fit.JAGS) )
          estEta <- popaGAM:::predict.popaGAM(lstRes$fit.JAGS, lstRes$cellsObj$covars, 
                                              "intensity", k, type="link")
      } else {
        warning(paste("Method", useMethods[m], "has not been predicted."))
        noMethod <- TRUE
      }
      if ( ! noMethod ) {
        corEta[simNum, useMethods[m], k] <- cor(trueEta,estEta)
        estLambda <- exp(estEta)
        corLambda[simNum, useMethods[m], k] <- cor(trueLambda, estLambda)
      }
    }
  }
  
  # plot comparison between methods.
  if ( doCompareSingleSimPlot ) {
    expDataFiles <- useMethods
    for ( m in 1:length(useMethods)) {
      expDataFiles[m] <- paste0(tolower(useMethods[m]),"Results-n",useDataMgcv,nSamps,"-S",simNum)
    }
    plotFileName <- paste0("estCompFuncs-n",useDataMgcv,nSamps, "-nocor-S",simNum)
    plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", includeBiasCovars = FALSE,
                       plotWidth=12, plotHeight = 15, plotFileName = plotFileName)
  }
}
saveFile <- makeFileName(paste0("corRes-n",useDataMgcv,nSamps,"-nSim",nSims),expDir,"RData")
save(corEta, corLambda, file=saveFile)
Total.time <- difftime(Sys.time(), start.time) 
message("Total time taken: ", Total.time)

# Box plots of similarity measures (i.e. correlations).
expType <- "mixed"
useDataMgcv <- "PA"
nSims <- 100
nSpecies <- 2
nSampsInExps <- c(2000,200)
expDir <- paste0(getwd(),"/Example-", expType,"-nocor")
nExps <- length(nSampsInExps)
expDataFiles <- vector("character",nExps)
plotColTitles <- list()
for ( i in 1:nExps  ) {
  nSamp <- nSampsInExps[i]
  expDataFiles[i] <- (paste0("corRes-n",useDataMgcv,nSamp,"-nSim",nSims))
  plotColTitles[[i]] <- substitute(n[samp] * " = " * ns, list(ns=nSamp))
}

plotFile <- paste0("similarityPlot-n",useDataMgcv,"-nSim",nSims)
plotSimilarity(expDir, expDataFiles, nSpecies, plotColTitles, plotDevice = "png",
               plotWidth = 12,  plotHeight = 10, plotFileName = plotFile)




#-----------------------------------------------------------------------------------------
#
# Mixed component function example - correlation of sample bias.
#
#-----------------------------------------------------------------------------------------

# Experiment run.
nSamps <- 2000
expType <- "mixed"
useDataMgcv <- "PA"
corBias <- c(1,1)
lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                        useMethods = c("MGCV","TMB"), useDataMgcv = useDataMgcv, 
                        quiet = TRUE, plotResults = FALSE)

# Experiment run (FYI: TMB and JAGS will produce same results as above experiment.)
nSamps <- 2000
expType <- "mixed"
useDataMgcv <- "PO"
corBias <- c(1,1)
lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias,
                        useMethods = c("MGCV"), useDataMgcv = useDataMgcv, 
                        quiet = TRUE, plotResults = FALSE)


# Mixed curves, correlation, nPA = 2000, mgcv uses PO.
expDir <- paste0(getwd(),"/Example-mixed-z1cor")
expDataFiles <- c("mgcvResults-nPO4800", 
                  # "jagsResults-nPA200", 
                  "tmbResults-nPA2000")
plotMethods <- c("MGCV","TMB") #c("MGCV","JAGS","TMB")
plotCompareMethods(expDir, expDataFiles, plotMethods, plotDevice="png", includeBiasCovars = TRUE,
                   plotWidth=12, plotHeight = 15, plotFileName = "estCompFuncs-nPA2000-nPO-V2")

# Experiment run.
nSamps <- 200
expType <- "mixed"
useDataMgcv <- "PA"
corBias <- c(1,1)
lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                        useMethods = c("MGCV","TMB"), useDataMgcv = useDataMgcv, 
                        quiet = TRUE, plotResults = FALSE)

# Mixed curves, correlation, nPA = 200, mgcv uses PO.
expDir <- paste0(getwd(),"/Example-mixed-z1cor")
expDataFiles <- c("mgcvResults-nPO4800", 
                  # "jagsResults-nPA200", 
                  "tmbResults-nPA200")
plotMethods <- c("MGCV","TMB") #c("MGCV","JAGS","TMB")
plotCompareMethods(expDir, expDataFiles, plotMethods, plotDevice="png", includeBiasCovars = TRUE,
                   plotWidth=12, plotHeight = 15, plotFileName = "estCompFuncs-nPA200-nPO-V2")


#-----------------------------------------------------------------------------------------
#
# Figure 4.7 - Mapped plots for correlated sampling bias
#
#-----------------------------------------------------------------------------------------


# What data (from above experiment runs) ...
expDir <- paste0(getwd(),"/Example-mixed-z1cor")
expDataFiles <- c("mgcvResults-nPO4800",
                  "mgcvResults-nPA2000",
                  "tmbResults-nPA2000")

# Plot settings.
plotToFile <- TRUE
plotDevice <- "png"
plotUnits <- "cm"
plotWidth <- 15.5
plotHeight <- 11.2
plotRes <- 600
mar.old <- par("mar")
mar.new <- c(2.1,2.1,1.3,1.1)
txtSize <- 0.66
lgdSize <- 0.66
posTitle <- c(3,57)
titleSize <- txtSize
plotNames <- c("True",  "MGCV (PO)", "MGCV (PA)", "TMB (POPA)")
speciesNames <- c("sp1","sp2")
corRes <- matrix(nrow=length(plotNames),ncol=length(speciesNames),
                 dimnames = list(plotNames,speciesNames))
legend.text <- expression(phantom(0)*phantom(0)*phantom(0)*ln(mu[c2]))  ## Hardwired "centering" of title


# For each species ...
for ( k in 1:length(speciesNames) ) {
  # Plot to file?
  if ( plotToFile ) {
    plotFile <- makeFileName(paste0("ResultsMapV3-sp",k), expDir, plotDevice)
    argList <- list(filename = plotFile, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
  }

  # Set up page for plotting.
  par(mfrow=c(2,2), mar=mar.new, cex=txtSize)
  
  # Plot true response.
  load(makeFileName(expDataFiles[1], expDir, "RData"))  # To get a cell object for true results.
  trueEta <- log(cellsObj$trueLambda[ ,k])
  rlPlot <- rasterFromXYZ(cbind(cellsObj$xy, trueEta), res = cellsObj$resCell)
  plot(rlPlot, asp=1)
#  plot(rlPlot, asp=1, legend=FALSE)
#  plot(rlPlot, asp=1, legend.only=TRUE, legend.args=list(text=legend.text, side=3, cex=lgdSize))  
  title(expression("True value of " * ln(mu[c2])), outer = FALSE, cex.main=titleSize)
  corRes[1,k] <- cor(trueEta,trueEta)
  
  # Plot mgcv estimated response from PO data (FYI: file loaded above).
  covars <- cbind(cellsObj$covars, cellsObj$biases)
  names(covars) <- c(names(cellsObj$covars), names(cellsObj$biases))
  estLambda <- mgcv::predict.gam(lstFit.MGCV[[k]], covars, type="link")
  par(cex=txtSize)
  rlPlot <- rasterFromXYZ(cbind(cellsObj$xy, estLambda), res = cellsObj$resCell)
  plot(rlPlot, asp=1)
  title(expression("MGCV (PO) estimate of " * ln(mu[c2])), outer = FALSE, cex.main=titleSize)
  corRes[2,k] <- cor(trueEta, estLambda)
  
  # Plot mgcv estimated response from PA data.
  load(makeFileName(expDataFiles[2], expDir, "RData"))
  estLambda <- mgcv::predict.gam(lstFit.MGCV[[k]], cellsObj$covars, type="link")
  par(cex=txtSize)
  rlPlot <- rasterFromXYZ(cbind(cellsObj$xy, estLambda), res = cellsObj$resCell)
  plot(rlPlot, asp=1)
  title(expression("MGCV (PA) estimate of " * ln(mu[c2])), outer = FALSE, cex.main=titleSize)
  corRes[3,k] <- cor(trueEta, estLambda)
  
  # # Plot JAGS response.
  # #load(makeFileName(expDataFiles[3], expDir, "RData"))
  # rlPlot[] <- NA
  # plot(rlPlot,asp=1)
  # text(posTitle[1],posTitle[2],"JAGS (POPA)", adj=c(0,1), cex=titleSize)
  # 
  # Plot TMB estimated response from PO and PA data.
  load(makeFileName(expDataFiles[3], expDir, "RData"))
  estLambda <- popaGAM:::predict.popaGAM(fit.TMB, cellsObj$covars, "intensity", k, type="link")
  par(cex=txtSize)
  rlPlot <- rasterFromXYZ(cbind(cellsObj$xy, estLambda), res = cellsObj$resCell)
  plot(rlPlot, asp=1)
  title(expression("TMB (POPA) estimate of " * ln(mu[c2])), outer = FALSE, cex.main=titleSize)
  corRes[4,k] <- cor(trueEta, estLambda)
  
  if ( plotToFile ) dev.off()
}
par(mfrow=c(1,1), mar=mar.old, cex=1.0)

#-----------------------------------------------------------------------------------------
#
# Figure 4.8 - Multiple simulations and correlation metric for sample bias example
#
#-----------------------------------------------------------------------------------------

# Experiment run for MGCV(PA) and TMB (POPA)
start.time <- Sys.time()
nSamps <- 2000  
expType <- "mixed"
corBias <- c(1,1)
expDir <- paste0(getwd(),"/Example-", expType,"-z1cor")
useMethods <- c("MGCV","MGCV","TMB")
useDataMgcv <- c("PA","PO","POPA")
nMethods <- length(useMethods)
nSims <- 100
nSpecies <- lstObjs$cells$numSpecies
expNames <- paste0(useMethods," (",useDataMgcv,")")
corEta <- array(dim=c(nSims, nMethods, nSpecies),
                dimnames = list(1:nSims, expNames, lstObjs$cells$namesSpecies))
corLambda <- corEta
doCompareSingleSimPlot <- FALSE
for (simNum in 1:nSims) {
  # Experiment run: Mixed curves, no correlation, nPA = 200.
  start.time.sim <- Sys.time()
  message("Simulation number: ", simNum)
  lstRes <- runExperiment(lstObjs$cells, lstObjs$domain, nSamps, expType, corBias, 
                          useMethods, gamma=gamma, delta=delta, useDataMgcv = useDataMgcv, sim = simNum-1, 
                          saveDir = NULL, quiet = TRUE, plotCheck =  FALSE, plotResults = FALSE)
  sim.time <- Sys.time() - start.time.sim
  message("Time taken: ", sim.time)
  
  # Calculate correlation to true ln(intensity).
  for ( k in 1:nSpecies ) {
    trueLambda <- lstRes$cellsObj$trueLambda[ ,k]
    trueEta <- log(trueLambda)
    
    for ( m in 1:nMethods ) {
      # Initialise inner loop.
      noMethod <- FALSE
      nameFit <- paste0("fit.", useMethods[m])
      if (  useMethods[m] == "MGCV" ) {
        if ( length(useDataMgcv) == 1 ) {
          nameFit <- paste0(nameFit, ".",useDataMgcv)
        } else {
          nameFit <- paste0(nameFit, ".",useDataMgcv[m])
        }
      }
      nameCor <- paste0(useMethods[m]," (",useDataMgcv[m],")")
      
      # Did this method converge?
      df <- cbind(lstRes$cellsObj$covars,lstRes$cellsObj$biases)
      
      if ( ! is.null(lstRes[[nameFit]]) ) {
        if ( useMethods[m] == "MGCV") {
          estEta <- mgcv::predict.gam(lstRes[[nameFit]][[k]], newdata = df, type="link")
        } else if ( ( useMethods[m] == "TMB" ) || ( useMethods[m] == "JAGS" ) ){
          estEta <- popaGAM:::predict.popaGAM(lstRes[[nameFit]], newdata = df, 
                                              "intensity", k, type="link")
        } else {
            warning(paste("Method", useMethods[m], "has no correlation to true code."))
            noMethod <- TRUE
        }
      
        if ( ! noMethod ) {
          corEta[simNum, nameCor, k] <- cor(trueEta,estEta)
          estLambda <- exp(estEta)
          corLambda[simNum, nameCor, k] <- cor(trueLambda, estLambda)
        }
      } else {
        # The method did not converge?
      }
    }
  }
  
  # plot comparison between methods.
  if ( doCompareSingleSimPlot ) {
    expDataFiles <- useMethods
    for ( m in 1:length(useMethods)) {
      if ( useMethods[m] == "MGCV") {
        useData <- useDataMgcv[m]
        numData <- ifelse(useDataMgcv[m]=="PO",cellsObj$numCells,nSamps)
      } else {
        useData <- "PA"
        numData <- nSamps
      }
      expDataFiles[m] <- paste0(tolower(useMethods[m]),"Results-n",useData,numData,"-S",simNum)
    }
    plotFileName <- paste0("estCompFuncs-n",useData,nSamps, "-z1cor-S",simNum)
    plotCompareMethods(expDir, expDataFiles, useMethods, plotDevice="png", includeBiasCovars = FALSE,
                       plotWidth=12, plotHeight = 15, plotFileName = plotFileName)
  }
}
saveFile <- makeFileName(paste0("corRes-sampBias-nSim",nSims),expDir,"RData")
save(corEta, corLambda, file=saveFile)
Total.time <- difftime(Sys.time(), start.time) 
message("Total time taken: ", Total.time)

# Box plots of similarity measures (i.e. correlations).
expType <- "mixed"
nSims <- 100
expDir <- paste0(getwd(),"/Example-", expType,"-z1cor")
expDataFiles <- paste0("corRes-sampBias-nSim",nSims)
nSpecies <- lstObjs$cells$numSpecies
plotColTitles <- ""
plotFile <- paste0("similarityPlot-sampBias-nSim",nSims)
plotSimilarity(expDir, expDataFiles, nSpecies, plotColTitles, plotDevice = "png",
               plotWidth = 10,  plotHeight = 10, plotFileName = plotFile)

