initSurveys <- function() {
  
  # Initialise or reset the survey object.
  
  # Create the object.
  survey <- list(numSamples = 0,        # The number of samples
                 numSurveys = 0,        # The number of surveys (i.e. clusters of samples)
                 widthSurveys = 0,      # The width of each survey (i.e. cluster width)
                 numGears = 0,          # The number of gear types used to collect data.
                 namesGears = NULL,     # A vector of the names of the gear types (used as column headings and/or identifiers).
                 xy = NULL,             # A two column data.frame that contains the location of the samples (random points within sample's cell)
                 cells = NULL,          # A vector of length numSamples that contains the sample cell number of the cells within the domain (cells of a raster).
                 rowsInCells = NULL,    # The row number of each sample's cell number in the cells object (vector of length numSamples).
                 surveys = NULL,        # The survey number for each sample (a vector with the same row order as cells in this list).
                 areas = NULL,          # The area of each sample (a vector with the same row order as cells in this list).
                 gears = NULL,          # The gear type used to collect each sample (a vector with the same row order as cells in this list).
                 N = NULL,              # The number of each species collected at each sample location (numSamples x numSpecies)
                 Y = NULL,              # The presence (1) or absence (0) of a species at a sample location (numSamples x numSpecies).
                 runNum = 0,            # The run number that has created the values in N.
                 isError = FALSE        # Error indicator.
  )
  
  # Return value.
  return(survey)
  
}
  
#-----------------------------------------------------------------------------------------

resetSimDataSurveys <- function(surveysObj) {
    
  # Reset the simulated data before each new run.
  
  surveysObj <- initSurveys()
  return(surveysObj)
  
}

#-----------------------------------------------------------------------------------------

makeSampleLocations <- function(domainMask, cellObj, numSamples, minArea, maxArea=minArea, 
                                numSurveys=0, myFormula="", coeffs=NULL, covars=NULL, 
                                widthSurveys=3, doPlotCheck = FALSE){
  
  # Simulates the survey sample locations using the given arguments.  The sample location is 
  # the cell of the domain in which the sample has happened.  Each sample has an area 
  # associated with it (random amount between the given limits) but no shape is assigned 
  # to this area.  Each sample area is assumed to occur wholely within its designated cell.  
  # Returns a survey object.
  #
  # Arguments ...
  # domainMask:   a raster layer that defines the projection and resolution of the domain 
  #               (as well as the cells that are included).  
  # cellObj:      a cell object that contains info about each cell in the domain.
  # numSamples:   number of sample locations to simulate
  # minArea:      minimum survey area allowed for simulated survey area
  # maxArea:      maximum survey area allowed for simulated survey area
  # numSurveys:   number of clusters the sample locations are grouped into
  #               numClusters=0 gives random locations (regardless of widthCluster value???).
  # myFormula:    formula for probability of observation to determine cluster locations
  # coeffs:       value of coefficients for this probability of observation.
  # covars:       data.frame of covariates for probability of observations (numCells x numCovars)
  # widthSurveys: width, in number of cells, of the rectangular area of each cluster.
  #               eg. widthClust=3, each cluster will include the nine cells around and 
  #               including its cluster focal cell (cell where the centre of the cluster is).
  # doPlotCheck:  whether or not to check with a simple plot.

  # Initialise the return value.
  survey <- initSurveys()

  # Check that numSamples has been set and is a positive number.
  if ( numSamples < 1 || is.null(numSamples) || is.na(numSamples) ) {
    survey$isError <- TRUE
    stop("The number of samples is not a positive integer.") 
  } else {
    survey$numSamples <- numSamples
  }
  
  # Check that the sample areas are non-zero and valid.
  if ( minArea <= 0.0 || maxArea <= 0.0 ) {
    survey$isError <- TRUE
    stop("Min and/or max area for surveys is not a positive number.")
  }
  
  # Check that maximum sample area is not less than the minimum sample area.
  if ( minArea > maxArea ) {
    tmp <- minArea
    minArea <- maxArea
    maxArea <- tmp
    warning("Minimum and maximum sample areas have been reversed.")
  }
  
  # Select the cells that are included in the clusters.
  if ( numSurveys > 0 ) {
    survey$numSurveys <- numSurveys
    # Check that the widthSurveys value is odd.
    if ( widthSurveys %% 2 == 0 ) {
      warning("Width of a surveys needs to be odd, 1 added to given value.")
      widthSurveys <- widthSurveys + 1
    }
    survey$widthSurveys <- widthSurveys
    
    # Get probability of observation and scale it so that it sums to 1. 
    # NB: it is irrelevant what the intercept is as it cancels out with scaling!
    prObsCell <- lambda.cell(myFormula, coeffs, covars)
    sumPrObs <- sum(prObsCell)
    prObsCell <- prObsCell / sumPrObs
    
    # Get location of surveys (parents of clusters)
    parentCells <- sample(cellObj$cells, numSurveys, replace=FALSE, prObsCell)
    
    # Get potential location of samples (children of clusters)
    children <- myAdjacent(domainMask, parentCells, widthSurveys, idCells = TRUE, include = TRUE)
    names(children) <- c("parentCells", "myCells") 
    numChildren <- dim(children)[1]

    # Figure out which surveys these belong to.
    children$mySurveys <- rep(NA, numChildren) 
    for ( surv in 1:numSurveys ) {
      # find all children belonging to a parent and assign a survey identity to them.
      parentCell <- parentCells[surv]
      indThisParent <- which(children$parentCells == parentCell)
      children$mySurveys[indThisParent] <- surv
    }

  } else if ( numSurveys == 0 ) {
    # Want whole domain and random sampling.  
    children <- data.frame(parentCells=cellObj$cells, 
                           myCells=cellObj$cells, 
                           mySurveys=rep(1, cellObj$numCells))
    numChildren <-cellObj$numCells

  } else {
    # Negative clusters are meaningless!
    survey$isError <- TRUE
    stop("Unable to make a negative number of clusters for surveying locations, check settings.")
  }
    
  # Generate the cells that will contain the required number of samples.
  indChosenChildren <- sample(1:numChildren, numSamples, replace=TRUE)   # uniform sampling of available cells
  survey$cells <- children$myCells[indChosenChildren]
  survey$surveys <- children$mySurveys[indChosenChildren]
  survey$rowsInCells <- match(survey$cells, cellObj$cells)  

  # Generate the random survey points.
  centreCells <- xyFromCell(domainMask, survey$cells)
  cellRes <- res(domainMask)
  x <- centreCells[ ,1] + runif(numSamples, -cellRes[1]/2, cellRes[1]/2)
  y <- centreCells[ ,2] + runif(numSamples, -cellRes[2]/2, cellRes[2]/2)
  survey$xy <- data.frame(x=x, y=y, stringsAsFactors = FALSE)
  
  # Get the area of each sample
  survey$areas <- runif(numSamples, min=minArea, max=maxArea)
  
  # Do plot check?
  if ( doPlotCheck ) {
    plot(domainMask, asp=1, main="Simulated sample locations", col=grey.colors(255))
    if ( numSurveys == 0 ) {
      sampleCol <- "red"
    } else {
      sampleCol <- rainbow(numSurveys)
    }
    points(survey$xy, pch="+", col=sampleCol[survey$surveys])
  }
  
  # Return value.
  return(survey)

}

#-----------------------------------------------------------------------------------------

assignGearUsedSamples <- function(surveyObj, numGears, namesGears,
                                  gearUseStrategy=c("rand","covar","survey"), 
                                  sampleCovar=NULL, rangeDomainCovar=NULL, meanGear=NULL, 
                                  sdGear=NULL, minPrevalence=10, maxit=5) {
  
  # Assign the gear type used to collect each sample.  There is a choice of the strategy 
  # used to assign gears.  Repeats assignment process until all gears are included to 
  # specified level of prevalence or maximum number of repeats are performed (error occurs
  # if maximum number of repeats performed without a successful set of gears being assigned).
  # Otherwise, returns the survey object with the gear info included.
  # 
  # Arguments ...
  # surveyObj:        the survey object that contains the number and location of the samples
  # numGears:         the number of gears available (set in surveyObj here).
  # nameGears:        the names of each gear (set in surveyObj here).
  # gearUseStrategy:  the strategy used to assign gears to samples, one of
  #                     rand   - uniform random assignment of one of the gears to each sample
  #                     covar  - use the given covariate value to give probability of each 
  #                              gear being used for a particular sample (pr(gear) ~ Norm(mean,sd))
  #                     survey - all samples in a survey use the same gear.  Survey gear is
  #                              randomly assigned (uniform dist).
  # sampleCovar:      a vector of covariate values for each sample (only necessary 
  #                   if gearUseStrategy = "covar").  
  # rangeDomainCovar: a vector of length 2 that contains the full range of the covariate in the 
  #                   domain (as those values in sampleCovar might not represent the full range).
  #                   Assumes the minimum is the first value and the maximum is the second value.
  #                   Only necessary if gearUseStrategy = "covar".
  # meanGear:         Sets the mean for each gear's normal distribution ( min(covar) < 
  #                   meanGear[gear] < max(covar) ).  When meanGear is NULL, the range of
  #                   the covariate is divided into numGears equal sections and meanGear
  #                   is the centre of each of these sections. Only necessary if 
  #                   gearUseStrategy = "covar".
  # sdGear:           Sets the standard deviation for each gear's normal distribution.  
  #                   When sdGear is NULL, all gears have an sd = half the width of the
  #                   sections (see meanGear). Only necessary if gearUseStrategy = "covar".
  # minPrevalence:    the minimum number of each gear type that needs to be present for the
  #                   assignment to be accepted (set to zero if don't care)
  # maxit:            the maximum number of iterations to perform to assign a set of gears
  #                   that includes all gears (set to one if want to keep whatever comes out).

  # Initialise values
  numSamples <- surveyObj$numSamples
  gears <- vector("integer",numSamples)
  
  # Check that numGears has been set and is a positive number.
  if ( numGears < 1 || is.null(numGears) || is.na(numGears)) 
    stop("The number of gear types is not a positive integer.") 

  # Check arguments ...
  if ( length(namesGears) != numGears ) 
    stop("Argument 'namesGears' is the wrong length.")
  gearUseStrategy <- match.arg(gearUseStrategy)
  
  # Check other arguments have been given when they are required.
  if ( gearUseStrategy == "covar" ) {
    
    if ( is.null(sampleCovar) ) 
      stop("Argument 'sampleCovar' is missing but is required when gearUseStrategy = covar.")
    if ( length(sampleCovar) != numSamples )
      stop("Argument 'sampleCovar' is the wrong length.")
    if ( is.null(rangeDomainCovar) ) 
      stop("Argument 'rangeDomainCovar' is missing but is required when gearUseStrategy = covar.")
    if ( length(rangeDomainCovar) != 2 )
      stop("Argument 'rangeDomainCovar' is the wrong length (should be = 2).")
    if ( rangeDomainCovar[1] > rangeDomainCovar[2] )
      rangeDomainCovar <- c(rangeDomainCovar[2],rangeDomainCovar[1])
    if ( ! is.null(meanGear) && length(meanGear) != numGears ) 
      stop("Argument 'meanGear' is the wrong length.")
    if ( ! is.null(meanGear) && length(sdGear) != numGears ) 
      stop("Argument 'sdGear' is the wrong length.")
    
  } else if ( gearUseStrategy == "survey") {
    
    if ( is.null(surveyObj$surveys) ) 
      stop("Please set survey for each sample before continuing.")
    if ( surveyObj$numSurveys < numGears ) 
      stop("The number of surveys that have provided samples to the PA data is less than the number of gears.")
  }
  
  if ( (minPrevalence * numGears) > numSamples ) 
    stop("Impossible for this number of each gear to be present for this number of samples.")
  if ( maxit <= 0 ) stop("Argument 'maxit' must be a positive integer.")
  
  # Assign gear type
  gearAssigned <- FALSE
  for ( it in 1:maxit ) {
    if ( gearUseStrategy == "rand" ) {
      # Gear types are assigned uniformly randomly to each observation.
      gears <- sample(1:numGears, size = numSamples, replace = TRUE)
      
    } else if ( gearUseStrategy == "covar" ) {
      # Gear types are assigned based on value of covariate at sample location.
      # Uses multiple overlapping normal distributions to assign probability for each gear
      # being assigned.
      
      # Create normal distributions.
      if ( is.null(meanGear) || is.null(sdGear) ) {
        # Equally spread normal curves.
        binBreaks <- seq(from=floor(rangeDomainCovar[1]), to=ceiling(rangeDomainCovar[2]), 
                         length.out=numGears+1)
        halfBin <- abs(binBreaks[2] - binBreaks[1])/2.0
        means <- binBreaks[1:numGears] + halfBin     # i.e. centre of bins.
        stdevs <- rep(halfBin*2.0, numGears)
      } else {
        # User specified means and sds.
        means <- meanGear
        stdevs <- sdGear
      }
      
      # Get probability of using each gear in each sample (have to select one!)
      gearProbs <- matrix(nrow=numSamples, ncol=numGears)
      for ( g in 1:numGears ) {
        gearProbs[ ,g] <- dnorm(sampleCovar, means[g], stdevs[g])
      }
      sumGearProbs <- matrix(apply(gearProbs, 1, sum), nrow=numSamples, ncol=numGears)
      gearProbs <- gearProbs / sumGearProbs
      
      # Select gear used in each sample.
      for ( j in 1:numSamples ) {
        gears[j] <- sample(1:numGears, size = 1, prob = gearProbs[j, ]) 
      }
  
    } else if ( gearUseStrategy == "survey" ) {
      # The same gear is used for all samples in the same survey.
      surveys <- unique(surveyObj$surveys)
      numSurveys <- length(surveys)
      gearInSurvey <- sample(1:numGears, numSurveys, replace = TRUE)
      for ( s in 1:numSurveys ) {
        # Which samples are in this survey?
        indSamplesInSurvey <- surveyObj$surveys == surveys[s]
        
        # Set these samples' gear.
        gears[indSamplesInSurvey] <- gearInSurvey[s]
      }
    }
    
    # Has prevalence been satisfied for each gear?
    prevGears <- table(gears)
    if ( all(prevGears >= minPrevalence) && (length(prevGears) == numGears) ) {
      surveyObj$numGears <- numGears
      surveyObj$namesGears <- namesGears
      surveyObj$gears <- gears
      gearAssigned <- TRUE
      break
    }
  }
  
  # Return value.
  if ( gearAssigned ) {
    return(surveyObj)
    
  } else {
    surveyObj$isError <- TRUE
    stop("Minimum requirement of ", minPrevalence," of each gear not able to be met in the given number of iterations.")
  }
  
}

#-----------------------------------------------------------------------------------------

is.surveys <- function(surveys) {
  
  # Test whether the argument is a valid surveys object.  Returns true if it is, false otherwise.
  # NB: only tests names of items at this stage, not classes of items!
  
  # Check argument is the right class (as far as we can!)
  if ( !is(surveys, "list") ) {
    # The argument is not even a list.  It is not a surveys object.
    return(FALSE)
  }
  
  # Get the expected names of the items for a surveys object.
  objectItemNames <- names(initsurveys())
  
  # Check the surveys argument has the same items.
  if ( all(names(surveys) %in% objectItemNames) ) {
    # The same item names, hence, a valid surveys object.
    return(TRUE)
    
  } else {
    # Not the same item names, hence, an invalid surveys object.
    return(FALSE)
  }
  
}

#-----------------------------------------------------------------------------------------

testPoissonCluster <- function(nPA, nClusts, childCellWidth=3,
                             nrows=64, ncols=96, simExt=extent(0,12,0,8), doPlots=TRUE) {
  
  # Test my own version of this spatstat function.  Need my own version so that I can
  # control the number of survey points.  Means I'll probably need a few more sim settings!
  #
  # nPA:            number of survey locations required
  # nClusts:        number of clusters of survey locations required (i.e. number of parent
  #                 points)
  # lambda.formula: string that specifies covariates and their relationship in intensity.
  # lambda.coeffs:  coefficients for above formula (including intercept).
  # childCellWidth: number of cells widths (around and including parent cell, i.e. odd nums) 
  #                 Really probably only works if cells are roughly square.  If, for example, 
  #                 they are long, thin rectangles then child points are not going to fall 
  #                 randomly around parent points but have a linear nature to their 
  #                 positioning.  But will do for test!
  # nrows:          number of rows in the raster that forms the domain (and cell structure).
  # ncols:          number of columns in the raster that forms the domain.
  # simExt:         extent of the domain.
  # doPlots:        whether or not to plot results.
  
  # # Numbers of things.
  # ncells <- nrows * ncols
  # 
  # # Set up domain.
  # rlBlank <- raster(simExt, nrows=nrows, ncols=ncols)
  # #rlMask <- setValues(rlBlank, TRUE)
  # #domainObj <- makeDomain(rlMask, NA)
  # 
  # # Set up the covariates
  # centreCells <- as.data.frame(xyFromCell(rlBlank, 1:ncells), stringsAsFactors=FALSE)
  # x <- centreCells$x
  # y <- centreCells$y
  # 
  # # Set up true intensity.
  # trueIntensity <- testLambda(x, y, lambda.formula, lambda.coeffs)
  # if (doPlots) {
  #   rlTrueIntensity <- setValues(rlBlank, trueIntensity)
  #   plot(rlTrueIntensity, asp=1)
  # }
  # 
  # # Make a probability of observation for each cell. (this could be dodgy! should be 1 - exp(-lambda)?)
  # sumTrueIntensity <- sum(trueIntensity)
  # probObsCell <- trueIntensity/sumTrueIntensity
  
  # Select the cells that are to be parent cells.
  parentCells <- sample(1:ncells, nClusts, replace=FALSE, probObsCell)
  if ( doPlots ) {
    # Add cells that contain parent points.
    points(x[parentCells], y[parentCells], pch=".")    
  }
  
  # Select the cells that can contain the children around and including the parent cells.
  if ( childCellWidth %% 2 == 0 ) {
    warning("childCellWidth needs to be odd, 1 added to given value.")
    childCellWidth <- childCellWidth + 1
  }
  neighbourMat <- matrix(1, nrow = childCellWidth, ncol = childCellWidth)
  neighbourMat[ceiling(childCellWidth/2), ceiling(childCellWidth/2)] <- 0    # parent cell!
  childCellPairs <- adjacent(rlBlank, parentCells, neighbourMat, pairs=TRUE, include=TRUE)
  childCellPairs <- as.data.frame(childCellPairs)
  childCells <- unique(childCellPairs$to)
  if ( doPlots ) {
    # Add cells that can contain child points.
    points(x[childCells], y[childCells], pch=".", col="red")
  }
  
  # Generate the cells that will contain the surveys.
  surveyCells <- sample(childCells, nPA, replace=TRUE)
  if ( doPlots ) {
    # Add plot of survey locations (assign x,y coordinates within each cell).
    cellRes <- res(rlBlank)
    xPA <- x[surveyCells] + runif(nPA, -cellRes[1]/2, cellRes[1]/2)
    yPA <- y[surveyCells] + runif(nPA, -cellRes[2]/2, cellRes[2]/2)
    points(xPA, yPA, pch="+")
  }
  
  # How many surveys in each cluster?
  clusters <- data.frame(cluster=1:nClusts, parentCell=parentCells, numChildCells=rep(0.0,nClusts),
                         numSamples=rep(0.0,nClusts), stringsAsFactors = FALSE)
  for ( cell in surveyCells ) {
    # Which is this cell's parent?
    indCell <- which(childCellPairs$to == cell)
    parentsOfCell <- childCellPairs$from[indCell]
    
    # Add proportion of child to each parent.
    numParents <- length(parentsOfCell)
    indParents <- which(parentCells %in% parentsOfCell)
    clusters$numSamples[indParents] <- clusters$numSamples[indParents] + (1.0/numParents)
  }
  
  # How many child cells in each cluster?
  for ( cell in childCells ) {
    # Which is this cell's parent?
    indCell <- which(childCellPairs$to == cell)
    parentsOfCell <- childCellPairs$from[indCell]
    
    # Add proportion of child to each parent.
    numParents <- length(parentsOfCell)
    indParents <- which(parentCells %in% parentsOfCell)
    clusters$numChildCells[indParents] <- clusters$numChildCells[indParents] + (1.0/numParents)
  }
  
  # Return value.
  return(list(cells=surveyCells, x=xPA, y=yPA, clusters=clusters))
  
}

#-----------------------------------------------------------------------------------------

myAdjacent <- function(domain, cells, width, idCells=FALSE, include=TRUE) {
  
  # My version of raster::adjacent function (which is too slow in my sim to use). 
  # Cells must contain cell numbers that are within the domain.
  #
  # Arguments
  # domain:  a raster layer that givens the cells that are included in the domain (i.e. a mask).
  # cells:   a vector containing the cells whose adjacent cells are to be found (but only 
  #          those in the domain)
  # width:   the width (in number of cells) of the box of cells around the given cells that
  #          are considered to be adjacent.  For example, if width = 3, then the eight cells
  #          (to the left, right, top, bottom, and diagonally NW, NE, SE and SW of a given 
  #          cell from cells) are adjacent to "cell" and potentially part of the answer.
  # idCells: return the cell number (from cells) that each returned cell is adjacent to. 
  #          If this is TRUE, a matrix with "cell" and "adjacent" columns is returned.
  # include: include the cell number (from cells) in the result.
  
  # Numbers of things.
  numCols <- domain@ncols
  numCells <- domain@nrows * numCols
  
  # Check cell numbers within argument cells are valid domain cell numbers.
  #domainCells <- 1:numCells 
  domainCells <- cellFromMask(domain, domain)
  if ( length(setdiff(cells, domainCells)) ) {
    stop("One or more values in argument 'cells' are not valid domain cell numbers.")
  }
  
  # Width needs to be odd.
  if ( width %% 2 == 0 ) {
    warning("width needs to be odd, 1 added to given value.")
    width <- width + 1
  }
  
  # What are the potential cells?
  potentialCells <- NULL     #data.frame(cell=NULL, adjacent=NULL, stringsAsFactors = FALSE)
  neighbours <- (1:width) - ((width %/% 2) + 1)
  for ( i in neighbours ) {
    for ( j in neighbours ) {
      # Work out what the cell number of this neighbour is for each cell in cells.
      neighbourCells <- cells + j + (numCols * i)
      potentialCells <- rbind(potentialCells, cbind(cells, neighbourCells))
    }
  }
  
  # Include the original cells in the result?
  if ( ! include ) {
    # No, remove them.
    indCells <- match(potentialCells[ ,2], cells)
    potentialCells <- potentialCells[!indCells, ]
  }
  
  #Get rid of cells that are not in the domain (including cells outside extent!)
  indValidCells <- potentialCells[ ,2] %in% domainCells
  adjacentCells <- as.data.frame(potentialCells[indValidCells, ], stringsAsFactors = FALSE)
  names(adjacentCells) <- c("cell", "adjacent")
  
  # Return value.
  if ( idCells ) {
    # Return the a matrix that also contains a column of cells that the adjacent cells are adjacent to.
    return(adjacentCells)
  } else {
    # Return a vector of the adjacent cells (this may include the original cells as adjacent cells!)
    return(adjacentCells$adjacent)
  }
  
}

#-----------------------------------------------------------------------------------------

makeNSampled <- function(survey, N, prDet=NULL, cellArea=1) {
  
  # Get the number of individuals counted in each sample.  Need to make sure that N is not
  # exceeded in any cell by considering all samples in a cell at the same time.  Returns 
  # the survey object with the number of individuals counted in each sample for each 
  # species, survey$N.  Also returns the presence-absence version of this, survey$Y.
  #
  # Arguments ...
  # survey:     a survey object that contains the cell, area, gear and rowInCells values for each survey
  # N:          the number of individuals in each cell for each species (numCells x numSpecies).
  #             (NB: same order of cells as given in cell object used to determine sample locations!)
  # prDet:      the probability of detecting an individual of a species with a gear type 
  #             (numGears x numSpecies).  If NULL, prDet = 1 for all gear x species.
  # cellArea:   a scalar containing the area of all cells in the domain (assumes all cells 
  #             have the same area).
  
  # Get the species.
  namesSpecies <- names(N)
  numSpecies <- dim(N)[2]
  if ( is.null(namesSpecies) ) {
    warning("Species names are not given by input. Names are being assigned.")
    namesSpecies <- paste0("sp",1:numSpecies)
  }
  
  # Check the prDet matrix dimensions.
  if ( is.null(prDet) ) {
    # Old model of 100% detection of individuals in sample area.
    prDet <- matrix(1.0, nrow=survey$numGears, ncol=numSpecies)
  } else {
    # Check dimensions are consistent with other input.
    numGears <- dim(prDet)[1]
    if ( numGears != survey$numGears ) 
      stop("Number of rows (gears) in 'prDet' matrix is not correct.")
    if ( numSpecies != dim(prDet)[2] ) 
      stop("Number of columns (species) in 'prDet' matrix is not correct.")
  }
  
  # Get the unique sample cells (those cells that contain one or more samples)
  uniqueSampleCells <- table(survey$cells)
  numUniqueSampleCells <- length(uniqueSampleCells)
  
  # Initialise return value.
  NSampled <- as.data.frame(matrix(nrow = survey$numSamples, ncol = numSpecies), stringsAsFactors=FALSE)
  names(NSampled) <- namesSpecies
  
  # Fill each element of the list vector.
  for ( i in 1:numUniqueSampleCells ) {
    # Initialise list of samples in each cell.
    sampleLst <- list(cell=as.integer(names(uniqueSampleCells)[i]), 
                      numSamples=as.integer(uniqueSampleCells[i]), 
                      rowInCells=NULL, 
                      rowsInSurvey=NULL,
                      areas=NULL,
                      gears=NULL)

    
    # What is the row in the survey object? (should be different for each sample in this cell)
    sampleLst$rowsInSurvey <- which(survey$cells %in% sampleLst$cell)
    if ( length(sampleLst$rowsInSurvey) != sampleLst$numSamples ) 
      stop("There should be a rowsInSurvey value for each of the sample in this cell.")
    
    # What is the row in the cell object? (should be the same for all samples in this cell)
    if ( length(unique(survey$rowsInCells[sampleLst$rowsInSurvey])) != 1 ) 
      stop("Same sample cells should have same row in cells object!")
    sampleLst$rowInCells <- survey$rowsInCells[sampleLst$rowsInSurvey[1]]

    # What is the area of each sample in the cell?
    sampleLst$areas <- survey$area[sampleLst$rowsInSurvey]
    
    # What is the gear type used by each sample in the cell?
    sampleLst$gears <- survey$gears[sampleLst$rowsInSurvey]
    
    # Probability of individual being sampled (being in each sample area and being collected by gear).
    probInSample <- matrix(sampleLst$areas/cellArea, nrow = sampleLst$numSamples, ncol = numSpecies,
                           dimnames = list(NULL,namesSpecies))   
    probSampleSpecies <- prDet[sampleLst$gears, ]
    prob <- probInSample * probSampleSpecies

    # How many of each species was collected for the samples in this cell?
    for ( species in namesSpecies ) {
      # Split number of individuals in this cell into those in each sample and those not sampled.
      probThisSpecies <- c(prob[ ,species], 1 - sum(prob[ ,species]))
      tmp <- rmultinom(1, N[sampleLst$rowInCells,species], probThisSpecies)
    
      # Save number sampled.
      NSampled[sampleLst$rowsInSurvey,species] <- tmp[1:sampleLst$numSamples,1]
    }
  }
  
  # FYI: not enough information provided here to do a plot check within this function.
  #      Plot check will need to be located externally.
  
  # Return value.
  survey$N <- NSampled
  survey$Y <- NSampled
  indPos <- which(survey$Y > 0, arr.ind=TRUE)
  survey$Y[indPos] <- 1                        # Convert to indicate presence.
  return(survey)
  
}

#-----------------------------------------------------------------------------------------

plotSampleVals <- function(survey, vals="N", cols=NULL, titles=NULL, lambda=NULL, 
                           sepGear=TRUE, ...) {
  
  # Plot the values in all samples within survey.  Uses xy values within survey and
  # plots the requested values as points on top of a plot of the given intensity. 
  # Plots one column of vals per page.
  # 
  # Arguments ...
  # survey: a survey object defined as a list of items.
  # vals:   name of the item within the survey object to plot against xy values OR
  #         data.frame of z values to plot (with numSample rows in same order as xy)
  # cols:   vector of name (string) or index (integer) of particular columns of data to plot 
  #         (value of NULL will plot all).
  # titles: vector of strings (same as number of columns to be used).  Uses column names  
  #         as titles if NULL.
  # lambda: a raster layer to use as a background to the points, if given.
  # sepGear: put each gear on a separate plot on the same page.
  
  # Get requested item's values from the list object.
  if ( is.character(vals) ) {
    allValues <- survey[[vals]]
  } else {
    allValues <- vals
  }
  
  # Restrict to just the required columns.
  if ( NCOL(allValues) == 1 ) {
    nPlots <- 1
  } else {
    if ( is.null(cols) ) cols <- 1:dim(allValues)[2]
    nPlots <- length(cols)
  }
  
  # Check there are xy values.
  if ( is.null(survey$xy) ) stop("There are no xy values in the given cell object.")
  
  # Titles for plots?
  if ( is.null(titles) ) {
    if ( NCOL(allValues) == 1 ) {
      if ( is.character(vals) ) {
        titles <- vals
      } else {
        titles <- "Sample values"
      }
    } else {
      titles <- names(allValues[ ,cols]) 
    }
  }
  if ( length(titles) != nPlots ) {
    if ( length(titles) == 1 ) {
      titles <- rep(titles, nPlots)
    } else {
      stop("The titles vector is the wrong length for the number of plots requested.")
    }
  }
  
  # Plot.
  if ( sepGear ) {
    opar <- par(mfrow=c(2,2), mar=c(2,4,1,2), oma=c(0,0,1,0))
    
    # Separate plot for each column requested.
    for ( i in 1:nPlots ) {
      thisPlot <- cols[i]
      
      # Separate sub-plot for each gear.
      for ( gt in 1:nGears ) {
        # Plot background intensity layer?
        if ( is.null(lambda) ) {
          plot(range(survey$xy[ ,1]), range(survey$xy[ ,2]), type="n", 
               ylab=paste0("gear = ", namesGears[gt]) )
        } else {
          plot(lambda, asp=1, ylab=paste0("gear = ", namesGears[gt]))
        }
        indThisGear <- which(surveyObj$gears == gt)
        points(survey$xy[indThisGear, ], pch=as.character(allValues[indThisGear,thisPlot]))
      }
      title(titles[i], outer=TRUE) 
      plot.new()
    }
    par(opar)
    
  } else {
    # Separate plot for each column requested.
    for ( i in 1:nPlots ) {
      thisPlot <- cols[i]
      
      # Plot background intensity layer?
      if ( is.null(lambda) ) {
        plot(range(survey$xy[ ,1]), range(survey$xy[ ,2]), type="n", xlab="", ylab="", ...)
      } else {
        plot(lambda, asp=1, ...)
      }
      points(survey$xy, pch=as.character(allValues[ ,thisPlot]), ...)
      title(titles[i]) 
    }
  }
    
}
