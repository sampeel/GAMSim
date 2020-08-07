initDomain <- function(){
  
  # Create a domain object.
  domain <- list(owin=NULL,             # an owin object giving the domain (extent with mask and cut to convex hull)
                 ext=NULL,              # an extent object giving the boundary of the domain as originally specified
                 mask=NULL,             # a RasterLayer object indicating which cells are to be included or excluded from the domain
                 maskValue=NULL,        # the value in mask that indicates this cell is to be excluded from the domain
                 proj=NULL,             # the projection (crs as a string) of the domain (of the mask)
                 poly=NULL,             # sp polygon version of the domain.
                 nrows=0,               # number of rows in the raster version of the domain.
                 ncols=0,               # number of columns in the raster version of the domain.
                 isError=FALSE          # error indicator
                 )
  
  # Return value.
  return(domain)
  
}

#-----------------------------------------------------------------------------------------

makeDomain <- function(x, maskValue=NULL, nrows=NULL, ncols=NULL, proj=NULL) {
  
  # Set the domain object.
  # 
  # Arguments ...
  # x:         either a raster extent (in which case provide nrows, ncols and proj too)
  #            or a mask raster (in which case provide maskValue)
  # maskValue:
  # nrows:
  # nrols:
  # proj:      a character or CRS 
  
  # Initialise the domain object.
  domain <- initDomain()
  
  # What class is the argument x?
  if ( inherits(x, "RasterLayer") ) {
    # It is a raster, assume it is a mask of some sort.
    domain$ext <- extent(x)
    domain$nrows <- nrow(x)
    domain$ncols <- ncol(x)
    domain$maskValue <- maskValue
    domain$mask <- x
    domain$proj <- crs(x)
  } else if ( inherits(x, "Extent")) {
    # It is an extent.
    domain$ext <- x
    domain$nrows <- nrows
    domain$ncols <- ncols
    domain$maskValue <- NA
    domain$mask <- raster(x, nrows=nrows, ncols=ncols, crs=proj)
    domain$mask[] <- TRUE    # All cells are in the domain.
    domain$proj <- proj
  } else {
    domain$isError = TRUE
    stop("Class of domain argument x is not valid.")
  }
  
  # Domain as an owin object.
  domain$owin <- as.owin.ext(domain$ext, domain$mask, domain$maskValue)

  
  # Domain as a polygon.
  if ( is.null(mask) ) {
    # This is faster, i think!
    domain$poly <- rasterToPolygons(raster(domain$ext, nrow=1, ncol=1))
  } else {
    domain$poly <- rasterToPolygons(domain$mask, dissolve = TRUE)
  }

  # Return value.
  return(domain)
  
}

#-----------------------------------------------------------------------------------------

is.domain <- function(domain) {
  
  # Test whether the argument is a valid domain object.  Returns true if it is, false otherwise.
  # NB: only tests names of items at this stage, not classes of items!
  
  # Check argument is the right class (as far as we can!)
  if ( !is(domain, "list") ) {
    # The argument is not even a list.  It is not a domain object.
    return(FALSE)
  }
  
  # Get the expected names of the items for a domain object.
  objectItemNames <- names(initDomain())
  
  # Check the domain argument has the same items.
  if ( all(names(domain) %in% objectItemNames) ) {
    # The same item names, hence, a valid domain object.
    return(TRUE)
    
  } else {
    # Not the same item names, hence, an invalid domain object.
    return(FALSE)
  }
  
}

#-----------------------------------------------------------------------------------------

