## Map class
setClass(
  Class="Map",
  representation=representation(
    MODE="integer", 
    XYZSTART="integer",  #store number in X, Y, Z of the first cell of the map as a 3-vector
    nXYZ="integer",       #store number of intervals along X, Y and Z as a 3-vector
    mapCRS="integer",     #store axis order as 3-vector (1=X, 2=Y, 3=Z) for writing back to a file	
    Cell="UnitCell",
    symmetry="Symmetry",        #Use the S4 Symmetry class
    LSKFLG="logical",
    SKWMAT="matrix",
    SKWTRN="numeric",
    LABEL="character",
    mapData="array") 
)

##Default methods
# Print map. S4 methods use 'show' rather than 'print'
setMethod(
  f="show",
  signature=signature(object = "Map"),
  definition=function(object)
  {
    cat("*** This is a Map object. ***\n")
    cat(" The cell parameters are:\n")
    cat(sprintf("       a = %7.2f\n       b = %7.2f\n       c = %7.2f\n   alpha =  %6.2f\n    beta =  %6.2f\n   gamma =  %6.2f\n",
      object@Cell@a,object@Cell@b,object@Cell@c,object@Cell@alpha@ang,object@Cell@beta@ang,object@Cell@gamma@ang))
      #note using '@' breaks encapsulation. It would be better to obtain values through accessor methods of the Cell class!
    cat(sprintf(" The space group extended Hermann-Maguin symbol is %s\n",object@symmetry@sym_xHM))
      #again, this breaks encapsulation. Better to use accessor methods of Symmetry.
    cat( " Density values:\n")
    if(typeof(object@mapData) != "complex") cat(sprintf("       min = %7.2f\n       mean = %7.2f\n       max = %7.2f\n",
      min(object@mapData), mean(object@mapData), max(object@mapData)))
  invisible(object)
  }
)

# plot map - currently just a stub to adapt showMap for the new S4 object
#setMethod(
#  f="plot",
#  signature="Map",
#  definition=function(x, ...)
#  {
#    # convert S4 map into old list version
#    slotnames <- slotNames(x) 
#    mapAsList <- vector("list", length(slotnames)) 
#    names(mapAsList) <- slotnames 
#    for(i in slotnames) mapAsList[[i]] <- slot(x, i)
#	
#	showMap(mapAsList, ...)
#  }
#)

##Generic methods

# Check Map object is correct
setGeneric(
  name="checkMap",
  def=function(object,message=FALSE){standardGeneric("checkMap")}
)

setMethod(
  f="checkMap",
  signature="Map",
  definition=function(object,message=FALSE)
  {
    # Check slots are correct objects of corresponding class
    slot_classes <- getSlots(class(object))
    slot_names <- names(slot_classes)
    if (length(slot_classes) != 11) stop("This is not an object of class Map, as it hasn't got 12 slots")
    
    # Check each slot is conformant in turn (see cell_S4.R for good example in checkUnitCell)
    if(slot_names[1] != "MODE") stop("First slot in 'Map' object does not appear to be named 'MODE'")
    if(!is.integer(object@MODE) | length(object@MODE) != 1)
      stop("Second slot in 'Map' object does not appear to be an integer")

    if(slot_names[2] != "XYZSTART") stop("Second slot in 'Map' object does not appear to be named 'XYZSTART'")
    if(!is.integer(object@XYZSTART) | length(object@XYZSTART) != 3)
      stop("Second slot in 'Map' object does not appear to be an integer 3-vector")

    if(slot_names[3] != "nXYZ") stop("Third slot in 'Map' object does not appear to be named 'nXYZ'")
    if(!is.integer(object@nXYZ) | length(object@nXYZ) != 3)
      stop("Third slot in 'Map' object does not appear to be an integer 3-vector")

	if(slot_names[4] != "mapCRS") stop("Fourth slot in 'Map' object does not appear to be named 'mapCRS'")
    if(!is.integer(object@mapCRS) | length(object@mapCRS) != 3)
      stop("Fourth slot in 'Map' object does not appear to be an integer 3-vector")
	  
    if(slot_names[5] != "Cell") stop("Fifth slot in 'Map' object does not appear to be named 'Cell'")
    if(slot_classes[5] != "UnitCell")
      stop("Fifth slot in 'Map' object does not appear to be a 'UnitCell' object")
    checkUnitCell(object@Cell)

    if(slot_names[6] != "symmetry") stop("Sixth slot in 'Map' object does not appear to be named 'symmetry'")
    if(slot_classes[6] != "Symmetry")
      stop("Seventh slot in 'Map' object does not appear to be a 'Symmetry' object")
    checkSymmetry(object@symmetry)

    if(slot_names[7] != "LSKFLG") stop("Seventh slot in 'Map' object does not appear to be named 'LSKFLG'")
    if(!is.logical(object@LSKFLG) | length(object@LSKFLG) != 1)
      stop("Seventh slot in 'Map' object does not appear to be of type 'logical'")

    if(slot_names[8] != "SKWMAT") stop("Eighth slot in 'Map' object does not appear to be named 'SKWMAT'")
    if(!is.numeric(object@SKWMAT) | !identical(dim(object@SKWMAT), as.integer(c(3,3))))
      stop("Eighth slot in 'Map' object does not appear to be a numeric matrix")

    if(slot_names[9] != "SKWTRN") stop("Ninth slot in 'Map' object does not appear to be named 'SKWTRN'")
    if(!is.numeric(object@SKWTRN) | length(object@SKWTRN) != 3)
      stop("Ninth slot in 'Map' object does not appear to be an numeric 3-vector")

    if(slot_names[10] != "LABEL") stop("Tenth slot in 'Map' object does not appear to be named 'LABEL'")
    if(!is.character(object@LABEL) | length(object@LABEL) != 10)
      stop("Tenth slot in 'Map' object does not appear to be a character vector")

    if(slot_names[11] != "mapData") stop("Eleventh slot in 'Map' object does not appear to be named 'mapData'")
    if(slot_classes[11] != "array" | !identical(length(dim(object@mapData)), 3))
      stop("Eleventh slot in 'Map' object does not appear to be a vector array of 3 dimensions")

    # Add other checks

    # All tests passed: green light!
    if (message) print("This object and its slots appear to be correct")          
  }
)

#May change this later to be a method of some other already existing generic function
#This requires rgl and misc3d. How best to do this? Could test for availability and return early if not available?
setGeneric(
  name="showMap",
  def=function(m, ...){standardGeneric("showMap")}
)

setMethod(
  f="showMap",
  signature="Map",
  definition=function(m, sigma = 1.5, colour = "blue", ...)
  {
    if(!require("rgl", quietly=TRUE)) stop("The rgl package must be installed to use this function.\n")
    if(!require("misc3d", quietly=TRUE)) stop("The misc3d package must be installed to use this function.\n")
    
    level <- mean(m@mapData) + sigma * sd(m@mapData)
  
    #get sampling interval size in fractional coordinates
    pitch <- 1 / m@nXYZ
  
    #find map origin (not actually origin, but upper edge in all directions of first cell?)
    o <- m@XYZSTART * pitch 

    # get sampling points in fractional coordinates
    x <- seq(from = o[[1]], by = pitch[[1]], length.out = dim(m@mapData)[1])
    y <- seq(from = o[[2]], by = pitch[[2]], length.out = dim(m@mapData)[2])
    z <- seq(from = o[[3]], by = pitch[[3]], length.out = dim(m@mapData)[3])  
  
    # Use function in the lattice module to get transformation matrix from fractional to
    # cartesian orthonormal coordinates
	cell <- extractCellParameters(m@Cell)
	
    mat <- .triclinic_to_orthogonal_01(cell[[1]], cell[[2]], cell[[3]],
                                    cell[[4]], cell[[5]], cell[[6]])                                    

    # build scene using fractional coordinate system
    scene <- contour3d(m@mapData, level = level, x=x, y=y, z=z, color = colour, draw = FALSE, ...)

    # convert crystallographic fractional coordinates (x,y,z) to Cartesian orthonormal (e1, e2, e3) for display.
    scene$v1 <- scene$v1 %*% mat
    scene$v2 <- scene$v2 %*% mat
    scene$v3 <- scene$v3 %*% mat
  
    # Transform scene by whatever skew transformation is present
    # THIS IS UNTESTED
    if(m@LSKFLG == 1){
      # N.B. the array of vertices is stored in the form
      # v1
      #      [,1]  [,2]  [,3]
      # [1,]   x1    y1    z1
      # [2,]   x2    y2    z2
      # [3,]   x3    y3    z3
      # [4,]   x4    y4    z4
      # ....  ...   ...   ...
      # while the skew matrix is stored in the form
      # SKWMAT
      #      [,1]  [,2]  [,3]
      # [1,]  S11   S21   S31
      # [2,]  S12   S22   S32
      # [3,]  S13   S23   S33
      #
      # therefore to apply the transformation
      # Xo(map) = S * (Xo(atoms) - t)
      # to each vector coordinate (xi, yi, zi)
      # it is easiest to calculate using the transpose of SKWMAT
      translation <- matrix(data = m@SKWTRN, nrow = dim(scene$v1)[1], ncol = 3, byrow=TRUE)
  
      scene$v1 <- scene$v1 - translation
      scene$v1 <- scene$v1 %*% t(m@SKWMAT)

      scene$v2 <- scene$v2 - translation  
      scene$v2 <- scene$v2 %*% t(m@SKWMAT)
  
      scene$v3 <- scene$v3 - translation  
      scene$v3 <- scene$v3 %*% t(m@SKWMAT)
    }
    # may need to experiment with lighting directions
    drawScene.rgl(scene)
    return(invisible(scene))
  }
)

## Friendly constructors
createMapFromCCP4MapFile <- function(filename){
  # read map header, symmetry records and data into a list
  lmap <- .readMap(filename, messages = FALSE)
  h <- lmap$header
  
  # extract cell and create UnitCell object
  cell <- createUnitCell(a=h$X,b=h$Y,c=h$Z,alpha=h$Alpha,beta=h$Beta,gamma=h$Gamma)
  checkUnitCell(cell)
  
  # extract symmetry number and create Symmetry object
  SG <- createSymmetry(h$ISPG)
  checkSymmetry(SG)
  
  # warn if this Symmetry object differs from lmap$symmetryRecords
  #~~~~~~TO DO~~~~~~~
  
  # transpose (if needed) to change column, row, section order to x, y, z order
  mapCRS=c(h$MAPC, h$MAPR, h$MAPS)
  first <- which(mapCRS == 1)
  second <- which(mapCRS == 2)
  third <- which(mapCRS == 3)  
  mapData <- aperm(lmap$mapData, c(first, second, third))
  
  # find the number in X, Y, Z of the first cell in the map
  CRSSTART <- c(h$NCSTART, h$NRSTART, h$NSSTART)
  XYZSTART <- c(CRSSTART[first], CRSSTART[second], CRSSTART[third])
  
  #create object
  object <- new(Class="Map",
    MODE=h$MODE, 
    XYZSTART=XYZSTART,  
    nXYZ=c(h$NX, h$NY, h$NZ),       
    mapCRS=mapCRS, 
	Cell=cell,    
    symmetry=SG,        
    LSKFLG=h$LSKFLG,
    SKWMAT=h$SKWMAT,
    SKWTRN=h$SKWTRN,
    LABEL=h$LABEL,
    mapData=mapData)
  
  return(object)
}
