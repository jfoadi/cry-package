.readMap <- function(filename, messages = TRUE){

  ################################################################################
  ## a function for reading a CCP4 format Map                                   ##
  ##                                                                            ##
  ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
  ################################################################################
  
  #read header
  header <- .readMapHeader(filename)
  
  # define machine stamps for little and big endian
  littleStamp <- as.raw(c(0x44, 0x41, 0x00, 0x00))
  bigStamp <- as.raw(c(0x11, 0x11, 0x00, 0x00))
  
  endian <- header$endian
  
  #check this is apparently a map file
  if(header$MAP != "MAP ") warning("This may not be a valid map file")

  # Create a connection to binary file and seek to the right point
  f <- file(filename, open="rb")
  seek(f, where = 1024)  
  
  # Read Symmetry records text
  symmRec <- readChar(f, header$NSYMBT)
  numLines <- header$NSYMBT / 80
  first <- 1 + 80*(0:(numLines - 1) )
  last <- (1:numLines)*80
  symmRec <- substring(symmRec, first, last)
  
  # Read map data - ONLY MODE == 2 HAS BEEN TESTED!
  numItems <- header$NC * header$NR * header$NS
  if(header$MODE == 0) map <- readBin(f, "integer", n=numItems, size=1, signed=TRUE, endian = endian)
  else if(header$MODE == 1) map <- readBin(f, "integer", n=numItems, size=2, signed=TRUE, endian = endian)
  else if(header$MODE == 2) map <- readBin(f, "numeric", n=numItems, size=4, endian = endian)
  else if(header$MODE == 3) { #complex number stored as integer pairs
    map <- readBin(f, "integer", n=2*numItems, size=2, endian = endian)
    map <- map * rep(c(1, 1i), times=numItems) #convert to complex
    odds <- subset(map, subset=rep(c(TRUE, FALSE), times=numItems)) #Re part
    evens <- subset(map, subset=rep(c(FALSE, TRUE), times=numItems)) #Im part
    map <- odds + evens #combine
  } else if(header$MODE == 4) map <- readBin(f, "complex", n=numItems, size=8, endian = endian)
  else if(header$MODE == 5) map <- readBin(f, "integer", n=numItems, size=1, signed=TRUE, endian = endian)
  else stop(paste("unrecognised MODE number",header$MODE))

  close(f)
  
  # reshape map data to 3 dimensions: columns, rows and sections.
  dim(map) <- c(header$NC, header$NR, header$NS)
  if(messages){
    cat("the map",filename,"has been read\n")
    cat("summary of header record\n")
    str(header)
    cat("summary of symmetry record\n")
    print(symmRec)
  }
  # Output data are packed in a list  
  data <- list(map,header,symmRec)
  names(data) <- c("mapData", "header", "symmetryRecords")
  return(data)
}

.readMapHeader <- function(filename, endian = "little"){

  ################################################################################
  ## a function for reading the header records of a CCP4 format Map             ##
  ##                                                                            ##
  ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
  ################################################################################

  # define machine stamps for little and big endian
  littleStamp <- as.raw(c(0x44, 0x41, 0x00, 0x00))
  bigStamp <- as.raw(c(0x11, 0x11, 0x00, 0x00))
  expectedStamp <- if(endian == "little") littleStamp else bigStamp
  alternateStamp <- if(endian == "little") bigStamp else littleStamp

  f <- file(filename, open="rb")

  #some function definitions for shorthand
  readInt <- function() readBin(f, "integer", n=1, size=4, endian = endian)
  readReal <- function() readBin(f, "numeric", n=1, size=4, endian = endian)    

  # Read 256 * 4 byte header
  header <- list()
  header$NC <- readInt()
  header$NR <- readInt()
  header$NS <- readInt()
  header$MODE <- readInt() #CCP4 uses mode 2 by default. 0=signed char,1=short,2=float
  header$NCSTART <- readInt()
  header$NRSTART <- readInt()
  header$NSSTART <- readInt()
  header$NX <- readInt()
  header$NY <- readInt()      
  header$NZ <- readInt()      
  header$X <- readReal()
  header$Y <- readReal()
  header$Z <- readReal()            
  header$Alpha <- readReal()
  header$Beta <- readReal()         
  header$Gamma <- readReal()     
  header$MAPC <- readInt()         
  header$MAPR <- readInt()         
  header$MAPS <- readInt()
  # I believe the following 3 items in the header are reals whatever the MODE.
  # See http://conventions.cnb.uam.es/References/maptest.c
  header$AMIN <- readReal()         
  header$AMAX <- readReal()         
  header$AMEAN <- readReal()
  header$ISPG <- readInt()         
  header$NSYMBT <- readInt()
  header$LSKFLG <- as.logical(readInt())
  header$SKWMAT <- matrix(data = readBin(f, "numeric", n=9, size=4, endian = endian), ncol=3, nrow=3)
  header$SKWTRN <- readBin(f, "numeric", n=3, size=4, endian = endian)
  header$FutureUse <- readBin(f, "integer", n=15, size=4, endian = endian)
  header$MAP <- readChar(f, 4)
  header$MACHST <- readBin(f, "raw", 4, endian = endian) #is 44 41 00 00 in hex on little endian systems
  header$ARMS <- readReal()
  header$NLABL <- readInt()
  header$LABEL <- readChar(f, 800)
  header$LABEL <- substring(header$LABEL, first=1+80*(0:9), last=(1:10)*80)
  
  close(f)
  
  # check byte order
  header$endian <- endian
  if(!identical(header$MACHST, expectedStamp)){
    
    # check byte order is at least recognised
    if(!identical(header$MACHST, alternateStamp)) stop(filename, " is not a readable map file")
    
    # try to read with the other byte order
    header <- if(endian == "little") .readMapHeader(filename, endian="big") else .readMapHeader(filename, endian="little")    
  }
  return(header)
}

#showMap <- function(map, sigma = 1.5, colour = "blue", ...){
#
#  ################################################################################
#  ## a function for displaying a CCP4 format Map with an interactive plot using ##
#  ## rgl and misc3d functions (3D data), or a contour plot (2D data)            ##
#  ##                                                                            ##
#  ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
#  ################################################################################
#
#  require(rgl)
#  require(misc3d)
#  
#  level <- map$header$AMEAN + sigma * map$header$ARMS
#  
#  mapData <- map$mapData
#  
#  # transpose to change column, row, section order to x, y, z order
#  # make look-up vector CRS
#  CRS <- c(map$header$MAPC, map$header$MAPR, map$header$MAPS)
#  # find order of permutation to get map data in x,y,z order
#  first <- which(CRS == 1)
#  second <- which(CRS == 2)
#  third <- which(CRS == 3)
#  mapData <- aperm(mapData, c(first, second, third))
#
#  #get sampling interval size in fractional coordinates
#  xPitch <- 1 / map$header$NX
#  yPitch <- 1 / map$header$NY
#  zPitch <- 1 / map$header$NZ
#  
#  #find map origin
#  CRSSTART <- c(map$header$NCSTART, map$header$NRSTART, map$header$NSSTART)
#  XYZSTART <- c(CRSSTART[first] * xPitch, CRSSTART[second] * yPitch, CRSSTART[third] * zPitch)
#
#  # get sampling points in fractional coordinates
#  x <- seq(from = XYZSTART[[1]], by = xPitch, length.out = dim(mapData)[1])
#  y <- seq(from = XYZSTART[[2]], by = yPitch, length.out = dim(mapData)[2])
#  z <- seq(from = XYZSTART[[3]], by = zPitch, length.out = dim(mapData)[3])  
#  
#  # Use function in the lattice module to get transformation matrix from fractional to
#  # cartesian orthonormal coordinates
#  mat <- triclinic_to_orthogonal_01(map$header$X, map$header$Y, map$header$Z,
#                                    map$header$Alpha, map$header$Beta, map$header$Gamma)                                    
#
#  #check for 2D slice in Y-Z plane
#  if(dim(mapData)[1] == 1){
#    slice <- mapData[1,,]
#    
#    #make 2D contour lines in fractional coords
#    # CAREFUL! In this case cl$x corresponds to y and cl$y to z
#    cl <- contourLines(y, z, slice)
#    
#    #transform coordinates for every contour
#    contourTransform <- function(cont){
#      # add x coordinate of this slice
#      coords <- cbind(rep(x, length.out = length(cont$x)), cont$x, cont$y)
#      # do transformation
#      coords <- coords %*% mat
#      cont$x <- coords[,2] # i.e. the y coord
#      cont$y <- coords[,3] # i.e. the z-coord
#      return(cont)
#    }
#    cl <- lapply(cl, contourTransform)
#    
#    # get xlim and ylim for plot (corresponding to y and z axes of the cell respectively)
#    corner1 <- c(x, y[[1]], z[[1]])
#    corner2 <- c(x, y[[length(y)]], z[[1]])
#    corner3 <- c(x, y[[length(y)]], z[[length(z)]])
#    corner4 <- c(x, y[[1]], z[[length(z)]])
#    edges <- rbind(corner1, corner2, corner3, corner4)
#    edges <- edges %*% mat
#    xlim <- range(edges[,2])
#    ylim <- range(edges[,3])
#
#    plot(NULL, xlim = xlim, ylim = ylim, asp=1, xlab="Y", ylab="Z",
#         main = "map slice plotted against orthonormal axes")
#    # draw plot
#    for(i in seq(along.with = cl)) lines(cl[[i]]$x, cl[[i]]$y, col = colour)
#    return(invisible(cl))
#  }
#    
#  #check for 2D slice in X-Z plane
#  if(dim(mapData)[2] == 1){
#    slice <- mapData[,1,]
#    
#    #make 2D contour lines in fractional coords
#    # CAREFUL! In this case cl$y corresponds to the real z axis
#    cl <- contourLines(x, z, slice)
#    
#    #transform coordinates for every contour
#    contourTransform <- function(cont){
#      # add y coordinate of this slice
#      coords <- cbind(cont$x, rep(y, length.out = length(cont$x)), cont$y)
#      # do transformation
#      coords <- coords %*% mat
#      cont$x <- coords[,1]
#      cont$y <- coords[,3] # i.e. the z-coord
#      return(cont)
#    }
#    cl <- lapply(cl, contourTransform)
#    
#    # get xlim and ylim for plot (ylim corresponds to the z axes of the cell)
#    corner1 <- c(x[[1]], y,  z[[1]])
#    corner2 <- c(x[[length(x)]], y, z[[1]])
#    corner3 <- c(x[[length(x)]], y, z[[length(z)]])
#    corner4 <- c(x[[1]], y, z[[length(z)]])
#    edges <- rbind(corner1, corner2, corner3, corner4)
#    edges <- edges %*% mat
#    xlim <- range(edges[,1])
#    ylim <- range(edges[,3])
#
#    plot(NULL, xlim = xlim, ylim = ylim, asp=1, xlab="X", ylab="Z",
#         main = "map slice plotted against orthonormal axes")
#    # draw plot
#    for(i in seq(along.with = cl)) lines(cl[[i]]$x, cl[[i]]$y, col = colour)
#    return(invisible(cl))
#  }
#  
#  #check for 2D slice in X-Y plane
#  if(dim(mapData)[3] == 1){
#    slice <- mapData[,,1]
#    
#    #make 2D contour lines in fractional coords
#    cl <- contourLines(x, y, slice)
#    
#    #transform coordinates for every contour
#    contourTransform <- function(cont){
#      # add z coordinate of this slice
#      coords <- cbind(cont$x, cont$y, rep(z, length.out = length(cont$x)))
#      # do transformation
#      coords <- coords %*% mat
#      cont$x <- coords[,1]
#      cont$y <- coords[,2]
#      return(cont)
#    }
#    cl <- lapply(cl, contourTransform)
#    
#    # get xlim and ylim for plot
#    corner1 <- c(x[[1]], y[[1]], z)
#    corner2 <- c(x[[length(x)]], y[[1]], z)
#    corner3 <- c(x[[length(x)]], y[[length(y)]], z)
#    corner4 <- c(x[[1]], y[[length(y)]], z)
#    edges <- rbind(corner1, corner2, corner3, corner4)
#    edges <- edges %*% mat
#    xlim <- range(edges[,1])
#    ylim <- range(edges[,2])
#
#    plot(NULL, xlim = xlim, ylim = ylim, asp=1, xlab="X", ylab="Y",
#         main = "map slice plotted against orthonormal axes")    
#    # draw plot
#    for(i in seq(along.with = cl)) lines(cl[[i]]$x, cl[[i]]$y, col = colour)
#    return(invisible(cl))
#  }
#  
#  # build scene using fractional coordinate system
#  scene <- contour3d(mapData, level = level, x=x, y=y, z=z, color = colour, draw = FALSE, ...)
#
#  # convert crystallographic fractional coordinates (x,y,z) to Cartesian orthonormal (e1, e2, e3) for display.
#  scene$v1 <- scene$v1 %*% mat
#  scene$v2 <- scene$v2 %*% mat
#  scene$v3 <- scene$v3 %*% mat
#  
#  # Transform scene by whatever skew transformation is present
#  # THIS IS UNTESTED
#  if(map$header$LSKFLG == 1){
#  # N.B. the array of vertices is stored in the form
#  # v1
#  #      [,1]  [,2]  [,3]
#  # [1,]   x1    y1    z1
#  # [2,]   x2    y2    z2
#  # [3,]   x3    y3    z3
#  # [4,]   x4    y4    z4
#  # ....  ...   ...   ...
#  # while the skew matrix is stored in the form
#  # SKWMAT
#  #      [,1]  [,2]  [,3]
#  # [1,]  S11   S21   S31
#  # [2,]  S12   S22   S32
#  # [3,]  S13   S23   S33
#  #
#  # therefore to apply the transformation
#  # Xo(map) = S * (Xo(atoms) - t)
#  # to each vector coordinate (xi, yi, zi)
#  # it is easiest to calculate using the transpose of SKWMAT
#  translation <- matrix(data = map$header$SKWTRN, nrow = dim(scene$v1)[1], ncol = 3, byrow=TRUE)
#  
#  scene$v1 <- scene$v1 - translation
#  scene$v1 <- scene$v1 %*% t(map$header$SKWMAT)
#
#  scene$v2 <- scene$v2 - translation  
#  scene$v2 <- scene$v2 %*% t(map$header$SKWMAT)
#  
#  scene$v3 <- scene$v3 - translation  
#  scene$v3 <- scene$v3 %*% t(map$header$SKWMAT)
#  }
#  # may need to experiment with lighting directions, esp. for 2D slices
#  drawScene.rgl(scene)
#  return(invisible(scene))
#}
#
#
#
#writeMap <- function(map, filename){
#
#  ################################################################################
#  ## a function for writing a CCP4 format Map                                   ##
#  ##                                                                            ##
#  ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
#  ################################################################################
#
#  # Create a connection to binary file
#  f <- file(filename, open="wb")
#
#  #some function definitions for shorthand.
#  writeInt <- function(obj){
#    if(is.null(obj)) stop(paste("map object does not contain "), deparse(substitute(obj)))
#    writeBin(as.integer(obj), f, size=4)
#  }
#  writeReal <- function(obj){
#    if(is.null(obj)) stop(paste("map object does not contain "), deparse(substitute(obj)))  
#    writeBin(as.double(obj), f, size=4)
#  }
# 
#  #start writing header
#  header <- map$header
#  writeInt(header$NC)
#  writeInt(header$NR)
#  writeInt(header$NS)
#  writeInt(header$MODE) #CCP4 uses mode 2 by default. 0=signed char,1=short,2=float
#  writeInt(header$NCSTART)
#  writeInt(header$NRSTART)
#  writeInt(header$NSSTART)
#  writeInt(header$NX)
#  writeInt(header$NY)
#  writeInt(header$NZ)
#  writeReal(header$X)
#  writeReal(header$Y)
#  writeReal(header$Z)    
#  writeReal(header$Alpha)
#  writeReal(header$Beta)    
#  writeReal(header$Gamma) 
#  writeInt(header$MAPC)   
#  writeInt(header$MAPR)  
#  writeInt(header$MAPS)
#  writeReal(header$AMIN)     
#  writeReal(header$AMAX)     
#  writeReal(header$AMEAN)
#  writeInt(header$ISPG)   
#  writeInt(header$NSYMBT)
#  writeInt(header$LSKFLG)
#  writeReal(header$SKWMAT)
#  writeReal(header$SKWTRN)
#  writeInt(header$FutureUse)
#  writeChar(header$MAP, f, eos=NULL)
#  #machine stamp
#  stamp <- if(.Platform$endian == "little") as.raw(c(0x44, 0x41, 0x00, 0x00)) else as.raw(c(0x11, 0x11, 0x00, 0x00))
#  writeBin(stamp, f)
#  writeReal(header$ARMS)
#  writeInt(header$NLABL)
#  writeChar(header$LABEL, f, eos=NULL)
#  
#  #Now write symmetry records
#  writeChar(map$symmetryRecords, f, eos=NULL)
#  
#  #Now write map data
#  if(header$MODE == 0) writeBin(as.integer(map$mapData), f, size=1)
#  else if(header$MODE == 1) writeBin(as.integer(map$mapData), f, size=2)
#  else if(header$MODE == 2) writeReal(map$mapData)
#    #complex number stored as integer pairs
#  else if(header$MODE == 3){
#    for (i in 1:length(map$mapData)){
#      re <- Re(map$mapData[[i]])
#      im <- Im(map$mapData[[i]])
#      writeBin(as.integer(re), f, size=2)
#      writeBin(as.integer(im), f, size=2)
#    }
#  } else if(header$MODE == 4) writeBin(as.complex(map$mapData), f, size=8)
#  else if(header$MODE == 5) map <- writeBin(as.integer(map$mapData), f, size=1)
#  else stop(paste("unrecognised MODE number",header$MODE))
#
#  close(f)
#}
#