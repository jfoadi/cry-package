# R code to implement the crystallographic ideas related to crystal lattices.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Default methods
#
# Print
setMethod(
          f="print",
          signature="BravaisType",
          definition=function(x,...)
                     {
                      # Extract UnitCell object
                      bl <- x@bl

                      # Print
                      cat(paste("*** This is a BravaisType object. The associated type is ",bl,". ***\n",sep=""))
                     }
         )
#
# Print
setMethod(
          f="print",
          signature="Lattice",
          definition=function(x,...)
                     {
                      # Extract UnitCell object
                      unit_cell <- x@cell

                      # Print
                      cat("*** This is a Lattice object. ***\n")
                      cat("    The unit cell has the following parameters:\n")
                      stringa <- sprintf("       a = %7.2f\n       b = %7.2f\n       c = %7.2f\n   alpha =  %6.2f\n    beta =  %6.2f\n   gamma =  %6.2f\n",
                                         unit_cell@a,unit_cell@b,unit_cell@c,unit_cell@alpha@ang,unit_cell@beta@ang,unit_cell@gamma@ang)
                      cat(stringa)

                      # Extract lattice information
                      latt_info <- extractLatticeStuff(x)

                      # Print
                      stringa <- sprintf("\n    The Bravais lattice is %s",latt_info[1])
                      cat(stringa)
                      stringa <- sprintf("\n    The crystal family is %s",latt_info[2])
                      cat(stringa)
                      stringa <- sprintf("\n    The crystal system is %s",latt_info[3])
                      cat(stringa)
                      stringa <- sprintf("\n    The lattice system is %s\n",latt_info[4])
                      cat(stringa)
                     }
         )

## Methods sharing generic functions with other objects
#
# Check Bravais lattice type
setMethod(
          f="isCorrect",
          signature="BravaisType",
          definition=function(object,name,message=FALSE)
                     {
                      # Stop if name is not "BravaisType"
                      if (name != "BravaisType")
                      {
                       msg <- paste("This is not an object of class ",name,". It is an object of class BravaisType.",sep="")
                       stop(msg)
                      }

                      # Check BravaisType character has length 2
                      if (nchar(object@bl) != 2) stop("Object of class 'BravaisType' has number of characters different from 2")

                      # Extract lattice and centering (also check on symbol validity)
                      latt <- substr(object@bl,1,1)
                      centring <- substr(object@bl,2,2)
                      if (latt != "a" & latt !="m" & latt != "o" & latt != "t" & latt != "h" & latt != "c") 
                                                          stop("Invalid lattice type")
                      if (centring != "P" & centring != "A" & centring != "B" & centring != "C" & centring != "F" &
                          centring != "I" & centring != "R") stop("Invalid centring type")

                      # Now check validity for specific families
                      # Triclinic family
                      if (latt == "a" & centring != "P") stop("Invalid Bravais lattice type")

                      # Monoclinic family
                      #if ((latt == "m" & centring == "I") | (latt == "m" & centring == "F") | (latt == "m" & centring == "R")) stop("Invalid Bravais lattice type")
                      if ((latt == "m" & centring == "F") | (latt == "m" & centring == "R")) stop("Invalid Bravais lattice type")

                      # Orthorombic family
                      if ((latt == "o" & centring != "P") & (latt == "o" & centring != "A") & (latt == "o" & centring != "B") & (latt == "o" & centring != "C") &
                          (latt == "o" & centring != "F") & (latt == "o" & centring != "I")) stop("Invalid Bravais lattice type")

                      # Tetragonal family
                      if ((latt == "t" & centring != "P") & (latt == "t" & centring != "I")) stop("Invalid Bravais lattice type")

                      # Hexagonal and Trigonal families
                      if (latt == "h"  & (centring != "P" & centring != "R")) stop("Invalid Bravais lattice type")

                      # Cubic family
                      if (latt == "c" & (centring == "A" | centring == "B" | centring == "C")) stop("Invalid Bravais lattice type")

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
#
# Check compatibility of cell with lattice and existence of Bravais lattice type
setMethod(
          f="isCorrect",
          signature="Lattice",
          definition=function(object,name,message=FALSE)
                     {
                      # Stop if name is not "Lattice"
                      if (name != "Lattice")
                      {
                       msg <- paste("This is not an object of class ",name,". It is an object of class Lattice.",sep="")
                       stop(msg)
                      }

                      # Extract slot names 
                      slot_names <- getSlots(class(object))

                      # Check object has 2 slots
                      if (length(slot_names) != 2) stop("Object 'Lattice' has a number of slots different from 2")
                      
                      # Check nature of first slot (it should be UnitCell)
                      if (names(slot_names[1]) != "cell") stop("First slot does not appear to be named 'cell'")
                      isCorrect(object@cell,"UnitCell") # Stop here if wrong UnitCell object

                      # Check nature of second slot (it should be BravaisType)
                      if (slot_names[2] != "BravaisType") stop("Second slot is not an object of class 'BravaisType'")
                      isCorrect(object@bl,"BravaisType")

                      # Now proceed with further checks
                      # Extract cell parameters
                      a <- object@cell@a
                      b <- object@cell@b
                      c <- object@cell@c
                      alpha <- object@cell@alpha
                      if (alpha@rad_flag) alpha <- radToDeg(alpha)
                      aa <- alpha@ang
                      beta <- object@cell@beta
                      if (beta@rad_flag) beta <- radToDeg(beta)
                      bb <- beta@ang
                      gamma <- object@cell@gamma
                      if (gamma@rad_flag) gamma <- radToDeg(gamma)
                      cc <- gamma@ang

                      # Extract lattice and centring
                      latt <- substr(object@bl@bl,1,1)
                      centring <- substr(object@bl@bl,2,2)

                      # Compatibility between unit cell and lattice type
                      if (latt == "m")                                    # Monoclinic family
                      {
                       erang <- 0.000001      # Finite accuracy of binary representation of numbers means we have to
                       diff_a <- abs(aa-90)   # test number == 90 in this way
                       diff_b <- abs(bb-90)
                       diff_c <- abs(cc-90)
                       if (diff_a > erang & (diff_b > erang | diff_c > erang)) stop ("Unit cell angles non compatible with lattice type")
                       if (diff_b > erang & (diff_a > erang | diff_c > erang)) stop ("Unit cell angles non compatible with lattice type")
                       if (diff_c > erang & (diff_a > erang | diff_b > erang)) stop ("Unit cell angles non compatible with lattice type")
                       if (diff_a > erang & centring == "A") stop ("Unit cell parameters non compatible with lattice type")
                       if (diff_b > erang & centring == "B") stop ("Unit cell parameters non compatible with lattice type")
                       if (diff_c > erang & centring == "C") stop ("Unit cell parameters non compatible with lattice type")
                      }
                      if (latt == "o")                                    # Orthorombic family
                      {
                       if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
                       {
                        stop ("Unit cell parameters non compatible with lattice type")
                       }
                      }
                      if (latt == "t")                                    # Tetragonal family
                      {
                       if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
                       {
                        stop ("Unit cell parameters non compatible with lattice type")
                       }
                       if (abs(a-b) > 0.000001 & abs(a-c) > 0.000001 & abs(b-c) > 0.000001)
                       {
                        stop ("Unit cell parameters non compatible with lattice type") 
                       }
                      }
                      if (latt == "c")                                    # Cubic family
                      {
                       if (abs(a-b) > 0.000001 | abs(a-c) > 0.000001 | abs(b-c) > 0.000001) 
                       {
                        stop ("Unit cell parameters non compatible with lattice type")
                       }
                       if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001) 
                       {
                        stop ("Unit cell parameters non compatible with lattice type")
                       }
                      }
                      if (latt == "h")                                    # Hexagonal family
                      {
                       if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type") 
                       if (abs(a-b) < 0.000001)
                       {
                        if (abs(cc-120) > 0.000001)
                        {
                         if (abs(b-c) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
                        }
                       }
                      }

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
#
# Extract BravaisType slots and return them in a list
setMethod(
          f="getFields",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Extract slots
                      bl <- object@bl
                      lista <- list(bl=bl)

                      return(lista)
                     } 
         )
#
# Extract Lattice slots and return them in a list
setMethod(
          f="getFields",
          signature="Lattice",
          definition=function(object)
                     {
                      # Extract slots
                      cell <- object@cell
                      bl <- object@bl
                      lista <- list(cell=cell,bl=bl)

                      return(lista)
                     } 
         )
#
# Extract parameters from BravaisType object
setMethod(
          f="getParameters",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Build list to be returned
                      lista <- list(bl=object@bl)

                      return(lista)
                     }
         )
#
# Extract parameters from Lattice object
setMethod(
          f="getParameters",
          signature="Lattice",
          definition=function(object)
                     {
                      # Call auxiliary function
                      tmp <- extractLatticeStuff(object)
                      latt_info <- c(crystal_family=tmp[2],crystal_system=tmp[3],lattice_system=tmp[4])

                      # Extract UnitCell object parameters
                      tmp2 <- object@cell

                      # Build list to be returned
                      lista <- list(a=tmp2@a,b=tmp2@b,c=tmp2@c,alpha=tmp2@alpha@ang,beta=tmp2@beta@ang,gamma=tmp2@gamma@ang,
                                    lattice_info=latt_info,bl=object@bl@bl)

                      return(lista)
                     }
         )
#
setMethod(
          f="getData",
          signature="BravaisType",
          definition=function(object)
                     {
                      # No data are contained in an object of class BravaisType, so this function returns NULL

                      return(NULL)
                     }
         )
#
setMethod(
          f="getData",
          signature="Lattice",
          definition=function(object)
                     {
                      # No data are contained in an object of class Lattice, so this function returns NULL

                      return(NULL)
                     }
         )
#
# Another way to extract BravaisType slots (similar to getFields)
setMethod(
          f="[",
          signature="BravaisType",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "bl") flist <- c(flist,list(bl=x@bl))
                       if (i[ii] == 1) flist <- c(flist,list(bl=x@bl))
                       if (i[ii] != "bl" & i[ii] != 1) flist <- c(flist,list(NA))
                      }

                      return(flist)
                     }
         )
#
# Change parameters values in object of class BravaisType
setReplaceMethod(
          f="[",
          signature="BravaisType",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: bl (character)
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "bl") x@bl <- value[ii]
                       if (i[ii] == 1) x@bl <- value[ii]
                       if (i[ii] != "bl" & i[ii] != 1) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      isCorrect(x,"BravaisType")

                      return(x)
                     }
         )
#
# Another way to extract Lattice slots (similar to getFields)
setMethod(
          f="[",
          signature="Lattice",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Create vector with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") flist <- c(flist,list(a=x@cell@a))
                       if (i[ii] == 1) flist <- c(flist,list(a=x@cell@a))
                       if (i[ii] == "b") flist <- c(flist,list(b=x@cell@b))
                       if (i[ii] == 2) flist <- c(flist,list(b=x@cell@b))
                       if (i[ii] == "c") flist <- c(flist,list(c=x@cell@c))
                       if (i[ii] == 3) flist <- c(flist,list(c=x@cell@c))
                       if (i[ii] == "alpha") flist <- c(flist,list(alpha=x@cell@alpha@ang))
                       if (i[ii] == 4) flist <- c(flist,list(alpha=x@cell@alpha@ang))
                       if (i[ii] == "beta") flist <- c(flist,list(beta=x@cell@beta@ang))
                       if (i[ii] == 5) flist <- c(flist,list(beta=x@cell@beta@ang))
                       if (i[ii] == "gamma") flist <- c(flist,list(gamma=x@cell@gamma@ang))
                       if (i[ii] == 6) flist <- c(flist,list(gamma=x@cell@gamma@ang))
                       if (i[ii] == "bl") flist <- c(flist,list(bl=x@bl@bl))
                       if (i[ii] == 7) flist <- c(flist,list(bl=x@bl@bl)) 
                       if (i[ii] != "a" & i[ii] != 1 &
                           i[ii] != "b" & i[ii] != 2 &
                           i[ii] != "c" & i[ii] != 3 &
                           i[ii] != "alpha" & i[ii] != 4 &
                           i[ii] != "beta" & i[ii] != 5 &
                           i[ii] != "gamma" & i[ii] != 6 &
                           i[ii] != "bl" & i[ii] != 7
                          ) flist <- c(flist,list(NA))
                      }

                      return(flist)
                     }
         )
#
# Change parameters values in object of class Lattice
setReplaceMethod(
          f="[",
          signature="Lattice",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: bl (character)
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") x@cell@a <- value[ii]
                       if (i[ii] == 1) x@cell@a <- value[ii]
                       if (i[ii] == "b") x@cell@b <- value[ii]
                       if (i[ii] == 2) x@cell@b <- value[ii]
                       if (i[ii] == "c") x@cell@c <- value[ii]
                       if (i[ii] == 3) x@cell@c <- value[ii]
                       if (i[ii] == "alpha") x@cell@alpha@ang <- value[ii]
                       if (i[ii] == 4) x@cell@alpha@ang <- value[ii]
                       if (i[ii] == "beta") x@cell@beta@ang <- value[ii]
                       if (i[ii] == 5) x@cell@beta@ang <- value[ii]
                       if (i[ii] == "gamma") x@cell@gamma@ang <- value[ii]
                       if (i[ii] == 6) x@cell@gamma@ang <- value[ii]
                       if (i[ii] == "bl") x@bl@bl <- value[ii]
                       if (i[ii] == 7) x@bl@bl <- value[ii]
                       if (i[ii] != "a" & i[ii] != 1 &
                           i[ii] != "b" & i[ii] != 2 &
                           i[ii] != "c" & i[ii] != 3 &
                           i[ii] != "alpha" & i[ii] != 4 &
                           i[ii] != "beta" & i[ii] != 5 &
                           i[ii] != "gamma" & i[ii] != 6 &
                           i[ii] != "bl" & i[ii] != 7
                          ) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      isCorrect(x,"Lattice")

                      return(x)
                     }
         )

## Methods specific for BravaisType and Lattice classes
#
# Returns centering operators
setMethod(
          f="getCentringOps",
          signature="Lattice",
          definition=function(object){
                                      # Centring operators are a list of 3D vectors

                                      # P
                                      if (substr(object@bl@bl,2,2) == "P") lop <- list(c(0,0,0))

                                      # A
                                      if (substr(object@bl@bl,2,2) == "A") lop <- list(c(0,0,0),c(0,0.5,0.5))

                                      # B
                                      if (substr(object@bl@bl,2,2) == "B") lop <- list(c(0,0,0),c(0.5,0,0.5))

                                      # C
                                      if (substr(object@bl@bl,2,2) == "C") lop <- list(c(0,0,0),c(0.5,0.5,0))

                                      # I
                                      if (substr(object@bl@bl,2,2) == "I") lop <- list(c(0,0,0),c(0.5,0.5,0.5))

                                      # F
                                      if (substr(object@bl@bl,2,2) == "F") lop <- list(c(0,0,0),c(0,0.5,0.5),c(0.5,0,0.5),c(0.5,0.5,0))

                                      # R
                                      if (substr(object@bl@bl,2,2) == "R") lop <- list(c(0,0,0),c(2/3,1/3,1/3),c(1/3,2/3,2/3))

                                      return(lop)
                                     }
         )

## Friendly constructors
#
# Create a BravaisType onject starting from a string
createBravaisType <- function(bl="aP")
{
 # Build object without thinking if it's correct
 object <- new(Class="BravaisType",bl=bl)

 # Now check if object just built is correct
 isCorrect(object,"BravaisType")

 # If things are OK, then return object
 return(object)
}
#
# Assign values to a Lattice object while creating it
createLattice <- function(a=1,b=1,c=1,alpha=90,beta=90,gamma=90,bl="cP")
{
 # Create UnitCell object
 unit_cell <- createUnitCell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

 # Create BravaisType object
 bl <- createBravaisType(bl=bl)

 # If everything ok merge two objects into one
 object <- new(Class="Lattice",cell=unit_cell,bl=bl)

 # Check compatibility between cell and lattice, and other small things
 isCorrect(object,"Lattice")

 return(object)
}

## Other or auxiliary functions, just used internally
#
# Extract crystal family, crystal system and lattice system
setGeneric(
           name="extractLatticeStuff",
           def=function(object){standardGeneric("extractLatticeStuff")}
          )
#
setMethod(
          f="extractLatticeStuff",
          signature="Lattice",
          definition=function(object)
                     {
                      # Extract UnitCell and BravaisType objects
                      unit_cell <- object@cell
                      bravais_type <- object@bl

                      # Extract cell parameters
                      a <- unit_cell@a
                      b <- unit_cell@b
                      c <- unit_cell@c
                      alpha <- unit_cell@alpha@ang
                      beta <- unit_cell@beta@ang
                      gamma <- unit_cell@gamma@ang

                      # Extract crystal family and centring
                      latt <- substr(bravais_type@bl,1,1)
                      centring <- substr(bravais_type@bl,2,2)

                      # Determine crystal family, crystal system and lattice system 
                      if (latt == "a")
                      {
                       cr_fam <- "triclinic"
                       cr_sys <- "triclinic"
                       lt_sys <- "triclinic"
                      }
                      if (latt == "m")
                      {
                       cr_fam <- "monoclinic"
                       cr_sys <- "monoclinic"
                       lt_sys <- "monoclinic"
                      }
                      if (latt == "o")
                      {
                       cr_fam <- "orthorombic"
                       cr_sys <- "orthorombic"
                       lt_sys <- "orthorombic"
                      }
                      if (latt == "t")
                      {
                       cr_fam <- "tetragonal"
                       cr_sys <- "tetragonal"
                       lt_sys <- "tetragonal"
                      }
                      if (latt == "c")
                      {
                       cr_fam <- "cubic"
                       cr_sys <- "cubic"
                       lt_sys <- "cubic"
                      }
                      if (latt == "h")
                      {
                       cr_fam <- "hexagonal"
                       if (centring == "R")
                       {
                        cr_sys <- "trigonal"
                        if (a == b & b == c & alpha == beta & beta == gamma) lt_sys <- "rhombohedral"
                        if (!(a == b & b == c & alpha == beta & beta == gamma)) lt_sys <- "hexagonal (centred)"
                       }
                       if (centring == "P")
                       {
                        cr_sys <- "hexagonal"
                        lt_sys <- "hexagonal"
                       }
                      }

                      # Store info in a 4D character vector
                      latt_info <- c(bravais_type@bl,cr_fam,cr_sys,lt_sys)

                      return(latt_info)
                     }
          )
#
# Correspondence between bl character of a BravaisType object and integer numbers
blConversionTable <- function(value)
{
 # If value is a number
 if (is.numeric(value))
 {
  if (!(value %in% 1:19)) stop("Only integers between 1 and 19 correspond to allowed Bravais types.")
  if (value == 1) rvalue <- "aP"
  if (value == 2) rvalue <- "mP"
  if (value == 3) rvalue <- "mA"
  if (value == 4) rvalue <- "mB"
  if (value == 5) rvalue <- "mC"
  if (value == 6) rvalue <- "mI"
  if (value == 7) rvalue <- "oP"
  if (value == 8) rvalue <- "oA"
  if (value == 9) rvalue <- "oB"
  if (value == 10) rvalue <- "oC" 
  if (value == 11) rvalue <- "oF"
  if (value == 12) rvalue <- "oI"
  if (value == 13) rvalue <- "tP"
  if (value == 14) rvalue <- "tI"
  if (value == 15) rvalue <- "hP"
  if (value == 16) rvalue <- "hR"
  if (value == 17) rvalue <- "cP"
  if (value == 18) rvalue <- "cF"
  if (value == 19) rvalue <- "cI"
 }

 # If value is a character
 if (!is.numeric(value))
 {
  if (!is.character(value)) stop("Input for this function is either a number or a character")
  if (value != "aP" & 
      value != "mP" & value != "mA" & value != "mB" & value != "mC" & value != "mI" &
      value != "oP" & value != "oA" & value != "oB" & value != "oC" & value != "oF" & value != "oI" &
      value != "tP" & value != "tI" &
      value != "hR" & value != "hP" &
      value != "cP" & value != "cF" & value != "cI") stop("This type of Bravais type character does not exists.")
  if (value == "aP") rvalue <- 1
  if (value == "mP") rvalue <- 2
  if (value == "mA") rvalue <- 3
  if (value == "mB") rvalue <- 4
  if (value == "mC") rvalue <- 5
  if (value == "mI") rvalue <- 6
  if (value == "oP") rvalue <- 7
  if (value == "oA") rvalue <- 8
  if (value == "oB") rvalue <- 9
  if (value == "oC") rvalue <- 10
  if (value == "oF") rvalue <- 11
  if (value == "oI") rvalue <- 12
  if (value == "tP") rvalue <- 13
  if (value == "tI") rvalue <- 14
  if (value == "hP") rvalue <- 15
  if (value == "hR") rvalue <- 16
  if (value == "cP") rvalue <- 17
  if (value == "cF") rvalue <- 18
  if (value == "cI") rvalue <- 19
  rvalue <- as.integer(rvalue)
 }

 return(rvalue)
}
