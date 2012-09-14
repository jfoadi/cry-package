# R code to implement the crystallographic ideas related to crystal lattices.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## BravaisType and Lattice classes
#
# BravaisType
setClass(
         Class="BravaisType",
         representation=representation(bl="character"),
         validity=function(object)
         {
          # Extract lattice and centering (also check on symbol validity)
          if (length(object@bl) != 0)
          {
           latt <- substr(object@bl,1,1)
           centring <- substr(object@bl,2,2)
           if (latt != "a" & latt !="m" & latt != "o" & latt != "t" & latt != "h" & latt != "c") stop("Invalid lattice type.") 
           if (centring == "S") stop("Please, specify if it is an A, B or C centring")
           if (centring != "P" & centring != "A" & centring != "B" & centring != "C" & centring != "F" & centring != "I" & centring != "R") stop("Invalid centring type")

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
          }

          # Green light: object appears to be correct
          return(TRUE)
         }
        )
#
# Lattice
setClass(
         Class="Lattice",
         representation=representation(cell="UnitCell",bl="BravaisType"),
         validity=function(object)
         {
          if (length(object@cell@a) != 0 & length(object@cell@b) != 0 & length(object@cell@c) != 0 &
              length(object@cell@alpha@ang) != 0 & length(object@cell@beta@ang) != 0 & length(object@cell@gamma@ang) != 0 &
              length(object@bl@bl) != 0)
          {
           # Extract cell parameters
           a <- object@cell@a
           b <- object@cell@b
           c <- object@cell@c
           aa <- object@cell@alpha@ang
           bb <- object@cell@beta@ang
           cc <- object@cell@gamma@ang

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
           if (latt == "h" & centring == "P")                    # Hexagonal family
           {
            if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type") 
            if (abs(cc-120) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
            if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
           }
           if (latt == "h" & centring == "R")                    # Rombohedral family
           {
            if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
            if (abs(a-b) < 0.000001)
            {
             if (abs(a-c) < 0.000001)
             {
              if (abs(aa-bb) > 0.000001 | abs(aa-cc) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
              if (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)
              {
               messaggio <- paste("Settings for this lattice are for a cubic, rather than a rombohedral system.",
                                  "Change BravaisType if you wish to carry on with this object")
               stop(messaggio)
              }
             }
             if (abs(a-c) > 0.000001)
             {
              if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 & abs(cc-120) > 0.000001) 
                  stop ("Unit cell parameters non compatible with lattice type")
             }
            }
           }
          }

          # Green light: object appears to be correct
          return(TRUE)
         }
        )

## Default methods
#
# Print
setMethod(
          f="print",
          signature="BravaisType",
          definition=function(x,...)
                     {
                      # Initial message
                      cat("*** This is a BravaisType object. ***\n")

                      # If object is empty print message
                      if (length(x@bl) == 0) cat(" The slot for parameter bl is empty.\n")

                      # Extract bl character
                      if (length(x@bl) != 0)
                      {
                       bl <- x@bl
                       cat(paste(" The slot for parameter bl contains character ",bl,".\n",sep=""))
                      }
                     }
         )
#
# Print
setMethod(
          f="print",
          signature="Lattice",
          definition=function(x,...)
                     {
                      # Initial message
                      cat("*** This is a Lattice object. ***\n")

                      # Print empty or cell values
                      cat(" The unit cell has the following parameters:\n")
                      if (length(x@cell@a) == 0) cat(" a is empty.\n")
                      if (length(x@cell@a) != 0) cat(sprintf(" a = %7.2f\n",x@cell@a))
                      if (length(x@cell@b) == 0) cat(" b is empty.\n")
                      if (length(x@cell@b) != 0) cat(sprintf(" b = %7.2f\n",x@cell@b))
                      if (length(x@cell@c) == 0) cat(" c is empty.\n")
                      if (length(x@cell@c) != 0) cat(sprintf(" c = %7.2f\n",x@cell@c))
                      if (length(x@cell@alpha@ang) == 0) cat(" alpha is empty.\n")
                      if (length(x@cell@alpha@ang) != 0) cat(sprintf(" alpha = %6.2f degrees.\n",x@cell@alpha@ang))
                      if (length(x@cell@beta@ang) == 0) cat(" beta is empty.\n")
                      if (length(x@cell@beta@ang) != 0) cat(sprintf(" beta = %6.2f degrees.\n",x@cell@beta@ang))
                      if (length(x@cell@gamma@ang) == 0) cat(" gamma is empty.\n")
                      if (length(x@cell@gamma@ang) != 0) cat(sprintf(" gamma = %6.2f degrees.\n",x@cell@gamma@ang))
                      if (length(x@bl@bl) == 0) cat(" Information on the lattice is not vailable.\n")

                      # Extract lattice information (if lattice has been defined)
                      if (length(x@bl@bl) != 0)
                      {
                       latt_info <- extractLatticeStuff(x)

                       # Print
                       stringa <- sprintf("\n The Bravais lattice is %s",latt_info[1])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal family is %s",latt_info[2])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal system is %s",latt_info[3])
                       cat(stringa)
                       stringa <- sprintf("\n The lattice system is %s\n",latt_info[4])
                       cat(stringa)
                      }
                     }
         )
#
# show
setMethod(
          f="show",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Initial message
                      cat("*** This is a BravaisType object. ***\n")

                      # If object is empty print message
                      if (length(object@bl) == 0) cat(" The slot for parameter bl is empty.\n")

                      # Extract bl character
                      if (length(object@bl) != 0)
                      {
                       bl <- object@bl
                       cat(paste(" The slot for parameter bl contains character ",bl,".\n",sep=""))
                      }
                     }
         )
#
# Show
setMethod(
          f="show",
          signature="Lattice",
          definition=function(object)
                     {
                      # Initial message
                      cat("*** This is a Lattice object. ***\n")

                      # Print empty or cell values
                      cat(" The unit cell has the following parameters:\n")
                      if (length(object@cell@a) == 0) cat(" a is empty.\n")
                      if (length(object@cell@a) != 0) cat(sprintf(" a = %7.2f\n",object@cell@a))
                      if (length(object@cell@b) == 0) cat(" b is empty.\n")
                      if (length(object@cell@b) != 0) cat(sprintf(" b = %7.2f\n",object@cell@b))
                      if (length(object@cell@c) == 0) cat(" c is empty.\n")
                      if (length(object@cell@c) != 0) cat(sprintf(" c = %7.2f\n",object@cell@c))
                      if (length(object@cell@alpha@ang) == 0) cat(" alpha is empty.\n")
                      if (length(object@cell@alpha@ang) != 0) cat(sprintf(" alpha = %6.2f degrees.\n",object@cell@alpha@ang))
                      if (length(object@cell@beta@ang) == 0) cat(" beta is empty.\n")
                      if (length(object@cell@beta@ang) != 0) cat(sprintf(" beta = %6.2f degrees.\n",object@cell@beta@ang))
                      if (length(object@cell@gamma@ang) == 0) cat(" gamma is empty.\n")
                      if (length(object@cell@gamma@ang) != 0) cat(sprintf(" gamma = %6.2f degrees.\n",object@cell@gamma@ang))
                      if (length(object@bl@bl) == 0) cat(" Information on the lattice is not vailable.\n")

                      # Extract lattice information (if lattice has been defined)
                      if (length(object@bl@bl) != 0)
                      {
                       latt_info <- extractLatticeStuff(object)

                       # Print
                       stringa <- sprintf("\n The Bravais lattice is %s",latt_info[1])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal family is %s",latt_info[2])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal system is %s",latt_info[3])
                       cat(stringa)
                       stringa <- sprintf("\n The lattice system is %s\n",latt_info[4])
                       cat(stringa)
                      }
                     }
         )

## Methods sharing generic functions with other objects
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
                      # Build list to be returned
                      lista <- list(a=object@cell@a,b=object@cell@b,c=object@cell@c,
                                    alpha=object@cell@alpha@ang,beta=object@cell@beta@ang,gamma=object@cell@gamma@ang,
                                    bl=object@bl@bl)

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

                      # Only 1 parameter for object of class "BravaisType")
                      if (n < 1 | n > 1) stop("This object does not include so many parameters")

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "bl") flist <- c(flist,list(bl=x@bl))
                       if (i[ii] == 1) flist <- c(flist,list(bl=x@bl))
                       if (i[ii] != "bl" & i[ii] != 1) stop("Parameter not included in list for this object")
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

                      # Only 1 parameter for object of class "BravaisType")
                      if (n < 1 | n > 1) stop("This object does not include so many parameters")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: bl (character)
                      wflag <- FALSE
                      if (is.list(value)) value <- value[[1]]
                      if (i == "bl") x@bl <- value
                      if (i == 1) x@bl <- value
                      if (i != "bl" & i != 1) wflag <- TRUE
                      if (wflag) warning("Requested slots does not exist")

                      if (validObject(x)) return(x)
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

                      # Only 7 parameters for object of class "Lattice")
                      if (n < 1 | n > 7) stop("This object does not include so many parameters")

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
                          ) stop("Parameter not included in list for this object")
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

                      # Only 7 parameters for object of class "Lattice")
                      if (n < 1 | n > 7) stop("This object does not include so many parameters")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Check if input is a list
                      if (n > 1 & !is.list(value)) stop("For multiple replacement, values need to be included in a list")

                      # Parameters types are: a (numeric), b (numeric), c (numeric), alpha (numeric), beta (numeric), gamma (numeric),
                      # bl (character)  
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") x@cell@a <- value[[ii]]
                       if (i[ii] == 1) x@cell@a <- value[[ii]]
                       if (i[ii] == "b") x@cell@b <- value[[ii]]
                       if (i[ii] == 2) x@cell@b <- value[[ii]]
                       if (i[ii] == "c") x@cell@c <- value[[ii]]
                       if (i[ii] == 3) x@cell@c <- value[[ii]]
                       if (i[ii] == "alpha") x@cell@alpha@ang <- value[[ii]]
                       if (i[ii] == 4) x@cell@alpha@ang <- value[[ii]]
                       if (i[ii] == "beta") x@cell@beta@ang <- value[[ii]]
                       if (i[ii] == 5) x@cell@beta@ang <- value[[ii]]
                       if (i[ii] == "gamma") x@cell@gamma@ang <- value[[ii]]
                       if (i[ii] == 6) x@cell@gamma@ang <- value[[ii]]
                       if (i[ii] == "bl") x@bl@bl <- value[[ii]]
                       if (i[ii] == 7) x@bl@bl <- value[[ii]]
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

                      if (validObject(x)) return(x)
                     }
         )

## Methods specific for BravaisType and Lattice classes
#
# Returns centering operators (shared with Symmetry class)
setMethod(
          f="getCentringOps",
          signature="Lattice",
          definition=function(object){
                                      # Centring operators are a list of 3D vectors

                                      if (length(object@bl@bl) != 0)
                                      {
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
                                      if (length(object@bl@bl) == 0) return(NULL)
                                     }
         )

## Friendly constructors and functions
#
setGeneric(
           name="bravaistype",
           def=function(bl,file,...){standardGeneric("bravaistype")}
          )
#
# Create a BravaisType object starting from a string (could be Bravais Type or file name)
setMethod(
          f="bravaistype",
          signature=c("character","missing"),
          function(bl,file,...)
          {
           # If character indicates file name call proper specification of generic function
           if (bl != "aP" & bl != "mS" & bl != "mA" & bl != "mB" & bl != "mC" & bl != "mI" & bl != "mP" &
               bl != "oP" & bl != "oS" & bl != "oA" & bl != "oB" & bl != "oC" & bl != "oI" & bl != "oF" &
               bl != "tP" & bl != "tI" &
               bl != "hP" & bl != "hR" &
               bl != "cP" & bl != "cI" & bl != "cF")
           {
            ans <- file.exists(bl)
            if (ans)
            {
             object <- createFromFile(bl,objectName="BravaisType")
             return(object)
            }
            if (!ans)
            {
             stop("Character input is neither a valid Bravais Type nor an existing file name")
            }
           }

           # Otherwise is a Bravais Type
           object <- new(Class="BravaisType",bl=bl)

           return(object)
          }
         )
#
# Create a BravaisType object starting from a file
setMethod(
          f="bravaistype",
          signature=c("missing","character"),
          function(bl,file,...)
          {
           ans <- file.exists(file)
           if (!ans)
           {
            stop("Character input is not an existing file name")
           }
           object <- createFromFile(file,objectName="BravaisType")

           return(object)
          }
         )
#
# Create a BravaisType object starting from nothing (default is cP)
setMethod(
          f="bravaistype",
          signature=c("missing","missing"),
          function(bl,file,...)
          {
           object <- new(Class="BravaisType",bl="cP")

           return(object)
          }
         )

#
setGeneric(
           name="lattice",
           def=function(cell,bl,file,cbl,...){standardGeneric("lattice")}
          )
#
# Starting from both UnitCell and BravaisType
setMethod(
          f="lattice",
          signature=c("UnitCell","BravaisType","missing","missing"),
          function(cell,bl,file,cbl,...)
          {
           object <- new(Class="Lattice",cell=cell,bl=bl)

           return(object)
          }
         )
#
# Starting from UnitCell only and assigning arbitrary, but compatible, BravaisType
setMethod(
          f="lattice",
          signature=c("UnitCell","missing","missing","missing"),
          function(cell,bl,file,cbl,...)
          {
           # Extract values from cell
           a <- cell@a
           b <- cell@b
           c <- cell@c
           aa <- cell@alpha@ang
           bb <- cell@beta@ang
           cc <- cell@gamma@ang

           # Default value for bbl is triclinic
           bbl <- "aP"

           # Cubic
           if (abs(a-b) < 0.000001 & abs(a-c) < 0.000001 & abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001) bbl <- "cP"

           # Tetragonal
           if ((abs(a-b) < 0.000001 & abs(a-c) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "tP"
           if ((abs(a-c) < 0.000001 & abs(a-b) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "tP"
           if ((abs(b-c) < 0.000001 & abs(b-a) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "tP"

           # Hexagonal
           if ((abs(a-b) < 0.000001 & abs(a-c) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-120) < 0.000001)) bbl <- "hP"
           if ((abs(a-c) < 0.000001 & abs(a-b) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-120) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "hP"
           if ((abs(b-c) < 0.000001 & abs(b-a) >= 0.000001) & (abs(aa-120) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "hP"

           # Rombohedral
           if ((abs(a-b) < 0.000001 & abs(a-c) < 0.000001) & (abs(aa-bb) < 0.000001 & abs(aa-cc) < 0.000001)) bbl <- "hR"

           # Orthorombic
           if ((abs(a-b) > 0.000001 & abs(a-c) > 0.000001 & abs(b-c) > 0.000001) &
               (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "oP"

           # Monoclinic
           if (abs(aa-120) > 0.000001 & abs(aa-90) > 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001) bbl <- "mP"
           if (abs(bb-120) > 0.000001 & abs(bb-90) > 0.000001 & abs(aa-90) < 0.000001 & abs(cc-90) < 0.000001) bbl <- "mP"
           if (abs(cc-120) > 0.000001 & abs(cc-90) > 0.000001 & abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001) bbl <- "mP"

           # Create BravaisType object
           bl <- bravaistype(bbl)

           # Create Lattice object
           object <- new(Class="Lattice",cell=cell,bl=bl)

           return(object)
          }
         )
#
# Starting from BravaisType only and assigning arbitrary, but compatible, unit cell
setMethod(
          f="lattice",
          signature=c("BravaisType","missing","missing","missing"),
          function(cell,bl,file,cbl,...)
          {
           # Extract character from BravaisType object (cell in this case, because first in the list)
           bbl <- cell@bl

           # Cubic (default 1 1 1 90 90 90)
           if (bbl == "cP" | bbl == "cI" | bbl == "cF") cella <- unitcell()

           # Hexagonal
           if (bbl == "hP") cella <- unitcell(a=1,b=1,c=2,gamma=120) 

           # Rombohedral
           if (bbl == "hR") cella <- unitcell(a=1,b=1,c=1,alpha=80,beta=80,gamma=80)

           # Tetragonal
           if (bbl == "tP" | bbl == "tI") cella <- unitcell(a=1,b=1,c=2)

           # Orthorombic
           if (bbl == "oP" | bbl == "oI" | bbl == "oF") cella <- unitcell(a=1,b=2,c=3)

           # Monoclinic
           if (bbl == "mP" | bbl == "mS" | bbl == "mI" | bbl == "mA" | bbl == "mC") cella <- unitcell(a=1,b=2,c=3,beta=110)
           if (bbl == "mB") cella <- unitcell(a=1,b=2,c=3,alpha=110)

           # Triclinic
           if (bbl == "aP") cella <- unitcell(a=1,b=2,c=3,alpha=80,beta=70,gamma=75)

           # Create lattice object
           object <- new(Class="Lattice",cell=cella,bl=cell)

           return(object)
          }
         )
#
# Starting from character. Could be either file name or lattice symbol. Cell is arbitrarily, but compatibly, assigned.
setMethod(
          f="lattice",
          signature=c("character","missing","missing","missing"),
          function(cell,bl,file,cbl,...)
          {
           # If character indicates file name call proper specification of generic function
           if (cell != "aP" & cell != "mS" & cell != "mA" & cell != "mB" & cell != "mC" & cell != "mI" & cell != "mP" &
               cell != "oP" & cell != "oS" & cell != "oA" & cell != "oB" & cell != "oC" & cell != "oI" & cell != "oF" &
               cell != "tP" & cell != "tI" &
               cell != "hP" & cell != "hR" &
               cell != "cP" & cell != "cI" & cell != "cF")
           {
            object <- createFromFile(cell,objectName="Lattice")

            return(object)
           }

           # Character indicates BravaisType object (cell in this case, because first in the list)
           bbl <- cell

           # Cubic (default 1 1 1 90 90 90)
           if (bbl == "cP" | bbl == "cI" | bbl == "cF") cella <- unitcell()

           # Hexagonal
           if (bbl == "hP") cella <- unitcell(a=1,b=1,c=2,gamma=120) 

           # Rombohedral
           if (bbl == "hR") cella <- unitcell(a=1,b=1,c=1,alpha=80,beta=80,gamma=80)

           # Tetragonal
           if (bbl == "tP" | bbl == "tI") cella <- unitcell(a=1,b=1,c=2)

           # Orthorombic
           if (bbl == "oP" | bbl == "oI" | bbl == "oF") cella <- unitcell(a=1,b=2,c=3)

           # Monoclinic
           if (bbl == "mP" | bbl == "mS" | bbl == "mI" | bbl == "mA" | bbl == "mC") cella <- unitcell(a=1,b=2,c=3,beta=110)
           if (bbl == "mB") cella <- unitcell(a=1,b=2,c=3,alpha=110)

           # Triclinic
           if (bbl == "aP") cella <- unitcell(a=1,b=2,c=3,alpha=80,beta=70,gamma=75)

           # Create BravaisType object
           bl <- bravaistype(bbl)

           # Create Lattice object
           object <- new(Class="Lattice",cell=cella,bl=bl)

           return(object)
          }
         )
#
# Starting from BravaisType only and assigning arbitrary, but compatible, unit cell
setMethod(
          f="lattice",
          signature=c("missing","BravaisType","missing","missing"),
          function(cell,bl,file,cbl,...)
          {
           # Extract character from BravaisType object
           bbl <- bl@bl

           # Cubic (default 1 1 1 90 90 90)
           if (bbl == "cP" | bbl == "cI" | bbl == "cF") cella <- unitcell()

           # Hexagonal
           if (bbl == "hP") cella <- unitcell(a=1,b=1,c=2,gamma=120) 

           # Rombohedral
           if (bbl == "hR") cella <- unitcell(a=1,b=1,c=1,alpha=80,beta=80,gamma=80)

           # Tetragonal
           if (bbl == "tP" | bbl == "tI") cella <- unitcell(a=1,b=1,c=2)

           # Orthorombic
           if (bbl == "oP" | bbl == "oI" | bbl == "oF") cella <- unitcell(a=1,b=2,c=3)

           # Monoclinic
           if (bbl == "mP" | bbl == "mS" | bbl == "mI" | bbl == "mA" | bbl == "mC") cella <- unitcell(a=1,b=2,c=3,beta=110)
           if (bbl == "mB") cella <- unitcell(a=1,b=2,c=3,alpha=110)

           # Triclinic
           if (bbl == "aP") cella <- unitcell(a=1,b=2,c=3,alpha=80,beta=70,gamma=75)

           # Create lattice object
           object <- new(Class="Lattice",cell=cella,bl=bl)

           return(object)
          }
         )
#
# Starting from character only (indicating BravaisType) and assigning arbitrary, but compatible, unit cell
setMethod(
          f="lattice",
          signature=c("missing","missing","missing","character"),
          function(cell,bl,file,cbl,...)
          {
           # Character indicates BravaisType object (cell in this case, because first in the list)
           bbl <- cbl

           # Cubic (default 1 1 1 90 90 90)
           if (bbl == "cP" | bbl == "cI" | bbl == "cF") cella <- unitcell()

           # Hexagonal
           if (bbl == "hP") cella <- unitcell(a=1,b=1,c=2,gamma=120) 

           # Rombohedral
           if (bbl == "hR") cella <- unitcell(a=1,b=1,c=1,alpha=80,beta=80,gamma=80)

           # Tetragonal
           if (bbl == "tP" | bbl == "tI") cella <- unitcell(a=1,b=1,c=2)

           # Orthorombic
           if (bbl == "oP" | bbl == "oI" | bbl == "oF") cella <- unitcell(a=1,b=2,c=3)

           # Monoclinic
           if (bbl == "mP" | bbl == "mS" | bbl == "mI" | bbl == "mA" | bbl == "mC") cella <- unitcell(a=1,b=2,c=3,beta=110)
           if (bbl == "mB") cella <- unitcell(a=1,b=2,c=3,alpha=110)

           # Triclinic
           if (bbl == "aP") cella <- unitcell(a=1,b=2,c=3,alpha=80,beta=70,gamma=75)

           # Create BravaisType object
           bl <- bravaistype(bbl)

           # Create Lattice object
           object <- new(Class="Lattice",cell=cella,bl=bl)

           return(object)
          }
         )
#
# Starting from a file
setMethod(
          f="lattice",
          signature=c("missing","missing","character","missing"),
          function(cell,bl,file,cbl,...)
          {
           object <- createFromFile(file,objectName="Lattice")

           return(object)
          }
         )
#
# Starting from nothing (cubic lattice with 1 1 1 cell sides)
setMethod(
          f="lattice",
          signature=c("missing","missing","missing","missing"),
          function(cell,bl,file,cbl,...)
          {
           bla <- bravaistype()
           cella <- unitcell()

           object <- new(Class="Lattice",cell=cella,bl=bla)

           return(object)
          }
         )

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
