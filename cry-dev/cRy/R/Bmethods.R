###################################################################
###################################################################
### This module is part of the cRy package for crystallography. ###
### Authors: J. Foadi & D. G. Waterman                          ###
###                                                             ###
### Methods collated later in the code.                         ###
###################################################################
###################################################################



######### For Angle class
#
# Degrees to radians conversion
setMethod(
          f="degToRad",
          signature="Angle",
          definition=function(object)
                     {
                      if (length(object@rad_flag) != 0)
                      {
                       if (object@rad_flag)
                       {
                        warning("Object of class Angle is already expressed in radians")
                        return(object)
                       }
                       if (length(object@ang) != 0)
                       {
                        object@ang <- (object@ang*pi/180)%%(2*pi)
                        object@rad_flag <- TRUE
                        return(object)
                       }
                       if (length(object@ang) == 0)
                       {
                        object@rad_flag <- TRUE
                        return(object)
                       }
                      }
                      if (length(object@rad_flag) == 0)
                      {
                       warning("This object of class Angle has not been expressed neither in degrees nor in radians.")
                       return(object)
                      }
                     }
         )
#
# Radians to degrees conversion
setGeneric(
           name="radToDeg",
           def=function(object){standardGeneric("radToDeg")}
          )
setMethod(
          f="radToDeg",
          signature="Angle",
          definition=function(object)
                     {
                      if (length(object@rad_flag) != 0)
                      {
                       if (!object@rad_flag)
                       {
                        warning("Object of class Angle is already expressed in degrees")
                        return(object)
                       }
                       if (length(object@ang) != 0)
                       {
                        object@ang <- (object@ang*180/pi)%%360
                        object@rad_flag <- FALSE
                        return(object)
                       }
                       if (length(object@ang) == 0)
                       {
                        object@rad_flag <- FALSE
                        return(object)
                       }
                      }
                      if (length(object@rad_flag) == 0)
                      {
                       warning("This object of class Angle has not been expressed neither in degrees nor in radians.")
                       return(object)
                      }
                     }
         )



######### For UnitCell class
#
# Compute reciprocal cell starting from direct cell
setMethod(
          f="computeReciprocalUnitCell",
          signature="UnitCell",
          definition=function(object)
                     {
                      # Return empty object if one or more slots are empty
                      if (length(object@a) == 0 | length(object@b) == 0 | length(object@c) == 0 |
                          length(object@alpha@ang) == 0 | length(object@beta@ang) == 0 | length(object@gamma@ang) == 0)
                      {
                       return(new("ReciprocalUnitCell"))
                      }

                      # Compute the 6-parameters list
                      lista <- .dcl(object@a,object@b,object@c,object@alpha@ang,object@beta@ang,object@gamma@ang)

                      # Calculate all 6 reciprocal unit cell parameters
                      ar <- lista[[7]]
                      br <- lista[[8]]
                      cr <- lista[[9]]
                      alphar <- atan2(lista[[10]],lista[[13]])*180/pi
                      betar  <- atan2(lista[[11]],lista[[14]])*180/pi
                      gammar <- atan2(lista[[12]],lista[[15]])*180/pi

                      # Create ReciprocalUnitCell object with above values
                      new_object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,
                                                     alphar=new(Class="Angle",ang=alphar,rad_flag=FALSE),
                                                     betar=new(Class="Angle",ang=betar,rad_flag=FALSE),
                                                     gammar=new(Class="Angle",ang=gammar,rad_flag=FALSE))

                      return(new_object)
                     }
         )
#
# Compute cell volume (valid for both direct and reciprocal cells)
setMethod(
          f="computeCellVolume",
          signature="UnitCell",
          definition=function(object)
                     {
                      # If object is empty return NULL
                      if (length(object@a) == 0 | length(object@b) == 0 | length(object@c) == 0 |
                          length(object@alpha@ang) == 0 | length(object@beta@ang) == 0 | length(object@gamma@ang) == 0)
                      {
                       return(NULL)
                      }

                      # Extract cell parameters
                      a <- object@a
                      b <- object@b
                      c <- object@c
                      alpha <- object@alpha@ang
                      beta <- object@beta@ang
                      gamma <- object@gamma@ang

                      # Compute the 6-parameters list
                      lista <- .dcl(a,b,c,alpha,beta,gamma)

                      return(lista[[16]])
                     }
         )



######### For ReciprocalUnitCell class
#
# Compute direct cell starting from reciprocal cell
setMethod(
          f="computeUnitCell",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # Return empty object if one or more slots are empty
                      if (length(object@ar) == 0 | length(object@br) == 0 | length(object@cr) == 0 |
                          length(object@alphar@ang) == 0 | length(object@betar@ang) == 0 | length(object@gammar@ang) == 0)
                      {
                       return(new("UnitCell"))
                      }

                      # Now compute the 6-parameters list
                      lista <- .dcl(object@ar,object@br,object@cr,object@alphar@ang,object@betar@ang,object@gammar@ang)

                      # Calculate all 6 unit cell parameters
                      a <- lista[[7]]
                      b <- lista[[8]]
                      c <- lista[[9]]
                      alpha <- atan2(lista[[10]],lista[[13]])*180/pi
                      beta  <- atan2(lista[[11]],lista[[14]])*180/pi
                      gamma <- atan2(lista[[12]],lista[[15]])*180/pi

                      # Create UnitCell object with above values
                      new_object <- new(Class="UnitCell",a=a,b=b,c=c,
                                                     alpha=new(Class="Angle",ang=alpha,rad_flag=FALSE),
                                                     beta=new(Class="Angle",ang=beta,rad_flag=FALSE),
                                                     gamma=new(Class="Angle",ang=gamma,rad_flag=FALSE))

                      return(new_object)
                     }
         )
#
# Compute cell volume (valid for both direct and reciprocal cells)
setMethod(
          f="computeCellVolume",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # If object is empty return NULL
                      if (length(object@ar) == 0 | length(object@br) == 0 | length(object@cr) == 0 |
                          length(object@alphar@ang) == 0 | length(object@betar@ang) == 0 | length(object@gammar@ang) == 0)
                      {
                       return(NULL)
                      }

                      # Extract cell parameters
                      a <- object@ar
                      b <- object@br
                      c <- object@cr
                      alpha <- object@alphar@ang
                      beta <- object@betar@ang
                      gamma <- object@gammar@ang

                      # Compute the 6-parameters list
                      lista <- .dcl(a,b,c,alpha,beta,gamma)

                      return(lista[[16]])
                     }
         )



######### For BravaisType class
#
# Return centring operators. Applies to BravaisType, Lattice and Symmetry classes
setMethod(
          f="getCentringOps",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Centring operators are a list of 3D vectors

                      if (length(object@bl) != 0)
                      {
                       # P
                       if (substr(object@bl,2,2) == "P") lop <- list(c(0,0,0))

                       # A
                       if (substr(object@bl,2,2) == "A") lop <- list(c(0,0,0),c(0,0.5,0.5))

                       # B
                       if (substr(object@bl,2,2) == "B") lop <- list(c(0,0,0),c(0.5,0,0.5))

                       # C
                       if (substr(object@bl,2,2) == "C") lop <- list(c(0,0,0),c(0.5,0.5,0))

                       # I
                       if (substr(object@bl,2,2) == "I") lop <- list(c(0,0,0),c(0.5,0.5,0.5))

                       # F
                       if (substr(object@bl,2,2) == "F") lop <- list(c(0,0,0),c(0,0.5,0.5),c(0.5,0,0.5),c(0.5,0.5,0))

                       # H or R (H does not appear in the Pearson symbol. It means a centred cell. But here we have an anbiguity and
                       # we will use no centring, i.e. will assume rombohedral. Things will be defined with Lattice objects
                       if (substr(object@bl,2,2) == "R") lop <- list(c(0,0,0))

                       return(lop)
                      }
                      if (length(object@bl) == 0) return(NULL)
                     }
         )



######### For Lattice class
#
# Return centring operators. Applies to BravaisType, Lattice and Symmetry classes
setMethod(
          f="getCentringOps",
          signature="Lattice",
          definition=function(object)
                     {
                      # Cell parameters are needed to discriminate between rombohedral (no centring) and hexagonal (centring) cells.
                      cell_flag <- FALSE
                      if (length(object@cell) != 0) 
                      {
                       cell_flag <- TRUE
                       a <- object@cell@a
                       b <- object@cell@b
                       c <- object@cell@c
                       alpha <- object@cell@alpha@ang
                       beta <- object@cell@beta@ang
                       gamma <- object@cell@gamma@ang
                      }

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

                       # H and R
                       if (substr(object@bl@bl,2,2) == "R") 
                       {
                        if (abs(a-b) < 0.000001 & abs(a-c) < 0.000001 & abs(alpha-beta) < 0.000001 & abs(alpha-gamma) < 0.000001) 
                        {
                         lop <- list(c(0,0,0))
                        }
                        else
                        {
                         lop <- list(c(0,0,0),c(2/3,1/3,1/3),c(1/3,2/3,2/3))
                        }
                       }

                       return(lop)
                      }
                      if (length(object@bl@bl) == 0) return(NULL)
                     }
         )



######### For Symmetry class
#
# Return centring operators. Applies to BravaisType, Lattice and Symmetry classes
setMethod(
          f="getCentringOps",
          signature="Symmetry",
          definition=function(object)
                     {
                      # If Symmetry object is empty, return NULL
                      if (length(object@sym_xHM) == 0)
                      {
                       return(NULL)
                      }

                      # Extract vectors from Symmetry object
                      lop <- .syminfo_to_matrix_list(object@sym_xHM)$C

                      return(lop)
                     }
         )



######### For UnitCell and Symmetry classes together
#
# Check symmetry is compatible with unit cell parameters
setMethod(
          f="checkSymmetryWithCell",
          signature=c("UnitCell","Symmetry"),
          definition=function(obUnitCell,obSymmetry)
                     {
                      # Space group number and settings
                      sgn <- .getSymmetryNumber(obSymmetry)

                      # Extract crystal system associated with symmetry
                      cr_sys <- .crystal_system(sgn[1])

                      # Extract centring associated with space group
                      sl <- substr(obSymmetry@sym_xHM,1,1)

                      # Build BravaisType object
                      if (cr_sys == "TRICLINIC") fl <- "a"
                      if (cr_sys == "MONOCLINIC") fl <- "m"
                      if (cr_sys == "ORTHOROMBIC") fl <- "o"
                      if (cr_sys == "TETRAGONAL") fl <- "t"
                      if (cr_sys == "CUBIC") fl <- "c"
                      if (cr_sys == "HEXAGONAL") fl <- "h"
                      if (cr_sys == "TRIGONAL") fl <- "h"
                      bl <- bravaistype(paste(fl,sl,sep=""))

                      # Now build Lattice object using UnitCell and BravaisType objects. If incompatible, check will fail
                      tmplatt <- lattice(obUnitCell,bl)

                      # Other checks, specific to various settings
                      cpar <- c(obUnitCell@a,obUnitCell@b,obUnitCell@c,obUnitCell@alpha@ang,obUnitCell@beta@ang,obUnitCell@gamma@ang) 
                      a <- cpar[1]
                      b <- cpar[2]
                      c <- cpar[3]
                      alpha <- cpar[4]
                      beta <- cpar[5]
                      gamma <- cpar[6]

                      # Monoclinic
                      if (sgn[1] %in% 3:15)
                      {
                       idx <- .getMonoclinicConstraints(sgn)
                       if (abs(cpar[idx]-90) < 0.000001) stop("Unit cell parameters not compatible with this space group")
                      }

                      # Trigonal
                      if (sgn[1] %in% c(146,148,155,160,161,166,167))
                      {
                       if (sgn[2] == 1)
                       {
                        if (abs(a-c) < 0.000001 | abs(alpha-gamma) < 0.000001) stop("Unit cell parameters not compatible with this space group")
                       }
                       if (sgn[2] == 2)
                       {
                        if (abs(a-b) > 0.000001 | abs(a-c) > 0.000001) stop("Unit cell parameters not compatible with this space group")
                        if (abs(alpha-beta) > 0.000001 | abs(alpha-gamma) > 0.000001) stop("Unit cell parameters not compatible with this space group")
                       }
                      }
                     }
         )
