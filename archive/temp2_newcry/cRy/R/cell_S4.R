# R code to implement the crystallographic ideas related to unit cell.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## UnitCell and ReciprocalUnitCell classes
#
# UnitCell (angles always stored in degrees)
setClass(
         Class="UnitCell",
         representation=representation(a="numeric",b="numeric",c="numeric",alpha="Angle",beta="Angle",gamma="Angle"),
         validity=function(object)
                  {
                   # Angles are always stored in degrees
                   if (length(object@alpha@ang) != 0)
                   {
                    if (object@alpha@rad_flag) stop("Cell angles always stored in degrees for object of class UnitCell.")
                   }
                   if (length(object@beta@ang) != 0)
                   {
                    if (object@beta@rad_flag) stop("Cell angles always stored in degrees for object of class UnitCell.")
                   }
                   if (length(object@gamma@ang) != 0)
                   {
                    if (object@gamma@rad_flag) stop("Cell angles always stored in degrees for object of class UnitCell.")
                   }

                   # Cell sides can be any number different from zero or infinity
                   if (length(object@a) != 0)
                   {
                    if (object@a <= 0 | !is.finite(object@a)) stop("Cell side a can be neither zero nor infinity.")
                   }
                   if (length(object@b) != 0)
                   {
                    if (object@b <= 0 | !is.finite(object@b)) stop("Cell side b can be neither zero nor infinity.")
                   }
                   if (length(object@c) != 0)
                   {
                    if (object@c <= 0 | !is.finite(object@c)) stop("Cell side c can be neither zero nor infinity.")
                   }

                   # Angles are objects of class Angles, so they all have values between 0 and 180 degrees. But combinations
                   # of alpha beta and gamma need to obey certain rules (see J. Foadi & G. Evans (2011))
                   if (length(object@alpha@ang) != 0 & length(object@beta@ang) != 0 & length(object@gamma@ang) != 0)
                   {
                    aa <- object@alpha@ang
                    bb <- object@beta@ang
                    cc <- object@gamma@ang
                    ss <- aa+bb+cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                    ss <- aa+bb-cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                    ss <- aa-bb+cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                    ss <- -aa+bb+cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                   }

                   # Green light: object appears to be correct
                   return(TRUE)
                  }
        )
#
# ReciprocalUnitCell (angles always stored in degrees)
setClass(
         Class="ReciprocalUnitCell",
         representation=representation(ar="numeric",br="numeric",cr="numeric",alphar="Angle",betar="Angle",gammar="Angle"),
         validity=function(object)
                  {
                   # Angles are always stored in degrees
                   if (length(object@alphar@ang) != 0)
                   {
                    if (object@alphar@rad_flag) stop("Cell angles always stored in degrees for object of class ReciprocalUnitCell.")
                   }
                   if (length(object@betar@ang) != 0)
                   {
                    if (object@betar@rad_flag) stop("Cell angles always stored in degrees for object of class ReciprocalUnitCell.")
                   }
                   if (length(object@gammar@ang) != 0)
                   {
                    if (object@gammar@rad_flag) stop("Cell angles always stored in degrees for object of class ReciprocalUnitCell.")
                   }

                   # Cell sides can be any number different from zero or infinity
                   if (length(object@ar) != 0)
                   {
                    if (object@ar <= 0 | !is.finite(object@ar)) stop("Cell side ar can be neither zero nor infinity.")
                   } 
                   if (length(object@br) != 0)
                   {
                    if (object@br <= 0 | !is.finite(object@br)) stop("Cell side br can be neither zero nor infinity.")
                   } 
                   if (length(object@cr) != 0)
                   {
                    if (object@cr <= 0 | !is.finite(object@cr)) stop("Cell side cr can be neither zero nor infinity.")
                   } 

                   # Angles are objects of class Angles, so they all have values between 0 and 180 degrees. But combinations
                   # of alphar betar and gammar need to obey certain rules
                   if (length(object@alphar@ang) != 0 & length(object@betar@ang) != 0 & length(object@gammar@ang) != 0)
                   {
                    aar <- object@alphar@ang
                    bbr <- object@betar@ang
                    ccr <- object@gammar@ang
                    ssr <- aar+bbr+ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                    ssr <- aar+bbr-ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                    ssr <- aar-bbr+ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                    ssr <- -aar+bbr+ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocalunit cell with these angles cannot exist")
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
          signature="UnitCell",
          definition=function(x,...)
                     {
                      cat("*** This is a UnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(x@a) == 0) cat(" Slot for parameter a is empty.\n")
                      if (length(x@a) != 0) cat(sprintf(" Parameter a = %8.3f\n",x@a))
                      if (length(x@b) == 0) cat(" Slot for parameter b is empty.\n")
                      if (length(x@b) != 0) cat(sprintf(" Parameter b = %8.3f\n",x@b))
                      if (length(x@c) == 0) cat(" Slot for parameter c is empty.\n")
                      if (length(x@c) != 0) cat(sprintf(" Parameter c = %8.3f\n",x@c))
                      if (length(x@alpha@ang) == 0) cat(" Slot for parameter alpha is empty.\n")
                      if (length(x@alpha@ang) != 0) cat(sprintf(" Parameter alpha = %6.2f degrees\n",x@alpha@ang))
                      if (length(x@beta@ang) == 0) cat(" Slot for parameter beta is empty.\n")
                      if (length(x@beta@ang) != 0) cat(sprintf(" Parameter beta = %6.2f degrees\n",x@beta@ang))
                      if (length(x@gamma@ang) == 0) cat(" Slot for parameter gamma is empty.\n")
                      if (length(x@gamma@ang) != 0) cat(sprintf(" Parameter gamma = %6.2f degrees\n",x@gamma@ang))
                     }
         )
#
# Print
setMethod(
          f="print",
          signature="ReciprocalUnitCell",
          definition=function(x,...)
                     {
                      cat("*** This is a ReciprocalUnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(x@ar) == 0) cat(" Slot for parameter ar is empty.\n")
                      if (length(x@ar) != 0) cat(sprintf(" Parameter ar = %8.3f\n",x@ar))
                      if (length(x@br) == 0) cat(" Slot for parameter br is empty.\n")
                      if (length(x@br) != 0) cat(sprintf(" Parameter br = %8.3f\n",x@br))
                      if (length(x@cr) == 0) cat(" Slot for parameter cr is empty.\n")
                      if (length(x@cr) != 0) cat(sprintf(" Parameter cr = %8.3f\n",x@cr))
                      if (length(x@alphar@ang) == 0) cat(" Slot for parameter alphar is empty.\n") 
                      if (length(x@alphar@ang) != 0) cat(sprintf(" Parameter alphar = %6.2f degrees\n",x@alphar@ang)) 
                      if (length(x@betar@ang) == 0) cat(" Slot for parameter betar is empty.\n") 
                      if (length(x@betar@ang) != 0) cat(sprintf(" Parameter betar = %6.2f degrees\n",x@betar@ang)) 
                      if (length(x@gammar@ang) == 0) cat(" Slot for parameter gammar is empty.\n") 
                      if (length(x@gammar@ang) != 0) cat(sprintf(" Parameter gammar = %6.2f degrees\n",x@gammar@ang)) 
                     }
         )
#
# Show
setMethod(
          f="show",
          signature="UnitCell",
          definition=function(object)
                     {
                      cat("*** This is a UnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(object@a) == 0) cat(" Slot for parameter a is empty.\n")
                      if (length(object@a) != 0) cat(sprintf(" Parameter a = %8.3f\n",object@a))
                      if (length(object@b) == 0) cat(" Slot for parameter b is empty.\n")
                      if (length(object@b) != 0) cat(sprintf(" Parameter b = %8.3f\n",object@b))
                      if (length(object@c) == 0) cat(" Slot for parameter c is empty.\n")
                      if (length(object@c) != 0) cat(sprintf(" Parameter c = %8.3f\n",object@c))
                      if (length(object@alpha@ang) == 0) cat(" Slot for parameter alpha is empty.\n")
                      if (length(object@alpha@ang) != 0) cat(sprintf(" Parameter alpha = %6.2f degrees\n",object@alpha@ang))
                      if (length(object@beta@ang) == 0) cat(" Slot for parameter beta is empty.\n")
                      if (length(object@beta@ang) != 0) cat(sprintf(" Parameter beta = %6.2f degrees\n",object@beta@ang))
                      if (length(object@gamma@ang) == 0) cat(" Slot for parameter gamma is empty.\n")
                      if (length(object@gamma@ang) != 0) cat(sprintf(" Parameter gamma = %6.2f degrees\n",object@gamma@ang))
                     }
         )
#
# Show
setMethod(
          f="show",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      cat("*** This is a ReciprocalUnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(object@ar) == 0) cat(" Slot for parameter ar is empty.\n")
                      if (length(object@ar) != 0) cat(sprintf(" Parameter ar = %8.3f\n",object@ar))
                      if (length(object@br) == 0) cat(" Slot for parameter br is empty.\n")
                      if (length(object@br) != 0) cat(sprintf(" Parameter br = %8.3f\n",object@br))
                      if (length(object@cr) == 0) cat(" Slot for parameter cr is empty.\n")
                      if (length(object@cr) != 0) cat(sprintf(" Parameter cr = %8.3f\n",object@cr))
                      if (length(object@alphar@ang) == 0) cat(" Slot for parameter alphar is empty.\n") 
                      if (length(object@alphar@ang) != 0) cat(sprintf(" Parameter alphar = %6.2f degrees\n",object@alphar@ang)) 
                      if (length(object@betar@ang) == 0) cat(" Slot for parameter betar is empty.\n") 
                      if (length(object@betar@ang) != 0) cat(sprintf(" Parameter betar = %6.2f degrees\n",object@betar@ang)) 
                      if (length(object@gammar@ang) == 0) cat(" Slot for parameter gammar is empty.\n") 
                      if (length(object@gammar@ang) != 0) cat(sprintf(" Parameter gammar = %6.2f degrees\n",object@gammar@ang)) 
                     }
         )

## Non default methods corresponding to generic functions common to many cRy methods
#
# Extract UnitCell slots and return them in a list
setMethod(
          f="getFields",
          signature="UnitCell",
          definition=function(object)
                     {
                      # Extract slots
                      a <- object@a
                      b <- object@b
                      c <- object@c
                      alpha <- object@alpha
                      beta <- object@beta
                      gamma <- object@gamma
                      lista <- list(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

                      return(lista)
                     } 
         )
#
# Extract ReciprocalUnitCell slots and return them in a list
setMethod(
          f="getFields",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # Extract slots
                      ar <- object@ar
                      br <- object@br
                      cr <- object@cr
                      alphar <- object@alphar
                      betar <- object@betar
                      gammar <- object@gammar
                      lista <- list(ar=ar,br=br,cr=cr,alphar=alphar,betar=betar,gammar=gammar)

                      return(lista)
                     } 
         )
#
# Extract UnitCell parameters
setMethod(
          f="getParameters",
          signature="UnitCell",
          definition=function(object)
                     {
                      # Extract parameters
                      a <- object@a
                      b <- object@b
                      c <- object@c
                      alpha <- object@alpha@ang
                      beta <- object@beta@ang
                      gamma <- object@gamma@ang
                      lista <- list(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

                      return(lista)
                     } 
         )
#
# Extract ReciprocalUnitCell parameters
setMethod(
          f="getParameters",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # Extract parameters
                      ar <- object@ar
                      br <- object@br
                      cr <- object@cr
                      alphar <- object@alphar@ang
                      betar <- object@betar@ang
                      gammar <- object@gammar@ang
                      lista <- list(ar=ar,br=br,cr=cr,alphar=alphar,betar=betar,gammar=gammar)

                      return(lista)
                     } 
         )
#
# Extract UnitCell data (none)
setMethod(
          f="getData",
          signature="UnitCell",
          definition=function(object)
                     {
                      # No data are contained in an object of class UnitCell, so this function returns NULL

                      return(NULL)
                     } 
         )
#
# Extract ReciprocalUnitCell data (none)
setMethod(
          f="getData",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # No data are contained in an object of class ReciprocalUnitCell, so this function returns NULL

                      return(NULL)
                     } 
         )
#
# Another way to extract UnitCell slots (similar to getFields)
setMethod(
          f="[",
          signature="UnitCell",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") flist <- c(flist,list(a=x@a))
                       if (i[ii] == "b") flist <- c(flist,list(b=x@b))
                       if (i[ii] == "c") flist <- c(flist,list(c=x@c))
                       if (i[ii] == "alpha") flist <- c(flist,list(alpha=x@alpha@ang))
                       if (i[ii] == "beta") flist <- c(flist,list(beta=x@beta@ang))
                       if (i[ii] == "gamma") flist <- c(flist,list(gamma=x@gamma@ang))
                       if (i[ii] == 1) flist <- c(flist,list(a=x@a))
                       if (i[ii] == 2) flist <- c(flist,list(b=x@b))
                       if (i[ii] == 3) flist <- c(flist,list(c=x@c))
                       if (i[ii] == 4) flist <- c(flist,list(alpha=x@alpha@ang))
                       if (i[ii] == 5) flist <- c(flist,list(beta=x@beta@ang))
                       if (i[ii] == 6) flist <- c(flist,list(gamma=x@gamma@ang))
                       if (i[ii] != "a" & i[ii] != "b" & i[ii] != "c" &
                           i[ii] != "alpha" & i[ii] != "beta" & i[ii] != "gamma" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6)
                           stop("Slot not included in object range")
                      }

                      return(flist)
                     }
         )
#
# Change parameters values in object of class UnitCell
setReplaceMethod(
          f="[",
          signature="UnitCell",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # If n > 1 and value is not a list, stop
                      if (n > 1 & !is.list(value)) stop("For multiple replacement, replacing values need to be included in a list")

                      # Valid slots are: a (numeric), b (numeric), c (numeric), alpha (Angle), beta (Angle) gamma (Angle)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") x@a <- value[[ii]]
                       if (i[ii] == 1) x@a <- value[[ii]]
                       if (i[ii] == "b") x@b <- value[[ii]]
                       if (i[ii] == 2) x@b <- value[[ii]]
                       if (i[ii] == "c") x@c <- value[[ii]]
                       if (i[ii] == 3) x@c <- value[[ii]]
                       if (i[ii] == "alpha") x@alpha <- angle(value[[ii]])
                       if (i[ii] == 4) x@alpha <- angle(value[[ii]])
                       if (i[ii] == "beta") x@beta <- angle(value[[ii]])
                       if (i[ii] == 5) x@beta <- angle(value[[ii]])
                       if (i[ii] == "gamma") x@gamma <- angle(value[[ii]])
                       if (i[ii] == 6) x@gamma <- angle(value[[ii]])
                       if (i[ii] != "a" & i[ii] != "b" & i[ii] != "c" &
                           i[ii] != "alpha" & i[ii] != "beta" & i[ii] != "gamma" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      if (validObject(x)) return(x)
                     }
         )
#
# Another way to extract ReciprocalUnitCell slots (similar to getFields)
setMethod(
          f="[",
          signature="ReciprocalUnitCell",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ar") flist <- c(flist,list(ar=x@ar))
                       if (i[ii] == "br") flist <- c(flist,list(br=x@br))
                       if (i[ii] == "cr") flist <- c(flist,list(cr=x@cr))
                       if (i[ii] == "alphar") flist <- c(flist,list(alphar=x@alphar@ang))
                       if (i[ii] == "betar") flist <- c(flist,list(betar=x@betar@ang))
                       if (i[ii] == "gammar") flist <- c(flist,list(gammar=x@gammar@ang))
                       if (i[ii] == 1) flist <- c(flist,list(ar=x@ar))
                       if (i[ii] == 2) flist <- c(flist,list(br=x@br))
                       if (i[ii] == 3) flist <- c(flist,list(cr=x@cr))
                       if (i[ii] == 4) flist <- c(flist,list(alphar=x@alphar@ang))
                       if (i[ii] == 5) flist <- c(flist,list(betar=x@betar@ang))
                       if (i[ii] == 6) flist <- c(flist,list(gammar=x@gammar@ang))
                       if (i[ii] != "ar" & i[ii] != "br" & i[ii] != "cr" &
                           i[ii] != "alphar" & i[ii] != "betar" & i[ii] != "gammar" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6)
                           stop("Slot not included in object range")
                      }

                      return(flist)
                     }
         )
#
# Change parameters values in object of class ReciprocalUnitCell
setReplaceMethod(
          f="[",
          signature="ReciprocalUnitCell",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # If n > 1 and value is not a list, stop
                      if (n > 1 & !is.list(value)) stop("For multiple replacement, replacing values need to be included in a list")

                      # Valid slots are: ar (numeric), br (numeric), cr (numeric), alphar (Angle), betar (Angle) gammar (Angle)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ar") x@ar <- value[[ii]]
                       if (i[ii] == 1) x@ar <- value[[ii]]
                       if (i[ii] == "br") x@br <- value[[ii]]
                       if (i[ii] == 2) x@br <- value[[ii]]
                       if (i[ii] == "cr") x@cr <- value[[ii]]
                       if (i[ii] == 3) x@cr <- value[[ii]]
                       if (i[ii] == "alphar") x@alphar <- angle(value[[ii]])
                       if (i[ii] == 4) x@alphar <- angle(value[[ii]])
                       if (i[ii] == "betar") x@betar <- angle(value[[ii]])
                       if (i[ii] == 5) x@betar <- angle(value[[ii]])
                       if (i[ii] == "gammar") x@gammar <- angle(value[[ii]])
                       if (i[ii] == 6) x@gammar <- angle(value[[ii]])
                       if (i[ii] != "ar" & i[ii] != "br" & i[ii] != "cr" &
                           i[ii] != "alphar" & i[ii] != "betar" & i[ii] != "gammar" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      if (validObject(x)) return(x)
                     }
         )

## Non default methods corresponding to generic functions common only to Angle methods
#
# Compute reciprocal cell starting from direct cell
setGeneric(
           name="computeReciprocalUnitCell",
           def=function(object){standardGeneric("computeReciprocalUnitCell")}
          )
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
# Compute direct cell starting from reciprocal cell
setGeneric(
           name="computeUnitCell",
           def=function(object){standardGeneric("computeUnitCell")}
          )
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
setGeneric(
           name="computeCellVolume",
           def=function(object){standardGeneric("computeCellVolume")}
          )
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
#
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

## Friendly constructors and functions
#
setGeneric(
           name="unitcell",
           def=function(a,b,c,alpha,beta,gamma,file,...){standardGeneric("unitcell")}
          )
# All 6 parameters given (no default assignment)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","numeric","numeric","numeric","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(alpha)
           bb <- angle(beta)
           cc <- angle(gamma)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Three cell sides and one angle (beta: default)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","numeric","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(alpha)
           bb <- angle(90)
           cc <- angle(90)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Three cell sides and one angle (beta)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","missing","numeric","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(beta)
           cc <- angle(90)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Three cell sides and one angle (gamma)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","missing","missing","numeric","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(gamma)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Three cell sides and no angles (orthorombic, tetragonal or cubic)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - a)
setMethod(
          f="unitcell",
          signature=c("numeric","missing","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           b <- a
           c <- a

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - b)
setMethod(
          f="unitcell",
          signature=c("missing","numeric","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           a <- b
           c <- b

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - c)
setMethod(
          f="unitcell",
          signature=c("missing","missing","numeric","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           a <- c
           b <- c

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# No cell sides and no angles (cubic with default side 1)
setMethod(
          f="unitcell",
          signature=c("missing","missing","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           a <- 1
           b <- 1
           c <- 1

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Values taken from file (two methods, with and without "file=")
setMethod(
          f="unitcell",
          signature=c("character","missing","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # In this case a is really containing a file name
           file <- a

           # Use createFromFile function for creating object
           object <- createFromFile(fileName=file,objectName="UnitCell")

           return(object)
          }
         )
setMethod(
          f="unitcell",
          signature=c("missing","missing","missing","missing","missing","missing","character"),
          function(a,b,c,alpha,beta,gamma,file,...)
          {
           # Use createFromFile function for creating object
           object <- createFromFile(fileName=file,objectName="UnitCell")

           return(object)
          }
         )
#
setGeneric(
           name="reciprocalunitcell",
           def=function(ar,br,cr,alphar,betar,gammar,file,...){standardGeneric("reciprocalunitcell")}
          )
# All 6 parameters given (no default assignment)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","numeric","numeric","numeric","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(alphar)
           bb <- angle(betar)
           cc <- angle(gammar)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Three cell sides and one angle (alphar)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","numeric","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(alphar)
           bb <- angle(90)
           cc <- angle(90)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Three cell sides and one angle (betar)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","missing","numeric","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(betar)
           cc <- angle(90)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Three cell sides and one angle (gammar)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","missing","missing","numeric","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(gammar)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Three cell sides and no angles (orthorombic, tetragonal or cubic)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - ar)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","missing","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           br <- ar
           cr <- ar

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - br)
setMethod(
          f="reciprocalunitcell",
          signature=c("missing","numeric","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           ar <- br
           cr <- br

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - cr)
setMethod(
          f="reciprocalunitcell",
          signature=c("missing","missing","numeric","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           ar <- cr
           br <- cr

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# No cell sides and no angles (cubic with default side 1)
setMethod(
          f="reciprocalunitcell",
          signature=c("missing","missing","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           ar <- 1
           br <- 1
           cr <- 1

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Values taken from file (two methods, with and without "file=")
setMethod(
          f="reciprocalunitcell",
          signature=c("ANY","missing","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # In this case a is really containing a file name
           file <- ar

           # Use createFromFile function for creating object
           object <- createFromFile(fileName=file,objectName="ReciprocalUnitCell")

           return(object)
          }
         )
setMethod(
          f="reciprocalunitcell",
          signature=c("missing","missing","missing","missing","missing","missing","character"),
          function(ar,br,cr,alphar,betar,gammar,file,...)
          {
           # Use createFromFile function for creating object
           object <- createFromFile(fileName=file,objectName="ReciprocalUnitCell")

           return(object)
          }
         )
