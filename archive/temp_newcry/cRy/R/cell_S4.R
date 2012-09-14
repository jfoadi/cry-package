# R code to implement the crystallographic ideas related to unit cell.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Default methods
#
# Print cell
setMethod(
          f="print",
          signature="UnitCell",
          definition=function(x,...)
                     {
                      cat("*** This is a UnitCell object. ***\n")
                      cat(" The cell parameters are:\n")
                      stringa <- sprintf("       a = %7.2f\n       b = %7.2f\n       c = %7.2f\n   alpha =  %6.2f\n    beta =  %6.2f\n   gamma =  %6.2f\n",
                                         x@a,x@b,x@c,x@alpha@ang,x@beta@ang,x@gamma@ang)
                      cat(stringa)
                     }
         )
#
# Print reciprocal cell
setMethod(
          f="print",
          signature="ReciprocalUnitCell",
          definition=function(x,...)
                     {
                      cat("*** This is a ReciprocalUnitCell object. ***\n")
                      cat(" The reciprocal cell parameters are:\n")
                      stringa <- sprintf("       ar =  %9.6f\n       br =  %9.6f\n       cr =  %9.6f\n   alphar = %6.2f\n    betar = %6.2f\n   gammar = %6.2f\n",
                                         x@ar,x@br,x@cr,x@alphar@ang,x@betar@ang,x@gammar@ang)
                      cat(stringa)
                     }
         )

## Methods sharing generic functions with other objects
#
# Check UnitCell object is correct
setMethod(
          f="isCorrect",
          signature="UnitCell",
          definition=function(object,name,message=FALSE)
                     {
                      # Stop if name is not "UnitCell"
                      if (name != "UnitCell")
                      {
                       msg <- paste("This is not an object of class ",name,". It is an object of class UnitCell.",sep="")
                       stop(msg)
                      }

                      # Extract UnitCell slots
                      slot_names <- getSlots(class(object))

                      # Check object has 6 slots (this should always be true)
                      if (length(slot_names) != 6) stop("Object UnitCell has a number of slots different from 6")

                      # Check nature of first slot
                      if (names(slot_names[1]) != "a") stop("First slot in 'UnitCell' object does not appear to be named 'a'")
                      if (!is.numeric(object@a)) stop("First slot in 'UnitCel' object does not appear to be a number")

                      # Check nature of second slot
                      if (names(slot_names[2]) != "b") stop("Second slot in 'UnitCell' object does not appear to be named 'b'")
                      if (!is.numeric(object@b)) stop("Second slot in 'UnitCel' object does not appear to be a number")

                      # Check nature of third slot
                      if (names(slot_names[3]) != "c") stop("Third slot in 'UnitCell' object does not appear to be named 'c'")
                      if (!is.numeric(object@c)) stop("Third slot in 'UnitCel' object does not appear to be a number")
                     
                      # Check nature of fourth slot
                      if (names(slot_names[4]) != "alpha") stop("Fourth slot in 'UnitCell' does not appear to be named 'alpha'")
                      isCorrect(object@alpha,name="Angle")
                     
                      # Check nature of fifth slot
                      if (names(slot_names[5]) != "beta") stop("Fifth slot in 'UnitCell' does not appear to be named 'beta'")
                      isCorrect(object@beta,name="Angle")
                     
                      # Check nature of sixth slot
                      if (names(slot_names[6]) != "gamma") stop("Sixth slot in 'UnitCell' does not appear to be named 'gamma'")
                      isCorrect(object@gamma,name="Angle")

                      # Angles are objects of class Angles, so they all have values between 0 and 180 degrees. But combinations
                      # of alpha beta and gamma need to obey certain rules
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

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
#
# Check ReciprocalUnitCell object is correct
setMethod(
          f="isCorrect",
          signature="ReciprocalUnitCell",
          definition=function(object,name,message=FALSE)
                     {
                      # Stop if name is not "ReciprocalUnitCell"
                      if (name != "ReciprocalUnitCell")
                      {
                       msg <- paste("This is not an object of class ",name,". It is an object of class ReciprocalUnitCell.",sep="")
                       stop(msg)
                      }

                      # Extract ReciprocalUnitCell slots
                      slot_names <- getSlots(class(object))

                      # Check object has 6 slots (this should always be true)
                      if (length(slot_names) != 6) stop("Object ReciprocalUnitCell has a number of slots different from 6")

                      # Check nature of first slot
                      if (names(slot_names[1]) != "ar") stop("First slot in 'ReciprocalUnitCell' object does not appear to be named 'ar'")
                      if (!is.numeric(object@ar)) stop("First slot in 'ReciprocalUnitCel' object does not appear to be a number")

                      # Check nature of second slot
                      if (names(slot_names[2]) != "br") stop("Second slot in 'ReciprocalUnitCell' object does not appear to be named 'br'")
                      if (!is.numeric(object@br)) stop("Second slot in 'ReciprocalUnitCel' object does not appear to be a number")

                      # Check nature of third slot
                      if (names(slot_names[3]) != "cr") stop("Third slot in 'ReciprocalUnitCell' object does not appear to be named 'cr'")
                      if (!is.numeric(object@cr)) stop("Third slot in 'ReciprocalUnitCel' object does not appear to be a number")
                     
                      # Check nature of fourth slot
                      if (names(slot_names[4]) != "alphar") stop("Fourth slot in 'ReciprocalUnitCell' does not appear to be named 'alphar'")
                      isCorrect(object@alphar,name="Angle")
                     
                      # Check nature of fifth slot
                      if (names(slot_names[5]) != "betar") stop("Fifth slot in 'ReciprocalUnitCell' does not appear to be named 'betar'")
                      isCorrect(object@betar,name="Angle")
                     
                      # Check nature of sixth slot
                      if (names(slot_names[6]) != "gammar") stop("Sixth slot in 'ReciprocalUnitCell' does not appear to be named 'gammar'")
                      isCorrect(object@gammar,name="Angle")

                      # Angles are objects of class Angles, so they all have values between 0 and 180 degrees. But combinations
                      # of alpha beta and gamma need to obey certain rules
                      aa <- object@alphar@ang
                      bb <- object@betar@ang
                      cc <- object@gammar@ang
                      ss <- aa+bb+cc
                      if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                      ss <- aa+bb-cc
                      if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                      ss <- aa-bb+cc
                      if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                      ss <- -aa+bb+cc
                      if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist")

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
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
# Extract cell parameters
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
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) flist <- c(flist,NA)
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

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: a (numeric), b (numeric), c (numeric), alpha (Angle), beta (Angle) gamma (Angle)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") x@a <- value[ii]
                       if (i[ii] == 1) x@a <- value[ii]
                       if (i[ii] == "b") x@b <- value[ii]
                       if (i[ii] == 2) x@b <- value[ii]
                       if (i[ii] == "c") x@c <- value[ii]
                       if (i[ii] == 3) x@c <- value[ii]
                       if (i[ii] == "alpha") x@alpha@ang <- value[ii]
                       if (i[ii] == 4) x@alpha@ang <- value[ii]
                       if (i[ii] == "beta") x@beta@ang <- value[ii]
                       if (i[ii] == 5) x@beta@ang <- value[ii]
                       if (i[ii] == "gamma") x@gamma@ang <- value[ii]
                       if (i[ii] == 6) x@gamma@ang <- value[ii]
                       if (i[ii] != "a" & i[ii] != "b" & i[ii] != "c" &
                           i[ii] != "alpha" & i[ii] != "beta" & i[ii] != "gamma" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      isCorrect(x,"UnitCell")

                      return(x)
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
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) flist <- c(flist,NA)
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

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: ar (numeric), br (numeric), cr (numeric), alphar (Angle), betar (Angle) gammar (Angle)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ar") x@ar <- value[ii]
                       if (i[ii] == 1) x@ar <- value[ii]
                       if (i[ii] == "br") x@br <- value[ii]
                       if (i[ii] == 2) x@br <- value[ii]
                       if (i[ii] == "cr") x@cr <- value[ii]
                       if (i[ii] == 3) x@cr <- value[ii]
                       if (i[ii] == "alphar") x@alphar@ang <- value[ii]
                       if (i[ii] == 4) x@alphar@ang <- value[ii]
                       if (i[ii] == "betar") x@betar@ang <- value[ii]
                       if (i[ii] == 5) x@betar@ang <- value[ii]
                       if (i[ii] == "gammar") x@gammar@ang <- value[ii]
                       if (i[ii] == 6) x@gammar@ang <- value[ii]
                       if (i[ii] != "ar" & i[ii] != "br" & i[ii] != "cr" &
                           i[ii] != "alphar" & i[ii] != "betar" & i[ii] != "gammar" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      isCorrect(x,"ReciprocalUnitCell")

                      return(x)
                     }
         )

## Methods specific for UnitCell and ReciprocalUnitCell classes
#
# Compute reciprocal cell starting from direct cell
setMethod(
          f="computeReciprocalUnitCell",
          signature="UnitCell",
          definition=function(object)
                     {
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
setMethod(
          f="computeUnitCell",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # If object is already a direct cell, then return with no further action
                      if (class(object) == "UnitCell") stop("Object is already of class UnitCell")
                       
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
          signature="UnitCell",
          definition=function(object)
                     {
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

## Friendly constructors
#
# Create a unit cell with default values or values specified by the user
#createUnitCell <- function(a=1,b=1,c=1,alpha=90,beta=90,gamma=90)
UnitCell <- function(a=1,b=1,c=1,alpha=90,beta=90,gamma=90)
{
 # Turn angles into Angle objects
 aa <- createAngle(alpha)
 bb <- createAngle(beta)
 cc <- createAngle(gamma)

 # Create UnitCell object
 object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

 # Check object is correct (execution stops if not)
 isCorrect(object,name="UnitCell")

 return(object)
}
#
# Create a reciprocal unit cell with default values or values specified by the user
createReciprocalUnitCell <- function(ar=1,br=1,cr=1,alphar=90,betar=90,gammar=90)
{
 # Turn angles into Angle objects
 aa <- createAngle(alphar)
 bb <- createAngle(betar)
 cc <- createAngle(gammar)

 # Create ReciprocalUnitCell object
 object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

 # Check object is correct (execution stops if not)
 isCorrect(object,name="ReciprocalUnitCell")

 return(object)
}
