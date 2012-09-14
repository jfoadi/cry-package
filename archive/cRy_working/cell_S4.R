# R code to implement the crystallographic ideas related to unit cell.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Classes
#
# UnitCell
setClass(
         Class="UnitCell",
         representation=representation(a="numeric",b="numeric",c="numeric",alpha="Angle",beta="Angle",gamma="Angle")
        )
#
# ReciprocalUnitCell
setClass(
         Class="ReciprocalUnitCell",
         representation=representation(ar="numeric",br="numeric",cr="numeric",alphar="Angle",betar="Angle",gammar="Angle")
        )

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

## Generic methods
#
# Check UnitCell object is correct
setGeneric(
           name="checkUnitCell",
           def=function(object,message=FALSE){standardGeneric("checkUnitCell")}
          )
#
setMethod(
          f="checkUnitCell",
          signature="UnitCell",
          definition=function(object,message=FALSE)
                     {
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
                      checkAngle(object@alpha)
                     
                      # Check nature of fifth slot
                      if (names(slot_names[5]) != "beta") stop("Fifth slot in 'UnitCell' does not appear to be named 'beta'")
                      checkAngle(object@beta)
                     
                      # Check nature of sixth slot
                      if (names(slot_names[6]) != "gamma") stop("Sixth slot in 'UnitCell' does not appear to be named 'gamma'")
                      checkAngle(object@gamma)

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
setGeneric(
           name="checkReciprocalUnitCell",
           def=function(object,message=FALSE){standardGeneric("checkReciprocalUnitCell")}
          )
#
setMethod(
          f="checkReciprocalUnitCell",
          signature="ReciprocalUnitCell",
          definition=function(object,message=FALSE)
                     {
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
                      checkAngle(object@alphar)
                     
                      # Check nature of fifth slot
                      if (names(slot_names[5]) != "betar") stop("Fifth slot in 'ReciprocalUnitCell' does not appear to be named 'betar'")
                      checkAngle(object@betar)
                     
                      # Check nature of sixth slot
                      if (names(slot_names[6]) != "gammar") stop("Sixth slot in 'ReciprocalUnitCell' does not appear to be named 'gammar'")
                      checkAngle(object@gammar)

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
# Change/Assign cell parameters
setGeneric(
           name="changeCellParameters",
           def=function(object,a=NA,b=NA,c=NA,alpha=NA,beta=NA,gamma=NA){standardGeneric("changeCellParameters")}
          )
#
setMethod(
          f="changeCellParameters",
          signature="UnitCell",
          definition=function(object,a=NA,b=NA,c=NA,alpha=NA,beta=NA,gamma=NA)
                     {
                      if (!is.na(a)) object@a <- a
                      if (!is.na(b)) object@b <- b
                      if (!is.na(c)) object@c <- c
                      if (!is.na(alpha)) object@alpha@ang <- alpha
                      if (!is.na(beta)) object@beta@ang <- beta
                      if (!is.na(gamma)) object@gamma@ang <- gamma

                      return(object)
                     } 
         )
#
# Change/Assign reciprocal cell parameters
#
setMethod(
          f="changeCellParameters",
          signature="ReciprocalUnitCell",
          definition=function(object,a=NA,b=NA,c=NA,alpha=NA,beta=NA,gamma=NA)
                     {
                      if (!is.na(a)) object@ar <- a
                      if (!is.na(b)) object@br <- b
                      if (!is.na(c)) object@cr <- c
                      if (!is.na(alpha)) object@alphar@ang <- alpha 
                      if (!is.na(beta)) object@betar@ang <- beta
                      if (!is.na(gamma)) object@gammar@ang <- gamma

                      return(object)
                     } 
         )
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
                      # If object is already a reciprocal cell, then return with no further action
                      if (class(object) == "ReciprocalUnitCell") stop("Object is already of class ReciprocalUnitCell")
                       
                      # Now compute the 6-parameters list
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
setGeneric(
           name="computeCellVolume",
           def=function(object){standardGeneric("computeCellVolume")}
          )
#
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
#
# Extract cell parameters into a length-6-vector (valid for both direct and reciprocal cells)
setGeneric(
           name="extractCellParameters",
           def=function(object){standardGeneric("extractCellParameters")}
          )
#
setMethod(
          f="extractCellParameters",
          signature="UnitCell",
          definition=function(object)
                     {
                      # Check object is of class UnitCell
                      checkUnitCell(object)

                      # Extract parameters
                      a <- object@a
                      b <- object@b
                      c <- object@c
                      alpha <- object@alpha@ang
                      beta <- object@beta@ang
                      gamma <- object@gamma@ang
                      cpars <- c(a,b,c,alpha,beta,gamma)

                      return(cpars)
                     } 
         )
#
setMethod(
          f="extractCellParameters",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # Check object is of class ReciprocalUnitCell
                      checkReciprocalUnitCell(object)

                      # Extract parameters
                      ar <- object@ar
                      br <- object@br
                      cr <- object@cr
                      alphar <- object@alphar@ang
                      betar <- object@betar@ang
                      gammar <- object@gammar@ang
                      cpars <- c(ar,br,cr,alphar,betar,gammar)

                      return(cpars)
                     } 
         )

## Friendly constructors
#
# Create a unit cell with default values or values specified by the user
createUnitCell <- function(a=1,b=1,c=1,alpha=90,beta=90,gamma=90)
{
 # Turn angles into Angle objects
 aa <- createAngle(alpha)
 bb <- createAngle(beta)
 cc <- createAngle(gamma)

 # Create UnitCell object
 object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

 # Check object is correct (execution stops if not)
 checkUnitCell(object)

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

 return(object)
}
#
# Create a UnitCell starting from a pdb file
#createUnitCellFromPDB <- function(file)
#{
 # Load pdb file into a named list
# lpdb <- .readPDB(file)
 
 # Cell parameters
# cpar <- lpdb$cryst1$cell_par

 # Create new UnitCell object
# object <- createUnitCell(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])

# return(object)
#}
#
# Create a UnitCell starting from an mtz file
createUnitCellFromMTZ <- function(file)
{
 # Load mtz file into a named list
 lmtz <- .readMTZ(file,messages=FALSE)
 
 # Extract cell parameters
 cpar <- lmtz$header$CELL

 # Create new UnitCell object
 object <- createUnitCell(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])

 return(object)
}
