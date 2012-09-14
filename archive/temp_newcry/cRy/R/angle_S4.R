# R code to implement the crystallographic ideas related to unit cell angles.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Default methods
#
# Print
setMethod(
          f="print",
          signature="Angle",
          definition=function(x,...)
                     {
                      if (!x@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.2f degrees. ***\n",x@ang)
                      if (x@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.4f radians. ***\n",x@ang)
                      cat(stringa)
                     }
         )

## Methods sharing generic functions with other objects
#
# Check Angle object is valid.
setMethod(
          f="isCorrect",
          signature="Angle",
          definition=function(object,name,message=FALSE)
                     {
                      # Stop if name is not "Angle"
                      if (name != "Angle")
                      {
                       msg <- paste("This is not an object of class ",name,". It is an object of class Angle.",sep="")
                       stop(msg)
                      }

                      # Extract Angle slots
                      slot_names <- getSlots(class(object))

                      # Check object has 2 slots (this should always be true)
                      #if (length(slot_names) != 2) stop("Object 'Angle' has a number of slots different from 2")

                      # Check nature of first slot
                      if (names(slot_names[1]) != "ang") stop("First slot in 'Angle' object does not appear to be named 'ang'")
                      if (!is.numeric(object@ang)) stop("First slot in 'Angle' object does not appear to be a number")

                      # Check nature of second slot
                      if (names(slot_names[2]) != "rad_flag") stop("Second slot in 'Angle' object does not appear to be named 'rad_flag'")
                      if (!is.logical(object@rad_flag)) stop("Second slot in 'Angle' object does not appear to be a logical variable")

                      # Check angle is a number between 0 and 180 (or 0 and pi)
                      value <- object@ang
                      if (!object@rad_flag & (value <= 0 | value >= 180)) stop("Numerical value for this angle is not in the 0 - 180 range")
                      if (object@rad_flag & (value <= 0 | value >= pi)) stop("Numerical value for this angle is not in the 0 - pi range")

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
#
setMethod(
          f="getFields",
          signature="Angle",
          definition=function(object)
                     {
                      lista <- list(ang=object@ang,rad_flag=object@rad_flag)

                      return(lista)
                     }
         )
#
setMethod(
          f="getParameters",
          signature="Angle",
          definition=function(object)
                     {
                      # For this object the parameters are simply the radians flag and the angle value
                      lista <- list(angle_value=object@ang)

                      return(lista)
                     }
         )
#
setMethod(
          f="getData",
          signature="Angle",
          definition=function(object)
                     {
                      # This class does not contain data and, thus, it returns NULL

                      return(NULL)
                     }
         )
#
# Another way to extract angle slots (similar to getFields)
setMethod(
          f="[",
          signature="Angle",
          definition=function(x,i,j,drop)
                     {
                      # One or two fields?
                      n <- length(i)

                      # Valid slots are: ang (numeric), rad_flag (logical)
                      fvect <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ang") fvect <- c(fvect,x@ang)
                       if (i[ii] == 1) fvect <- c(fvect,x@ang)
                       if (i[ii] == "rad_flag") fvect <- c(fvect,x@rad_flag)
                       if (i[ii] == 2) fvect <- c(fvect,x@rad_flag)
                       if (i[ii] != "ang" & i[ii] != "rad_flag" &
                           i[ii] != 1 & i[ii] != 2) fvect <- c(fvect,NA)
                      }

                      return(fvect)
                     }
         )
#
# Change ang or rad_flag value in object of class Angle
setReplaceMethod(
          f="[",
          signature="Angle",
          definition=function(x,i,j,value)
                     {
                      # One or two slots?
                      n <- length(i)

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: ang (numeric), rad_flag (logical)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ang") x@ang <- value
                       if (i[ii] == 1) x@ang <- value
                       if (i[ii] == "rad_flag") {
                                                 if (!x@rad_flag & value) x <- degToRad(x)
                                                 if (x@rad_flag & !value) x <- radToDeg(x)
                                                }
                       if (i[ii] == 2) {
                                        if (!x@rad_flag & value) x <- degToRad(x)
                                        if (x@rad_flag & !value) x <- radToDeg(x)
                                       }
                       if (i[ii] != "ang" & i[ii] != "rad_flag" &
                           i[ii] != 1 & i[ii] != 2) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      isCorrect(x,"Angle")

                      return(x)
                     }
         )

## Methods specific for UnitCell and ReciprocalUnitCell classes
#
# Change from degrees to radians
setMethod(
          f="degToRad",
          signature="Angle",
          definition=function(object){
                                      if (object@rad_flag) 
                                      {
                                       warning("Object of class Angle is already expressed in radians")
                                       return(object)
                                      }
                                      object@ang <- (object@ang*pi/180)%%(2*pi)
                                      object@rad_flag <- TRUE
                                      return(object)
                                     }
         )
#
# Change from radians to degrees
setMethod(
          f="radToDeg",
          signature="Angle",
          definition=function(object){
                                      if (!object@rad_flag) 
                                      {
                                       warning("Object of class Angle is already expressed in degrees") 
                                       return(object)
                                      }
                                      object@ang <- (object@ang*180/pi)%%360
                                      object@rad_flag <- FALSE
                                      return(object)
                                     }
         )

## Friendly constructor and functions
#
# Constructor
createAngle <- function(value=90.0,rad_flag=FALSE)
{
 # If |value| smaller than 6.28..., it might be in radians. Send a warning if rad_flag is FALSE
 if (abs(value) < 2*pi & !rad_flag) warning("Your input might be an angle in radians, but you haven't flagged it as such. Correct if necessary.")

 # If |value| greater or equal than 2*pi, but rad_flag is TRUE, maybe value is in degrees, but user has accidentally set flag to TRUE. Send a warning. 
 if (abs(value) >= 2*pi & rad_flag) warning("Your input might be in degrees, but you have flagged it as radians. Correct if necessary.")

 # Turn value in the 0 - 180 range
 if (rad_flag)
 {
  value <- value%%(2*pi)
  if (value > pi) value <- 2*pi-value
 }
 if (!rad_flag)
 {
  value <- value%%360
  if (value > 180) value <- 360-value
 }

 # Now create an instance of Angle class
 object <- new(Class="Angle",ang=value,rad_flag=rad_flag)

 # Check object is correct (execution stops if it is not)
 isCorrect(object,name="Angle")

 return(object)
}
