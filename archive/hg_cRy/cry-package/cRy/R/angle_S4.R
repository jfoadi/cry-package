# R code to implement the crystallographic ideas related to unit cell angles.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Classes
#
# Angle
setClass(
         Class="Angle",
         representation=representation(ang="numeric",rad_flag="logical")
        )

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

## Generic methods
#
# Check Angle object is valid.
setGeneric(
           name="checkAngle",
           def=function(object,message=FALSE){standardGeneric("checkAngle")} 
          )
#
setMethod(
          f="checkAngle",
          signature="Angle",
          definition=function(object,message=FALSE)
                     {
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
# Change Angle value and/or flag
setGeneric(
           name="changeAngle",
           def=function(object,value){standardGeneric("changeAngle")}
          )
#
setMethod(
          f="changeAngle",
          signature="Angle",
          definition=function(object,value)
                     {
                      # If value is a character ("degrees" or "radians") change rad_flag to FALSE or TRUE
                      if (is.character(value))
                      {
                       if (value == "degrees" & !object@rad_flag) return(object)
                       if (value == "degrees" & object@rad_flag) {object <- radToDeg(object); return(object)}
                       if (value == "radians" & object@rad_flag) return(object)
                       if (value == "radians" & !object@rad_flag) {object <- degToRad(object); return(object)}
                       if (value != "degrees" & value != "radians") stop("Invalid input value.")
                       return(object)
                      }
                      # If object angle is in degrees and outside 0 - 360 range, cast it back into that range
                      if (is.numeric(value) & !object@rad_flag)
                      {
                       value <- value%%360
                       object@ang <- value
                       return(object)
                      }

                      # If object angle is in radians and outside 0 - 2pi range, cast it back into that range
                      if (is.numeric(value) & object@rad_flag)
                      {
                       value <- value%%(2*pi)
                       object@ang <- value
                       return(object)
                      }
                     }
         )
#
# Degrees to radians conversion
setGeneric(
           name="degToRad",
           def=function(object){standardGeneric("degToRad")} 
          )
#
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
# Radians to degrees conversion
setGeneric(
           name="radToDeg",
           def=function(object){standardGeneric("radToDeg")} 
          )
#
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

## Friendly constructors
#
# Assign value to angle object while creating it
createAngle <- function(value,rad_flag=FALSE)
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
 checkAngle(object)

 return(object)
}
