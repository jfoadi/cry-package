# R code to implement the crystallographic ideas related to unit cell angles.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#


## Angle class
#
setClass(
         Class="Angle",
         representation=representation(ang="numeric",rad_flag="logical"),
         validity=function(object)
                  {
                   # Extract slots' content, even if they're empty
                   value <- object@ang
                   flag <- object@rad_flag

                   # Check angle is a number between 0 and 180 (or 0 and pi)
                   if (length(value) != 0 & length(flag) == 0)
                   {
                    if (object@ang < 0) return("Numerical value for this angle is not allowed because negative")
                   }
                   if (length(value) != 0 & length(flag) != 0)
                   {
                    if (!object@rad_flag & (value <= 0 | value >= 180)) return("Numerical value for this angle is not in the 0 - 180 range")
                    if (object@rad_flag & (value <= 0 | value >= pi)) return("Numerical value for this angle is not in the 0 - pi range")
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
          signature="Angle",
          definition=function(x,...)
                     {
                      if (length(x@rad_flag) != 0 & length(x@rad_flag) != 0)
                      {
                       if (!x@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.2f degrees. ***\n",x@ang)
                       if (x@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.4f radians. ***\n",x@ang)
                       cat(stringa)
                      }
                      if (length(x@rad_flag) == 0 & length(x@ang) != 0)
                      {
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is %6.2f, but it is unclear",x@ang)
                       stringa2 <- sprintf("whether this value is in degrees or radians. ***\n")
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(x@rad_flag) != 0 & length(x@ang) == 0)
                      {
                       if (x@rad_flag) fchar <- "radians"
                       if (!x@rad_flag) fchar <- "degrees"
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is not known, but")
                       stringa2 <- sprintf("it is measured in %s. ***\n",fchar) 
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(x@rad_flag) == 0 & length(x@ang) == 0)
                      {
                       stringa <- sprintf("*** This is an object of class Angle. Its slots are both empty. ***\n")
                       cat(stringa)
                      } 
                     }
         )
#
# Show
setMethod(
          f="show",
          signature="Angle",
          definition=function(object)
                     {
                      if (length(object@rad_flag) != 0 & length(object@rad_flag) != 0)
                      {
                       if (!object@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.2f degrees. ***\n",object@ang)
                       if (object@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.4f radians. ***\n",object@ang)
                       cat(stringa)
                      }
                      if (length(object@rad_flag) == 0 & length(object@ang) != 0)
                      {
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is %6.2f, but it is unclear",object@ang)
                       stringa2 <- sprintf("whether this value is in degrees or radians. ***\n")
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(object@rad_flag) != 0 & length(object@ang) == 0)
                      {
                       if (object@rad_flag) fchar <- "radians"
                       if (!object@rad_flag) fchar <- "degrees"
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is not known, but")
                       stringa2 <- sprintf("it is measured in %s. ***\n",fchar) 
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(object@rad_flag) == 0 & length(object@ang) == 0)
                      {
                       stringa <- sprintf("*** This is an object of class Angle. Its slots are both empty. ***\n")
                       cat(stringa)
                      } 
                     }
         )

## Non default methods corresponding to generic functions common to many cRy methods
#
# Extract Angle slots and return them in a list
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
# Extract Angle parameters
setMethod(
          f="getParameters",
          signature="Angle",
          definition=function(object)
                     {
                      # For this object the parameters are simply the radians flag and the angle value
                      lista <- list(ang=object@ang,rad_flag=object@rad_flag)

                      return(lista)
                     }
         )
#
# Extract Angle data (none)
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
# Another way to extract Angle slots (equivalent to getFields)
setMethod(
          f="[",
          signature="Angle",
          definition=function(x,i,j,drop)
                     {
                      # One or two fields?
                      n <- length(i)

                      # If n < 1 or n > 2 stop
                      if (n < 1 | n > 2) stop("This object includes only two values in its slots")

                      # Valid slots are: ang (numeric), rad_flag (logical)
                      flist <- list()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ang") flist <- c(flist,list(ang=x@ang))
                       if (i[ii] == 1) flist <- c(flist,list(ang=x@ang))
                       if (i[ii] == "rad_flag") flist <- c(flist,list(rad_flag=x@rad_flag))
                       if (i[ii] == 2) flist <- c(flist,list(rad_flag=x@rad_flag))
                       if (i[ii] != "ang" & i[ii] != "rad_flag" &
                           i[ii] != 1 & i[ii] != 2) stop("Slot is not included in object range")
                      }

                      return(flist)
                     }
         )
#
# Change parameters values in object of class Angle
setReplaceMethod(
          f="[",
          signature="Angle",
          definition=function(x,i,j,value)
                     {
                      # One or two slots?
                      n <- length(i)

                      # Only one or two slots accepted to replace values
                      if (n != 1 & n != 2) stop("This object of class Angle has only two slots")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Only one value is required to change
                      if (n == 1)
                      {
                       if (i == "ang") x@ang <- value
                       if (i == 1) x@ang <- value
                       if (i == "rad_flag")
                       {
                        if (!is.logical(value)) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value
                       }
                       if (i == 2)
                       {
                        if (!is.logical(value)) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value
                       }
                       if (i != 1 & i != 2 & i != "ang" & i != "rad_flag") stop("Slot to be replaced is not included in object range")
                      }

                      # Two values required to change
                      if (n == 2)
                      {
                       # Check that vlues are included in a list
                       if (!is.list(value)) stop("When replacing more than one value for object of class Angle, use a list")

                       if (i[1] == "ang") x@ang <- value[[1]]
                       if (i[1] == 1) x@ang <- value[[1]]
                       if (i[2] == "rad_flag")
                       {
                        if (!is.logical(value[[2]])) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value[[2]]
                       } 
                       if (i[2] == 2)
                       {
                        if (!is.logical(value[[2]])) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value[[2]]
                       } 
                       if (i[1] != 1 & i[1] != 2 & i[1] != "ang" & i[1] != "rad_flag") stop("Slot to be replaced is not included in object range")
                       if (i[2] != 1 & i[2] != 2 & i[2] != "ang" & i[2] != "rad_flag") stop("Slot to be replaced is not included in object range")
                      }

                      # Check changed object is correct
                      validObject(x)

                      return(x)
                     }
         )

## Non default methods corresponding to generic functions common only to Angle methods
#
# Degrees to radians conversion
setGeneric(
           name="degToRad",
           def=function(object){standardGeneric("degToRad")} 
          )
setMethod(
          f="degToRad",
          signature="Angle",
          definition=function(object){
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
          definition=function(object){
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

## Friendly constructors and functions
#
setGeneric(
           name="angle",
           def=function(x,y,...){standardGeneric("angle")} 
          )
# Both ang and rad_flag (no default assignment)
setMethod(
          f="angle",
          signature=c("numeric","logical"),
          function(x,y,...)
          {
           # Assign values to old code's variables
           value <- x
           rad_flag <- y

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

           return(object)
          }
         )
# Only ang (default assignment - FALSE - for rad_flag)
setMethod(
          f="angle",
          signature=c("numeric","missing"),
          function(x,y,...)
          {
           # Assign values to old code's variables
           value <- x

           # Values between 0 and 180 degrees
           value <- value%%360
           if (value > 180) value <- 360-value

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=value,rad_flag=FALSE)

           return(object)
          }
         )
# Only rad_flag (default assignment - 90 - for ang)
setMethod(
          f="angle",
          signature=c("logical","missing"),
          function(x,y,...)
          {
           # Assign values to old code's variables
           rad_flag <- x

           # ang value is fixed at 90.0 degrees 
           value <- 90.0

           # If rad_flag is TRUE, turn value in radians
           if (rad_flag) value <- 0.5*pi

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=value,rad_flag=rad_flag)

           return(object)
          }
         )
# No values (default assignment - 90 - for ang and FALSE for rad_flag)
setMethod(
          f="angle",
          signature=c("missing","missing"),
          function(x,y,...)
          {
           # ang value is fixed at 90.0 degrees 
           # rad_flag is fixed at FALSE

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=90,rad_flag=FALSE)

           return(object)
          }
         )
