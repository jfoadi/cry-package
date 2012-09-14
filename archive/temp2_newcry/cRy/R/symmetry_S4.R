# R code to implement the crystallographic ideas related to symmetry.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Symmetry class
#
# Symmetry
setClass(
         Class="Symmetry",
         representation=representation(sym_xHM="character"),
         validity=function(object)
         {
          # If object has been created empty, or with some empty slots, fine.
          if (length(object@sym_xHM) == 0)
          {
           # Green light: object appears to be correct
           return(TRUE)
          }
 
          # Check extended Herman-Maguin symbol is correct
          if (length(object@sym_xHM) != 0)
          {
           # To be used throughout this block
           tmp <- .translate_SG(object@sym_xHM)

           lista <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
           if (!lista$ans) stop(lista$msg)
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
          signature="Symmetry",
          definition=function(x,...)
                     {
                      # Initial message
                      cat("*** This is a Symmetry object ***\n")

                      # Empty object
                      if (length(x@sym_xHM) == 0)
                      {
                       stringa <- sprintf(" The slot for its extended Herman-Maguin symbol is empty.\n")
                       cat(stringa)
                      }

                      if (length(x@sym_xHM) != 0)
                      {
                       # Find space group number
                       lista <- .translate_SG(value=x@sym_xHM,SG_in="xHM",SG_out="number")
                       sgn <- lista$msg

                       # Find crystal system
                       cr_sys <- .crystal_system(sgn)

                       # Find symmetry and centring operations
                       lista <- .syminfo_to_op_xyz_list(x@sym_xHM)
                       xyz <- lista[[1]]
                       cenop <- lista[[2]]

                       # Print
                       stringa <- sprintf("    The space group extended Hermann-Maguin symbol is %s\n",x@sym_xHM)
                       cat(stringa)
                       stringa <- sprintf("    The space group number from the International Tables is %d\n",sgn)
                       cat(stringa)
                       stringa <- sprintf("    This space group belongs to the %s crystal system\n",cr_sys)
                       cat(stringa)
                       stringa <- sprintf("    This space group includes the following symmetry operations:\n")
                       cat(stringa)
                       for (i in 1:length(xyz))
                       {
                        stringa <- sprintf("    %2d) %s\n",i,xyz[i])
                        cat(stringa)
                       }
                       if (length(cenop) == 1)
                       {
                        stringa <- sprintf("    There are no centring operators associated with this group\n")
                        cat(stringa)
                       }
                       if (length(cenop) != 1)
                       {
                        stringa <- sprintf("    This space group is associated with a centred cell. The centring operators are:\n")
                        cat(stringa)
                        for (i in 2:length(cenop))
                        {
                         stringa <- sprintf("    %2d) %s\n",(i-1),cenop[i])
                         cat(stringa)
                        }
                       }
                      }
                     }
         )
#
# Show
setMethod(
          f="show",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Initial message
                      cat("*** This is a Symmetry object ***\n")

                      # Empty object
                      if (length(object@sym_xHM) == 0)
                      {
                       stringa <- sprintf(" The slot for its extended Herman-Maguin symbol is empty.\n")
                       cat(stringa)
                      }

                      if (length(object@sym_xHM) != 0)
                      {
                       # Find space group number
                       lista <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
                       sgn <- lista$msg

                       # Find crystal system
                       cr_sys <- .crystal_system(sgn)

                       # Find symmetry and centring operations
                       lista <- .syminfo_to_op_xyz_list(object@sym_xHM)
                       xyz <- lista[[1]]
                       cenop <- lista[[2]]

                       # Print
                       stringa <- sprintf("    The space group extended Hermann-Maguin symbol is %s\n",object@sym_xHM)
                       cat(stringa)
                       stringa <- sprintf("    The space group number from the International Tables is %d\n",sgn)
                       cat(stringa)
                       stringa <- sprintf("    This space group belongs to the %s crystal system\n",cr_sys)
                       cat(stringa)
                       stringa <- sprintf("    This space group includes the following symmetry operations:\n")
                       cat(stringa)
                       for (i in 1:length(xyz))
                       {
                        stringa <- sprintf("    %2d) %s\n",i,xyz[i])
                        cat(stringa)
                       }
                       if (length(cenop) == 1)
                       {
                        stringa <- sprintf("    There are no centring operators associated with this group\n")
                        cat(stringa)
                       }
                       if (length(cenop) != 1)
                       {
                        stringa <- sprintf("    This space group is associated with a centred cell. The centring operators are:\n")
                        cat(stringa)
                        for (i in 2:length(cenop))
                        {
                         stringa <- sprintf("    %2d) %s\n",(i-1),cenop[i])
                         cat(stringa)
                        }
                       }
                      }
                     }
         )

## Methods sharing generic functions with other objects
# 
# Extract Symmetry slots and return them in a list
setMethod(
          f="getFields",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Extract slot
                      sym_xHM <- object@sym_xHM
                      lista <- list(sym_xHM=sym_xHM)

                      return(lista)
                     }
         )
# 
# Extract parameters from Lattice object
setMethod(
          f="getParameters",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Build list to be returned
                      lista <- list(sym_xHM=object@sym_xHM)

                      return(lista)
                     }
         )
#
# No data returned for Symmetry object
setMethod(
          f="getData",
          signature="Symmetry",
          definition=function(object)
                     {
                      # No data are contained in an object of class Symmetry

                      return(NULL)
                     }
         )
#
# Another way to extract Symmetry slots (similar to getFields)
setMethod(
          f="[",
          signature="Symmetry",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Only 1 parameter for object of class Symmetry
                      if (n != 1) stop("Object of class Symmetry has only 1 parameter")

                      # Create list with required content in it
                      if (i[1] == "sym_xHM") flist <- list(sym_xHM=x@sym_xHM) 
                      if (i[1] == 1) flist <- list(sym_xHM=x@sym_xHM) 
                      if (i[1] != "sym_xHM" & i[1] != 1) stop("Parameter not included in list for this object")

                      return(flist)
                     }
         ) 
#
# Change parameter value in object of class Symmetry
setReplaceMethod(
                 f="[",
                 signature="Symmetry",
                 definition=function(x,i,j,value)
                            {
                             # Find out length of i
                             n <- length(i)

                             # Only 1 parameter for object of class Symmetry
                             if (n != 1) stop("Object of class Symmetry has only 1 parameter")

                             # Parameter type is: sym_xHM (character)
                             if (i[1] == "sym_xHM")
                             {
                              if (is.list(value))
                              {
                               value[[1]] <- .findHM(value[[1]])
                               x@sym_xHM <- value[[1]]
                              }
                              if (!is.list(value)) 
                              {
                               value <- .findHM(value)
                               x@sym_xHM <- value
                              }
                             }
                             if (i[1] == 1)
                             {
                              if (is.list(value))
                              {
                               value[[1]] <- .findHM(value[[1]])
                               x@sym_xHM <- value[[1]]
                              }
                              if (!is.list(value)) 
                              {
                               value <- .findHM(value)
                               x@sym_xHM <- value
                              }
                             }
                             if (i[1] != "sym_xHM" & i[1] != 1) stop("One or more requested slots do not exist")

                             if (validObject(x)) return(x)
                            }
                )

## Methods specific for Symmetry class
#
# Return centering operators (shared with Lattice)
setMethod(
          f="getCentringOps",
          signature="Symmetry",
          definition=function(object)
                     {
                      ltmp <- .syminfo_to_matrix_list(object@sym_xHM)
                      lista <- ltmp$C

                      return(lista)
                     }
         )
#
# Return number of settings for a same space group number
setGeneric(
           name="getNumberSettings",
           def=function(object){standardGeneric("getNumberSettings")}
          )
setMethod(
          f="getNumberSettings",
          signature="Symmetry",
          definition=function(object)
                     {
                      # First get space group number
                      lista <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
                      sgn <- lista$msg

                      # Now cycle until no more settings for that specific space group number are found
                      nsett <- 0
                      while (lista$ans)
                      {
                       nsett <- nsett+1 
                       lista <- .translate_SG(value=sgn,SG_in="number",SG_out="xHM",setting=nsett)
                      }
                      nsett <- nsett-1

                      return(nsett)
                     }
         )
#
# Return space group number and setting number related to an extended H-M symbol
setGeneric(
           name="getSymmetryNumber",
           def=function(object){standardGeneric("getSymmetryNumber")}
          )
setMethod(
          f="getSymmetryNumber",
          signature="Symmetry",
          definition=function(object)
                     {
                      # First space group number
                      lista <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
                      sgn <- lista$msg

                      # Now count how many settings
                      nsett <- 0
                      while (lista$ans)
                      {
                       nsett <- nsett+1 
                       lista <- .translate_SG(value=sgn,SG_in="number",SG_out="xHM",setting=nsett)
                      }
                      nsett <- nsett-1

                      # To finish, check which xHM symbol coincides with the one related to number and setting
                      setting <- NULL
                      for (i in 1:nsett)
                      {
                       lista <- .translate_SG(value=sgn,SG_in="number",SG_out="xHM",setting=i)
                       if (lista$msg == object@sym_xHM) setting <- i
                      }

                      return(c(sgn,setting))
                     }
         )
#
# Return space group extended H-M symbol
setGeneric(
           name="getSymmetryHM",
           def=function(object){standardGeneric("getSymmetryHM")}
          )
setMethod(
          f="getSymmetryHM",
          signature="Symmetry",
          definition=function(object)
                     {
                      return(object@sym_xHM)
                     }
         )
#
# Compute symmetry operators
setGeneric(
           name="computeSymmetryOps",
           def=function(object){standardGeneric("computeSymmetryOps")}
          )
setMethod(
          f="computeSymmetryOps",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Collect operators in matrix form in a list with 3 fields:
                      #  1) $PG a list whose members are 3X3 matrices        (rotation/inversion part of symmetry)
                      #  2) $T  a list whose members are vectors of length 3 (translational part of symmetry)
                      #  3) $C  a list whose members are vectors of length 3 (cell centring)
                      ltmp <- .syminfo_to_matrix_list(object@sym_xHM)

                      return(ltmp)
                     }
         )
#
# Check symmetry is compatible with unit cell parameters
setGeneric(
           name="checkSymmetryWithCell",
           def=function(object,cell){standardGeneric("checkSymmetryWithCell")}
          )
setMethod(
          f="checkSymmetryWithCell",
          signature=c("Symmetry","UnitCell"),
          definition=function(object,cell)
                     {
                     }
         )


## Friendly constructors and functions
#
setGeneric(
           name="symmetry",
           def=function(sym_xHM,sym_number,setting,file,...){standardGeneric("symmetry")}
          )
#
# Starting from character. Could be either extended Herman-Maguin symbol, or a file. If not an xHM, it tries to interpret it
setMethod(
          f="symmetry",
          signature=c("character","missing","missing","missing"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           # First try for standard symbols, extend them to xHM and create object
           sym_xHM <- .findHM(sym_xHM)

           # Check if sym_xHM is valid
           lista <- .translate_SG(value=sym_xHM,SG_in="xHM",SG_out="number")
           if (lista$ans)
           {
            object <- new(Class="Symmetry",sym_xHM=sym_xHM)
            return(object)
           }

           # If not a valid xHM symbol, it could be a file name
           if (!lista$ans)
           {
            ans <- file.exists(sym_xHM)
            if (ans)
            {
             object <- createFromFile(sym_xHM,objectName="Symmetry")
             return(object)
            }
            if (!ans)
            {
             cat("*** Input is neither a valid space group symbol, nor an existing file name ***\n")
             return(NULL)
            }
           }
          }
         )
#
# Starting from space group number and setting
setMethod(
          f="symmetry",
          signature=c("numeric","numeric","missing","missing"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           # Find extended HM symbol
           lista <- .translate_SG(sym_xHM,setting=sym_number)
           if (!lista$ans) stop(lista$msg)
           xHM <- lista$msg

           # Create object of class Symmetry
           object <- new(Class="Symmetry",sym_xHM=xHM)

           return(object)
          }
         )
#
# Starting from space group number (setting is fixed to 1)
setMethod(
          f="symmetry",
          signature=c("numeric","missing","missing","missing"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           # Find extended HM symbol
           lista <- .translate_SG(sym_xHM)
           if (!lista$ans) stop(lista$msg)
           xHM <- lista$msg

           # Create object of class Symmetry
           object <- new(Class="Symmetry",sym_xHM=xHM)

           return(object)
          }
         )
#
# Starting from space group number and setting, defined as sym_number and setting
setMethod(
          f="symmetry",
          signature=c("missing","numeric","numeric","missing"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           # Find extended HM symbol
           lista <- .translate_SG(sym_number,setting=setting)
           if (!lista$ans) stop(lista$msg)
           xHM <- lista$msg

           # Create object of class Symmetry
           object <- new(Class="Symmetry",sym_xHM=xHM)

           return(object)
          }
         )
#
# Starting from space group number, sym_number (setting is fixed to 1)
setMethod(
          f="symmetry",
          signature=c("missing","numeric","missing","missing"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           # Find extended HM symbol
           lista <- .translate_SG(sym_number)
           if (!lista$ans) stop(lista$msg)
           xHM <- lista$msg

           # Create object of class Symmetry
           object <- new(Class="Symmetry",sym_xHM=xHM)

           return(object)
          }
         )
#
# Starting from character, file.
setMethod(
          f="symmetry",
          signature=c("missing","missing","missing","character"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           # Check the file does exist
           ans <- file.exists(file)
           if (ans)
           {
            object <- createFromFile(file,objectName="Symmetry")
            return(object)
           }
           if (!ans)
           {
            stop("*** Input is not an existing file name ***\n")
           }
          }
         )
#
# Starting from nothing. Default symmetry id cubic P 2 3
setMethod(
          f="symmetry",
          signature=c("missing","missing","missing","missing"),
          function(sym_xHM,sym_number,setting,file,...)
          {
           object <- new(Class="Symmetry",sym_xHM="P 2 3")

           return(object)
          }
         )
