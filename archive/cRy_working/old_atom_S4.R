# R code to implement the crystallographic ideas related to atoms.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#
# !!! IMPORTANT !!! The symbol we use for sodium is "Na", rather than "NA", as the last one clashes
# with the R missing data type.
#

## Classes
#
# Ccoord
setClass(
         Class="Ccoord",
         representation=representation(xyz="matrix")
        )
#
# Fcoord
setClass(
         Class="Fcoord",
         representation=representation(xyzf="matrix",cell="UnitCell")
        )
#
# AtomSymbol
setClass(
         Class="AtomSymbol",
         representation=representation(symbol="character")
        )
#
# Atom
setClass(
         Class="Atom",
         representation=representation(symbol="AtomSymbol",xyz="Ccoord",bfac="list")
        )
#
# AtomInCell
setClass(
         Class="AtomInCell",
         representation=representation(symbol="AtomSymbol",xyzf="Fcoord",bfac="matrix")
        )

## Default methods
#
# Print Atom
setMethod(
          f="print",
          signature="Atom",
          definition=function(x,...)
                     {
                      # Check Atom object is correct
                      checkAtom(x)

                      cat("*** This is an object of class Atom ***\n")
                      idx <- which(x@symbol@symbol == .ATOMS_data.frame)
                      element <- .ATOMS_data.frame[idx,2]
                      stringa <- sprintf("    The chemical element is %s\n",element)
                      cat(stringa)
                      stringa <- sprintf("    Its location is given by coordinates (%6.2f,%6.2f,%6.2f)\n",x@xyz@xyz[1],x@xyz@xyz[2],x@xyz@xyz[3])
                      cat(stringa)
                      if (abs(x@bfac[1,2]-0) < 0.000001 & abs(x@bfac[1,3]-0) < 0.000001 & abs(x@bfac[2,3]-0) < 0.000001 &
                          abs(x@bfac[1,1]-x@bfac[2,2]) < 0.000001 & abs(x@bfac[2,2]-x@bfac[3,3]) < 0.000001)
                      {
                       stringa <- sprintf("    Its thermic factor is isotropic, of value %6.2f\n",x@bfac[1,1]) 
                       cat(stringa)
                      }
                      if (abs(x@bfac[1,1]-x@bfac[2,2]) > 0.000001 | abs(x@bfac[2,2]-x@bfac[3,3]) > 0.000001 | abs(x@bfac[1,2]-0) > 0.000001 |
                          abs(x@bfac[1,3]-0) > 0.000001 | abs(x@bfac[2,3]-0) > 0.000001)
                      {
                       stringa <- sprintf("    Its thermic factor is anisotropic, of matrix ||%6.2f %6.2f %6.2f||\n",
                                                                                                        x@bfac[1,1],x@bfac[1,2],x@bfac[1,3])
                       cat(stringa)
                       stringa <- sprintf("                                                 ||%6.2f %6.2f %6.2f||\n",
                                                                                                        x@bfac[2,1],x@bfac[2,2],x@bfac[2,3])
                       cat(stringa)
                       stringa <- sprintf("                                                 ||%6.2f %6.2f %6.2f||\n",
                                                                                                        x@bfac[3,1],x@bfac[3,2],x@bfac[3,3])
                       cat(stringa)
                      }
                     }
         )
#
# Print AtomInCell
setMethod(
          f="print",
          signature="AtomInCell",
          definition=function(x,...)
                     {
                      # Check AtomInCell object is correct
                      checkAtomInCell(x)

                      cat("*** This is an object of class AtomInCell ***\n")
                      idx <- which(x@symbol@symbol == .ATOMS_data.frame)
                      element <- .ATOMS_data.frame[idx,2]
                      stringa <- sprintf("    The chemical element is %s\n",element)
                      cat(stringa)
                      cat("    The unit cell in which it is located has the following parameters:\n")
                      stringa <- sprintf("        a = %7.2f\n",x@xyzf@cell@a)
                      cat(stringa)
                      stringa <- sprintf("        b = %7.2f\n",x@xyzf@cell@b)
                      cat(stringa)
                      stringa <- sprintf("        c = %7.2f\n",x@xyzf@cell@c)
                      cat(stringa)
                      stringa <- sprintf("    alpha =  %6.2f\n",x@xyzf@cell@alpha@ang)
                      cat(stringa)
                      stringa <- sprintf("     beta =  %6.2f\n",x@xyzf@cell@beta@ang)
                      cat(stringa)
                      stringa <- sprintf("    gamma =  %6.2f\n",x@xyzf@cell@gamma@ang)
                      cat(stringa)
                      stringa <- sprintf("    Its location is given by fractional coordinates (%7.5f,%7.5f,%7.5f)\n",
                                                                                              x@xyzf@xyzf[1],x@xyzf@xyzf[2],x@xyzf@xyzf[3])
                      cat(stringa)
                      if (abs(x@bfac[1,2]-0) < 0.000001 & abs(x@bfac[1,3]-0) < 0.000001 & abs(x@bfac[2,3]-0) < 0.000001 &
                          abs(x@bfac[1,1]-x@bfac[2,2]) < 0.000001 & abs(x@bfac[2,2]-x@bfac[3,3]) < 0.000001)
                      {
                       stringa <- sprintf("    Its thermic factor is isotropic, of value %6.2f\n",x@bfac[1,1]) 
                       cat(stringa)
                      }
                      if (abs(x@bfac[1,1]-x@bfac[2,2]) > 0.000001 | abs(x@bfac[2,2]-x@bfac[3,3]) > 0.000001 | abs(x@bfac[1,2]-0) > 0.000001 |
                          abs(x@bfac[1,3]-0) > 0.000001 | abs(x@bfac[2,3]-0) > 0.000001)
                      {
                       stringa <- sprintf("    Its thermic factor is anisotropic, of matrix ||%6.2f %6.2f %6.2f||\n",
                                                                                                        x@bfac[1,1],x@bfac[1,2],x@bfac[1,3])
                       cat(stringa)
                       stringa <- sprintf("                                                 ||%6.2f %6.2f %6.2f||\n",
                                                                                                        x@bfac[2,1],x@bfac[2,2],x@bfac[2,3])
                       cat(stringa)
                       stringa <- sprintf("                                                 ||%6.2f %6.2f %6.2f||\n",
                                                                                                        x@bfac[3,1],x@bfac[3,2],x@bfac[3,3])
                       cat(stringa)
                      }
                     }
         )

## Generic methods
#
# Check Ccoord object is correct
setGeneric(
           name="checkCcoord",
           def=function(object,message=FALSE){standardGeneric("checkCcoord")}
          )
#
setMethod(
          f="checkCcoord",
          signature="Ccoord",
          definition=function(object,message=FALSE)
                     {
                      # Stop if input vector has length different from 3
                      if (length(object@xyz) != 3) stop("This object is not of class Ccoord as it does not appear to have 3 components")

                      # All checks passed: green light
                      if (message == TRUE) print("This object and its slots appear to be correct")
                     }
         )
#
# Check Fcoord object is correct
setGeneric(
           name="checkFcoord",
           def=function(object,message=FALSE){standardGeneric("checkFcoord")}
          )
#
setMethod(
          f="checkFcoord",
          signature="Fcoord",
          definition=function(object,message=FALSE)
                     {
                      # Stop if input slot xyzf has length different from 3
                      if (length(object@xyzf) != 3) stop("This object of class Fcoord has slot xyzf with a number of components different from 3")
                      
                      # Stop if UnitCell object is not correct 
                      checkUnitCell(object@cell)

                      # All checks passed: green light
                      if (message == TRUE) print("This object and its slots appear to be correct")
                     }
         )
#
# Check AtomSymbol object is correct
setGeneric(
           name="checkAtomSymbol",
           def=function(object,message=FALSE){standardGeneric("checkAtomSymbol")}
          )
#
setMethod(
          f="checkAtomSymbol",
          signature="AtomSymbol",
          definition=function(object,message=FALSE)
                     {
                      # Stop if input character has more than 2 letters
                      if (length(object@symbol) > 2) stop("The symbol slot of this object can only be a string with a maximum of 2 characters")
           
                      # Check symbol is of a recognised element
                      stmp <- object@symbol
                      if (nchar(strsplit(stmp," ")[[1]][1]) == 0) stmp <- substr(stmp,2,2)
                      idx <- which(stmp == .ATOMS_data.frame)
                      if (length(idx) == 0) stop("Atomic symbol not recognised in object of class AtomSymbol")

                      # All checks passed: green light
                      if (message == TRUE) print("This object and its slots appear to be correct")
                     }
         )
#
# Check Atom object is correct
setGeneric(
           name="checkAtom",
           def=function(object,message=FALSE){standardGeneric("checkAtom")}
          )
#
setMethod(
          f="checkAtom",
          signature="Atom",
          definition=function(object,message=FALSE)
                     {
                      # Check symbol slot is a correct AtomSymbol object
                      checkAtomSymbol(object@symbol)

                      # Check bfac is a numeric 3X3 matrix
                      if (!is.numeric(object@bfac)) stop("bfac slot of Atom object needs to be a 3X3 matrix of numbers")
                      if (is.null(dim(object@bfac))) stop("bfac slot of Atom object needs to be a 3X3 matrix of numbers")
                      if (dim(object@bfac)[1] != 3 | dim(object@bfac)[2] != 3) stop("bfac slot of Atom object needs to be a 3X3 matrix of numbers") 

                      # Check bfac is a symmetric matrix
                      if (object@bfac[1,2] != object@bfac[2,1] | object@bfac[1,3] != object@bfac[3,1] | object@bfac[2,3] != object@bfac[3,2])
                                       stop("bfac slot of Atom object needs to be a symmetric matrix")

                      # All checks passed: green light
                      if (message == TRUE) print("This object and its slots appear to be correct")
                     }
         )
#
# Check AtomInCell object is correct
setGeneric(
           name="checkAtomInCell",
           def=function(object,message=FALSE){standardGeneric("checkAtomInCell")}
          )
#
setMethod(
          f="checkAtomInCell",
          signature="AtomInCell",
          definition=function(object,message=FALSE)
                     {
                      # Check symbol slot is a correct AtomSymbol object
                      checkAtomSymbol(object@symbol)

                      # Check xyzf is a correct Fcoord object
                      checkFcoord(object@xyzf)

                      # Check bfac is a numeric 3X3 matrix
                      if (!is.numeric(object@bfac)) stop("bfac slot of AtomInCell object needs to be a 3X3 matrix of numbers")
                      if (is.null(dim(object@bfac))) stop("bfac slot of AtomInCell object needs to be a 3X3 matrix of numbers")
                      if (dim(object@bfac)[1] != 3 | dim(object@bfac)[2] != 3) 
                                                                         stop("bfac slot of AtomInCell object needs to be a 3X3 matrix of numbers") 

                      # Check bfac is a symmetric matrix
                      if (object@bfac[1,2] != object@bfac[2,1] | object@bfac[1,3] != object@bfac[3,1] | object@bfac[2,3] != object@bfac[3,2])
                                       stop("bfac slot of AtomInCell object needs to be a symmetric matrix")

                      # All checks passed: green light
                      if (message == TRUE) print("This object and its slots appear to be correct")
                     }
         )
#
# Transform Fcoord in Ccoord
setGeneric(
           name="fromFcoordToCcoord",
           def=function(object,ochoice=1){standardGeneric("fromFcoordToCcoord")}
          )
#
setMethod(
          f="fromFcoordToCcoord",
          signature="Fcoord",
          definition=function(object,ochoice=1)
                     {
                      # Check object is of class Fcoord
                      checkFcoord(object)

                      # Check ochoice is 1 or 2
                      if (ochoice != 1 & ochoice != 2) stop("Only values 1 and 2 accepted for ochoice")

                      # Use function .frac_to_orth
                      tmp <- .frac_to_orth(object@xyzf,object@cell@a,object@cell@b,object@cell@c,
                                           object@cell@alpha@ang,object@cell@beta@ang,object@cell@gamma@ang,ochoice=ochoice)
      
                      # Now create Ccoord object
                      new_object <- createCcoord(as.numeric(tmp))

                      return(new_object)
                     }
         )
#
# Transform Ccoord in Fcoord
setGeneric(
           name="fromCcoordToFcoord",
           def=function(object,cell,ochoice=1){standardGeneric("fromCcoordToFcoord")}
          )
#
setMethod(
          f="fromCcoordToFcoord",
          signature="Ccoord",
          definition=function(object,cell,ochoice=1)
                     {
                      # Check object is of class Ccoord
                      checkCcoord(object)

                      # Check cell is an object of class UnitCell
                      checkUnitCell(cell)

                      # Check ochoice is 1 or 2
                      if (ochoice != 1 & ochoice != 2) stop("Only values 1 and 2 accepted for ochoice")

                      # Use function .frac_to_orth
                      tmp <- .orth_to_frac(object@xyz,cell@a,cell@b,cell@c,cell@alpha@ang,cell@beta@ang,cell@gamma@ang,ochoice=ochoice)
      
                      # Now create Fcoord object
                      new_object <- createFcoord(as.numeric(tmp),cell)

                      return(new_object)
                     }
         )

## Friendly constructors
#
# Create a Ccoord object starting from a 3D vector 
createCcoord <- function(xyz=c(0,0,0))
{
 # Create object
 object <- new(Class="Ccoord",xyz=xyz)

 # Check object is correct (execution stops if not)
 checkCcoord(object)

 return(object)
}
#
# Create a Fcoord object starting from a 3D vector 
createFcoord <- function(xyzf=c(0,0,0),cell=createUnitCell())
{
 # Check UnitCell object is correct
 checkUnitCell(cell)

 # Create object
 object <- new(Class="Fcoord",xyzf=xyzf,cell=cell)

 # Check object is correct (execution stops if not)
 checkFcoord(object)

 return(object)
}
#
# Create a AtomSymbol object starting from a character string
createAtomSymbol <- function(symbol=" C")
{
 # Create object
 object <- new(Class="AtomSymbol",symbol=symbol)

 # Check object is correct (execution stops if not)
 checkAtomSymbol(object)

 return(object)
}
#
# Create an Atom object starting from a character string, a 3D vector and a number
# (bfac iso) or a 3X3 symmetric matrix (bfac aniso)
createAtom <- function(symbol=" C",xyz=c(0,0,0),bfac=0.0)
{
 # Create AtomSymbol object
 symbol <- createAtomSymbol(symbol)

 # Create Ccoord object
 xyz <- createCcoord(xyz=xyz)

 # Create 3X3 diagonal matrix if bfac is a number
 if (is.null(dim(bfac))) bfac <- matrix(c(bfac,0,0,0,bfac,0,0,0,bfac),ncol=3)

 # Create object
 object <- new(Class="Atom",symbol=symbol,xyz=xyz,bfac=bfac)

 # Check object is correct (execution stops if not)
 checkAtom(object)

 return(object)
}
#
# Create an AtomInCell object starting from a character string, a UnitCell object,
# a 3D vector of fractional coordinates and a number
# (bfac iso) or a 3X3 symmetric matrix (bfac aniso)
createAtomInCell <- function(symbol=" C",cell=createUnitCell(),xyzf=c(0,0,0),bfac=0.0)
{
 # Create AtomSymbol object
 symbol <- createAtomSymbol(symbol)

 # Create Fcoord object
 xyzf <- createFcoord(xyzf=xyzf,cell=cell)

 # Create 3X3 diagonal matrix if bfac is a number
 if (is.null(dim(bfac))) bfac <- matrix(c(bfac,0,0,0,bfac,0,0,0,bfac),ncol=3)

 # Create object
 object <- new(Class="AtomInCell",symbol=symbol,xyzf=xyzf,bfac=bfac)

 # Check object is correct (execution stops if not)
 checkAtomInCell(object)

 return(object)
}
