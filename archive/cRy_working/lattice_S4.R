# R code to implement the crystallographic ideas related to crystal lattices.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Classes
#
# BravaisType
setClass(
         Class="BravaisType",
         representation=representation(bl="character")
        )
#
# Lattice
setClass(
         Class="Lattice",
         representation=representation(cell="UnitCell",bl="BravaisType")
        )

## Default methods
#
# Print
setMethod(
          f="print",
          signature="Lattice",
          definition=function(x,...)
                     {
                      # Check Lattice object is correct
                      checkLattice(x)

                      # Extract UnitCell object
                      unit_cell <- x@cell

                      # Print
                      cat("*** This is a Lattice object. ***\n")
                      cat("    The unit cell has the following parameters:\n")
                      stringa <- sprintf("       a = %7.2f\n       b = %7.2f\n       c = %7.2f\n   alpha =  %6.2f\n    beta =  %6.2f\n   gamma =  %6.2f\n",
                                         unit_cell@a,unit_cell@b,unit_cell@c,unit_cell@alpha@ang,unit_cell@beta@ang,unit_cell@gamma@ang)
                      cat(stringa)

                      # Extract lattice information
                      latt_info <- extractLatticeStuff(x)

                      # Print
                      stringa <- sprintf("\n    The Bravais lattice is %s",latt_info[1])
                      cat(stringa)
                      stringa <- sprintf("\n    The crystal family is %s",latt_info[2])
                      cat(stringa)
                      stringa <- sprintf("\n    The crystal system is %s",latt_info[3])
                      cat(stringa)
                      stringa <- sprintf("\n    The lattice system is %s\n",latt_info[4])
                      cat(stringa)
                     }
         )

## Generic methods
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
# Check Bravais lattice type
setGeneric(
           name="checkBravaisType",
           def=function(object,message=FALSE){standardGeneric("checkBravaisType")}
          )
#
setMethod(
          f="checkBravaisType",
          signature="BravaisType",
          definition=function(object,message=FALSE)
                     {
                      # Check BravaisType character has length 2
                      if (nchar(object@bl) != 2) stop("Object of class 'BravaisType' has number of characters different from 2")

                      # Extract lattice and centering (also check on symbol validity)
                      latt <- substr(object@bl,1,1)
                      centring <- substr(object@bl,2,2)
                      if (latt != "a" & latt !="m" & latt != "o" & latt != "t" & latt != "h" & latt != "c") 
                                                          stop("Invalid lattice type")
                      if (centring != "P" & centring != "A" & centring != "B" & centring != "C" & centring != "F" &
                          centring != "I" & centring != "R") stop("Invalid centring type")

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

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
#
# Check compatibility of cell with lattice and existence of Bravais lattice type
setGeneric(
           name="checkLattice",
           def=function(object,message=FALSE){standardGeneric("checkLattice")}
          )
#
setMethod(
          f="checkLattice",
          signature="Lattice",
          definition=function(object,message=FALSE)
                     {
                      # Extract slot names 
                      slot_names <- getSlots(class(object))

                      # Check object has 2 slots
                      if (length(slot_names) != 2) stop("Object 'Lattice' has a number of slots different from 2")
                      
                      # Check nature of first slot (it should be UnitCell)
                      if (names(slot_names[1]) != "cell") stop("First slot does not appear to be named 'cell'")
                      checkUnitCell(object@cell) # Stop here if wrong UnitCell object

                      # Check nature of second slot (it should be BravaisType)
                      if (slot_names[2] != "BravaisType") stop("Second slot is not an object of class 'BravaisType'")
                      checkBravaisType(object@bl)

                      # Now proceed with further checks
                      # Extract cell parameters
                      a <- object@cell@a
                      b <- object@cell@b
                      c <- object@cell@c
                      alpha <- object@cell@alpha
                      if (alpha@rad_flag) alpha <- radToDeg(alpha)
                      aa <- alpha@ang
                      beta <- object@cell@beta
                      if (beta@rad_flag) beta <- radToDeg(beta)
                      bb <- beta@ang
                      gamma <- object@cell@gamma
                      if (gamma@rad_flag) gamma <- radToDeg(gamma)
                      cc <- gamma@ang

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
                      if (latt == "h")                                    # Hexagonal family
                      {
                       if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type") 
                       if (abs(a-b) < 0.000001)
                       {
                        if (abs(cc-120) > 0.000001)
                        {
                         if (abs(b-c) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
                        }
                       }
                      }

                      # Green light: object appears to be correct
                      if (message) print("This object and its slots appear to be correct")
                     }
          )
#
# Change/Assign lattice parameters
setGeneric(
           name="changeLatticeParameters",
           def=function(object,cell=NULL,bl=NULL){standardGeneric("changeLatticeParameters")}
          )
#
setMethod(
          f="changeLatticeParameters",
          signature="Lattice",
          definition=function(object,cell=NULL,bl=NULL)
                     {
                      # Change unit cell parameters
                      if (class(cell) != "NULL")
                      {
                       # Check new unit cell is OK
                       checkUnitCell(cell)

                       # Now replace old unit cell with new unit cell
                       object@cell <- cell
                      }

                      # Change Bravais lattice
                      if (class(bl) != "NULL")
                      {
                       # Check new Bravais lattice is OK
                       checkBravaisType(bl)

                       # Now replace old Bravais lattice with new one
                       object@bl <- bl
                      } 

                      # Check modified lattice is OK
                      checkLattice(object)

                      # Green light: everything OK. Return modified object
                      return(object)
                     }
          )
#
# Returns centering operators
setGeneric(
           name="getCentringOps",
           def=function(object){standardGeneric("getCentringOps")} 
          )
#
setMethod(
          f="getCentringOps",
          signature="Lattice",
          definition=function(object){
                                      # Check object is of class Lattice (execution halted if it is not)
                                      checkLattice(object)

                                      # Centring operators are a list of 3D vectors

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
         )

## Friendly constructors
#
# Create a BravaisType onject starting from a string
createBravaisType <- function(bl="aP")
{
 # Build object without thinking if it's correct
 object <- new(Class="BravaisType",bl=bl)

 # Now check if object just built is correct
 checkBravaisType(object)

 # If things are OK, then return object
 return(object)
}
#
# Assign values to a Lattice object while creating it
createLattice <- function(unit_cell,bl=new(Class="BravaisType",bl="aP"))
{
 # Check input objectis of class UnitCell
 checkUnitCell(unit_cell)

 # Check Bravais lattice type
 checkBravaisType(bl)

 # If everything ok merge two objects into one
 object <- new(Class="Lattice",cell=unit_cell,bl=bl)

 # Check compatibility between cell and lattice, and other small things
 checkLattice(object)

 return(object)
}
