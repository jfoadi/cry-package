# R code to implement ideas related to crystallographic symmetry.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Classes
#
# Symmetry
setClass(
         Class="Symmetry",
         representation=representation(sym_number="numeric",     # We don't use "integer" because, for example, is.integer(1) is FALSE!
                                       setting="numeric",
                                       sym_xHM="character",
                                       cr_sys="character",
                                       xyz="character",
                                       cenop="character",
                                       PG="list",
                                       TR="list",
                                       CT="list")
        )

## Default methods
#
# Print
setMethod(
          f="print",
          signature="Symmetry",
          definition=function(x,...)
                     {
                      # Check Symmetry object is correct
                      checkSymmetry(x)

                      # Print
                      cat("*** This is a Symmetry object ***\n")
                      stringa <- sprintf("    The space group extended Hermann-Maguin symbol is %s\n",x@sym_xHM)
                      cat(stringa)
                      stringa <- sprintf("    The space group number from the International Tables is %d\n",x@sym_number)
                      cat(stringa)
                      stringa <- sprintf("    This space group belongs to the %s crystal system\n",x@cr_sys)
                      cat(stringa)
                      stringa <- sprintf("    This space group includes the following symmetry operations:\n")
                      cat(stringa)
                      for (i in 1:length(x@xyz))
                      {
                       stringa <- sprintf("    %2d) %s\n",i,x@xyz[i])
                       cat(stringa)
                      }
                      if (length(x@cenop) == 1)
                      {
                       stringa <- sprintf("    There are no centring operators associated with this group\n")
                       cat(stringa)
                      }
                      if (length(x@cenop) != 1)
                      {
                       stringa <- sprintf("    This space group is associated with a centred cell. The centring operators are:\n")
                       cat(stringa)
                       for (i in 2:length(x@cenop))
                       {
                        stringa <- sprintf("    %2d) %s\n",(i-1),x@cenop[i])
                        cat(stringa)
                       }
                      }
                     }
         )

## Generic methods
#
# Check Symmetry object is correct
setGeneric(
           name="checkSymmetry",
           def=function(object,message=FALSE){standardGeneric("checkSymmetry")}
          )
#
setMethod(
          f="checkSymmetry",
          signature="Symmetry",
          definition=function(object,message=FALSE)
                     {
                      # If object has been created empty, fine.
                      if (length(object@sym_number) == 0 | length(object@setting) == 0 | length(object@sym_xHM) == 0 | length(object@cr_sys) == 0 | length(object@xyz) == 0 |
                          length(object@cenop) == 0 | length(object@PG) == 0 | length(object@TR) == 0 | length(object@CT) == 0)
                      {
                       warning("This object of class Symmetry is fully or partly empty")

                       # All tests passed: green light!
                       if (message) print("This object and its slots appear to be correct")

                       return()
                      }

                      # Using .translate_SG check extended Hermann-Maguin symbol agrees with number
                      #if (object@setting != 0 & .translate_SG(object@sym_number,setting=object@setting) != object@sym_xHM)
                      # stop("Number and extended Hermann-Maguin symbol for this space group object do not agree")

                      # Check crystal system is correct
                      if(object@cr_sys != .crystal_system(object@sym_number)) stop("The crystal system associated to space group of this object is wrong")

                      # Generate vector of symmetry operations in xyz format and check they coincide with those in object. The same for cenop
                      tmp <- .syminfo_to_op_xyz_list(object@sym_xHM)
                      if (length(tmp[[1]]) != length(object@xyz)) stop("List of symmetry operators in character format do not appear to belong to the correct space group")
                      if (length(tmp[[1]]) == length(object@xyz))
                      { 
                       if (sum(tmp[[1]] == object@xyz) != length(object@xyz))
                           stop("List of symmetry operators in character format do not appear to belong to the correct space group")
                      }
                      if (length(tmp[[2]]) != length(object@cenop)) stop("List of centring operators in character format do not appear to belong to the correct space group")
                      if (length(tmp[[2]]) == length(object@cenop))
                      {
                       if (sum(tmp[[2]] == object@cenop) != length(object@cenop))
                           stop("List of centring operators in character format do not appear to belong to the correct space group")
                      }

                      # Check matrices of point group operators are correct
                      tmp <- .syminfo_to_matrix_list(object@sym_xHM)
                      if (length(tmp$PG) != length(object@PG)) stop("List of point group matrices does not seem to match object symmetry")
                      if (length(tmp$PG) == length(object@PG))
                      {
                       for (i in 1:length(object@PG))
                       {
                        if (sum(tmp$PG[[i]] == object@PG[[i]]) != 9) stop("List of point group matrices does not seem to match object symmetry")
                       }
                      }

                      # Check translation vectors are correct
                      if (length(tmp$T) != length(object@TR)) stop("List of vectors for translational symmetry do not seem to match object symmetry")
                      if (length(tmp$T) == length(object@TR))
                      {
                       for (i in 1:length(object@TR))
                       {
                        if (sum(object@TR[[i]] == tmp$T[[i]]) != 3) stop("List of vectors for translational symmetry do not seem to match object symmetry")
                       }
                      }

                      # Check centring vectors are correct
                      if (length(tmp$C) != length(object@CT)) stop("List of vectors for centring do not seem to match object symmetry")
                      if (length(tmp$C) == length(object@CT))
                      {
                       for (i in 1:length(object@CT))
                       {
                        if (sum(object@CT[[i]] == tmp$C[[i]]) != 3) stop("List of vectors for centring do not seem to match object symmetry")
                       }
                      }

                      # All tests passed: green light!
                      if (message) print("This object and its slots appear to be correct")
                     }
         )
#
# Get symmetry number and extended Herman-Maguin symbol
setGeneric(
           name="getSymmetryNN",
           def=function(object){standardGeneric("getSymmetryNN")}
          )
#
setMethod(
          f="getSymmetryNN",
          signature="Symmetry",
          definition=function(object)
                     {
                      NN <- list(sym_number=object@sym_number,sym_xHM=object@sym_xHM)

                      return(NN)
                     }
         )
#
# Extract symmetry matrices
setGeneric(
           name="extractSymmetryOperators",
           def=function(object){standardGeneric("extractSymmetryOperators")}
          )
#
setMethod(
          f="extractSymmetryOperators",
          signature="Symmetry",
          definition=function(object)
                     {
                      # First check object is of class Symmetry
                      checkSymmetry(object)

                      # Extract interested slots 
                      PG <- object@PG
                      TR <- object@TR

                      return(list(PG=PG,TR=TR))
                     }
         )
#
# Check symmetry is compatible with given cell (return Lattice object)
setGeneric(
           name="checkSymmetryWithCell",
           def=function(object,cell){standardGeneric("checkSymmetryWithCell")}
          )
#
setMethod(
          f="checkSymmetryWithCell",
          signature="Symmetry",
          definition=function(object,cell)
                     {
                      # Check Symmetry object is correct
                      checkSymmetry(object)

                      # Check UnitCell object is correct
                      checkUnitCell(cell)

                      # Extract crystal system associated with symmetry
                      cr_sys <- object@cr_sys

                      # Extract centring associated with space group
                      sl <- substr(object@sym_xHM,1,1) 

                      # Build BravaisType object
                      if (cr_sys == "TRICLINIC") fl <- "a"
                      if (cr_sys == "MONOCLINIC") fl <- "m"
                      if (cr_sys == "ORTHOROMBIC") fl <- "o"
                      if (cr_sys == "TETRAGONAL") fl <- "t"
                      if (cr_sys == "CUBIC") fl <- "c"
                      if (cr_sys == "HEXAGONAL") fl <- "h"
                      if (cr_sys == "TRIGONAL") fl <- "h"
                      bl <- createBravaisType(paste(fl,sl,sep=""))

                      # Now build Lattice object using UnitCell and BravaisType objects. If incompatible, check will fail
                      tmplatt <- createLattice(cell,bl) 

                      # A few more incompatibility conditions remain, connected with the differences in symmetry settings

                      # Extract cell parameters
                      a <- cell@a
                      b <- cell@b
                      c <- cell@c
                      alpha <- cell@alpha
                      if (alpha@rad_flag) alpha <- radToDeg(alpha)
                      aa <- alpha@ang
                      beta <- cell@beta
                      if (beta@rad_flag) beta <- radToDeg(beta)
                      bb <- beta@ang
                      gamma <- cell@gamma
                      if (gamma@rad_flag) gamma <- radToDeg(gamma)
                      cc <- gamma@ang
                      erang <- 0.000001      # Finite accuracy of binary representation of numbers means we have to
                      diff_a <- abs(aa-90)   # test number == 90 in this way
                      diff_b <- abs(bb-90)
                      diff_c <- abs(cc-90)
                      
                      # Extract sym_number and setting
                      sym_number <- object@sym_number
                      setting <- object@setting

                      # Monoclinic
                      if (sym_number == 3 & setting == 1) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 3 & setting == 2) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 3 & setting == 3) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 4 & setting == 1) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 4 & setting == 2) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 4 & setting == 3) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 5 & setting >= 1 & setting <= 3) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 5 & setting >= 10 & setting <= 11) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 5 & setting >= 4 & setting <= 6) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 5 & setting >= 7 & setting <= 9) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 6 & setting == 1) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 6 & setting == 2) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 6 & setting == 3) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 7 & setting >= 1 & setting <= 3) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 7 & setting >= 4 & setting <= 6) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 7 & setting >= 7 & setting <= 9) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 8 & setting >= 1 & setting <= 3) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 8 & setting >= 4 & setting <= 6) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 8 & setting == 10) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 8 & setting >= 7 & setting <= 9) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 9 & setting >= 1 & setting <= 6) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 9 & setting >= 7 & setting <= 12) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 9 & setting >= 13 & setting <= 18) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 10 & setting == 1) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 10 & setting == 2) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 10 & setting == 3) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 11 & setting == 1) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 11 & setting == 2) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 11 & setting == 3) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 12 & setting >= 1 & setting <= 3) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 12 & setting >= 4 & setting <= 6) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 12 & setting >= 7 & setting <= 9) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 13 & setting >= 1 & setting <= 3) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 13 & setting >= 4 & setting <= 6) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 13 & setting >= 7 & setting <= 9) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 14 & setting >= 1 & setting <= 3) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 14 & setting >= 4 & setting <= 6) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 14 & setting >= 7 & setting <= 9) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      if (sym_number == 15 & setting >= 1 & setting <= 6) if (diff_a > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 15 & setting >= 7 & setting <= 12) if (diff_a > erang | diff_b > erang) stop("Symmetry setting incompatible with cell angles")
                      if (sym_number == 15 & setting >= 13 & setting <= 18) if (diff_b > erang | diff_c > erang) stop("Symmetry setting incompatible with cell angles")

                      # Orthorombic

                      return(tmplatt)
                     }
         )


## Friendly constructors
#
# Create a Symmetry object starting from space group number or xHM character symbol
createSymmetry <- function(SG=1,setting=1)
{
 # If SG is an xHM character symbol and also setting different from one has been input
 # just use setting 1
 if (is.character(SG) & setting != 1) warning("Input is an xHM symbol; there is only one setting for this case")

 # If SG is a number, turn it into an xHM character symbol
 if (is.numeric(SG)) SG <- .translate_SG(SG,setting=setting)

 # If SG is an xHM character symbol, check it is properly formatted
 if (is.character(SG))
 {
  tmp <- .translate_SG(SG,SG_in="xHM",SG_out="xHM")
  if (SG != tmp) stop("Something wrong in the derivation of this space group. Please report to J. Foadi")
  sym_number <- .translate_SG(SG,SG_in="xHM",SG_out="number")
  sym_xHM <- SG
 }
 
 # Stop if input is not appropriate
 if (!is.character(SG) & !is.numeric(SG)) stop("Inappropriate input. It can only be the space group number or extended Hermann-Maguin symbol")

 # Find out crystal system
 cr_sys <- .crystal_system(sym_number)

 # Extract symmetry operators in xyz format
 symxyz <- .syminfo_to_op_xyz_list(SG)

 # Extract symmetry operators
 symops <- .syminfo_to_matrix_list(SG)

 # Create new Symmetry object
 object <- new(Class="Symmetry",sym_number=sym_number,setting=setting,sym_xHM=sym_xHM,cr_sys=cr_sys,
                                xyz=symxyz[[1]],cenop=symxyz[[2]],PG=symops$PG,TR=symops$T,CT=symops$C)

 # Check object is correct

 return(object)
}
#
# Create Symmetry object fromm pdb file
createSymmetryFromPDB <- function(file)
{
 # Load pdb file in named list structure
 lpdb <- .readPDB(file)

 # Extract xHM symbol
 SG <- lpdb$cryst1$SG

 # Case by case modifications to SG to match correct extended Hermann-Maguin symbol
 if (SG == "H 3") SG <- "R 3 :H"

 # Unless something else is told, setting is fixed to 1
 setting <- 1

 # Create Symmetry object
 object <- createSymmetry(SG=SG,setting=setting)

 return(object)
}
#
# Create Symmetry object fromm mtz file
createSymmetryFromMTZ <- function(file)
{
 # Load MTZ file in named list structure
 lmtz <- .readMTZ(file,messages=FALSE)

 # Extract space group old symbol
 stmp <- lmtz$header$SYMINF[[5]]
 nstmp <- nchar(stmp)
 pr <- 2
 se <- nstmp-1
 SG_old <- substring(stmp,pr,se)
 SG_xHM <- .translate_SG(SG_old,SG_in="old",SG_out="xHM")

 # Create Symmetry object using createSymmetry
 object <- createSymmetry(SG=SG_xHM)

 return(object)
}
