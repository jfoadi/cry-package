## Part of cRy package
#
## J. Foadi & D. G. Waterman - 2011
#
## Functions not immediately associated with those in other files

#
# This function creates cRy objects from specific files content.
# Valid files are: PDB, MTZ
# Valid classes are: UnitCell, ReciprocalUnitCell, BravaisType, Lattice, Symmetry
#
createFromFile <- function(fileName,fileType=NULL,objectName=NULL)
{
 # Determine type of file from name extension or fileType
 ltmp <- strsplit(fileName,".",fixed=TRUE)
 ext <- tolower(ltmp[[1]][length(ltmp[[1]])])
 if (ext != "pdb" & ext != "mtz" &
     is.null(fileType)
     )
 {
  messaggio <- paste("This file cannot be recognised from its extension. If possible, use argument [fileType].",
                     "\n  Alternatively, change its name to one with a valid extension.")
  stop(messaggio)
 }
 if (!is.null(fileType))
 {
  ext2 <- tolower(fileType)
  if (ext2 != "pdb" & ext2 != "mtz") stop("fileType not recognised.")
  if (ext2 != ext) stop("Contrasting information between fileName and fileType.")
  ext <- ext2
 }
 # ext correspond to an accepted file extension, now. Can carry on

 # Only following classes accepted:
 # "UnitCell", "ReciprocalUnitCell", "BravaisType", "Lattice", "Symmetry"
 if (!is.null(objectName)) 
 {
  if (objectName != "UnitCell" & objectName != "ReciprocalUnitCell" & objectName != "BravaisType" & objectName != "Lattice" &
      objectName != "Symmetry"
     ) stop("This function does not know how to create this input object.")
 }

 #
 # PDB case
 if (ext == "pdb")
 {
  if (is.null(objectName)) objectName <- "Structure"
  
  # Create UnitCell object
  if (objectName == "UnitCell" | objectName == "ReciprocalUnitCell") object <- createUnitCellFromPDB(fileName)
  if (objectName == "ReciprocalUnitCell") object <- computeReciprocalUnitCell(object)

  # Create BravaisType object
  if (objectName == "BravaisType") object <- createBravaisTypeFromPDB(fileName)

  # Create Lattice object
  if (objectName == "Lattice") object <- createLatticeFromPDB(fileName)

  # Create Symmetry object
  if (objectName == "Symmetry") object <- createSymmetryFromPDB(fileName)

  # Here all cases for which code is still not available
  if (objectName != "UnitCell" & objectName != "ReciprocalUnitCell" & objectName != "BravaisType" & objectName != "Lattice" &
      objectName != "Symmetry"
     )
  {
   warning(paste("No code available to create an objet of class ",objectName,". NULL returned instead.",sep=""))
   object <- NULL
  }
 }
 #
 # MTZ case
 if (ext == "mtz")
 {
  if (is.null(objectName)) objectName <- "MergedReflection"
  
  # Create UnitCell object
  if (objectName == "UnitCell" | objectName == "ReciprocalUnitCell") object <- createUnitCellFromMTZ(fileName)
  if (objectName == "ReciprocalUnitCell") object <- computeReciprocalUnitCell(object)

  # Create BravaisType object
  if (objectName == "BravaisType") object <- createBravaisTypeFromMTZ(fileName)

  # Create Lattice object
  if (objectName == "Lattice") object <- createLatticeFromMTZ(fileName)

  # Create Symmetry object
  if (objectName == "Symmetry") object <- createSymmetryFromMTZ(fileName)

  # Here all cases for which code is still not available
  if (objectName != "UnitCell" & objectName != "ReciprocalUnitCell" & objectName != "BravaisType" & objectName != "Lattice" &
      objectName != "Symmetry"
     )
  {
   warning(paste("No code available to create an objet of class ",objectName,". NULL returned instead.",sep=""))
   object <- NULL
  }
 }

 # Return object
 return(object)
}

## Old cRy functions are now used as auxiliary ones
# 
createUnitCellFromPDB <- function(file)
{
 # Load pdb file into a named list
 lpdb <- .readPDB(file)
 
 # Cell parameters
 cpar <- lpdb$cryst1$cell_par

 # Create new UnitCell object
 object <- unitcell(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])

 return(object)
}
#
# Create a UnitCell starting from an mtz file
createUnitCellFromMTZ <- function(file)
{
 # Load mtz file into a named list
 lmtz <- .readMTZ(file,messages=FALSE)
 
 # Extract cell parameters
 cpar <- lmtz$header$CELL

 # Create new UnitCell object
 object <- unitcell(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])

 return(object)
}
#
# Create a BravaisType object starting from a PDB file
createBravaisTypeFromPDB <- function(file)
{
 # Load pdb file into a named list
 lpdb <- .readPDB(file)

 # Extract SG name in xHM format and translate it in space group number
 SG_name <- lpdb$cryst1$SG
 # There are cases where xHM in syminfo is different from xHM in the PDB file
 SG_name <- .findHM(SG_name)
 lista <- .translate_SG(SG_name,SG_in="xHM",SG_out="number")
 if (!lista$ans) stop("Extended Herman-Maguin symbol from PDB file appears to be wrongly formatted")
 sgn <- lista$msg

 # Build bl string in two steps
 b <- .crystal_system(sgn)
 if (b == "TRICLINIC") b <- "a"
 if (b == "MONOCLINIC") b <- "m"
 if (b == "ORTHOROMBIC") b <- "o"
 if (b == "TETRAGONAL") b <- "t"
 if (b == "TRIGONAL") b <- "h"
 if (b == "HEXAGONAL") b <- "h"
 if (b == "CUBIC") b <- "c"
 l <- substr(SG_name,1,1)
 bl <- paste(b,l,sep="")
 object <- bravaistype(bl)

 return(object)
}
#
# Create a Lattice object starting from a PDB file
createLatticeFromPDB <- function(file)
{
 # Load pdb file into a named list
 lpdb <- .readPDB(file)

 # Extract SG name in xHM format and translate it in space group number
 SG_name <- lpdb$cryst1$SG
 # There are cases where xHM in syminfo is different from xHM in the PDB file
 SG_name <- .findHM(SG_name)
 lista <- .translate_SG(SG_name,SG_in="xHM",SG_out="number")
 if (!lista$ans) stop("Extended Herman-Maguin symbol from PDB file appears to be wrongly formatted")
 sgn <- lista$msg

 # Build bl string in two steps
 b <- .crystal_system(sgn)
 if (b == "TRICLINIC") b <- "a"
 if (b == "MONOCLINIC") b <- "m"
 if (b == "ORTHOROMBIC") b <- "o"
 if (b == "TETRAGONAL") b <- "t"
 if (b == "TRIGONAL") b <- "h"
 if (b == "HEXAGONAL") b <- "h"
 if (b == "CUBIC") b <- "c"
 l <- substr(SG_name,1,1)
 bl <- paste(b,l,sep="")

 # Extract cell parameters
 cpar <- lpdb$cryst1$cell_par

 # Create object of class UnitCell
 cella <- unitcell(a=cpar[1],b=cpar[2],c=cpar[3],alpha=cpar[4],beta=cpar[5],gamma=cpar[6])

 # Create object of class BravaisType
 bl <- bravaistype(bl)

 # Create object of class Lattice
 object <- lattice(cell=cella,bl=bl)

 return(object)
}
#
# Create a Symmetry object starting from a PDB file
createSymmetryFromPDB <- function(file)
{
 # Load pdb file into a named list
 lpdb <- .readPDB(file)

 # Extract SG name in xHM format and translate it in space group number
 SG_name <- lpdb$cryst1$SG
 # There are cases where xHM in syminfo is different from xHM in the PDB file
 SG_name <- .findHM(SG_name)

 # Create object of class Symmetry
 object <- new(Class="Symmetry",sym_xHM=SG_name)

 return(object)
}
#
# Create a BravaisType object starting from an MTZ file
createBravaisTypeFromMTZ <- function(file)
{
 # Load mtz file into a named list
 lmtz <- .readMTZ(file,message=FALSE)

 # Extract SG number
 sgn <- lmtz$header$SYMINF[[4]]
 lista <- .translate_SG(sgn)
 if (!lista$ans)
 {
  cat("*** Space Group symbol read in MTZ file header has formatting problems ***\n")
  stop(lista$msg)
 }
 SG_name <- lista$msg

 # Build bl string in two steps
 b <- .crystal_system(sgn)
 if (b == "TRICLINIC") b <- "a"
 if (b == "MONOCLINIC") b <- "m"
 if (b == "ORTHOROMBIC") b <- "o"
 if (b == "TETRAGONAL") b <- "t"
 if (b == "TRIGONAL") b <- "h"
 if (b == "HEXAGONAL") b <- "h"
 if (b == "CUBIC") b <- "c"
 l <- substr(SG_name,1,1)
 bl <- paste(b,l,sep="")
 object <- bravaistype(bl)

 return(object)
}
#
# Create a Lattice object starting from an MTZ file
createLatticeFromMTZ <- function(file)
{
 # Load mtz file into a named list
 lmtz <- .readMTZ(file,message=FALSE)

 # Extract SG number
 sgn <- lmtz$header$SYMINF[[4]]
 lista <- .translate_SG(sgn)
 if (!lista$ans)
 {
  cat("*** Space Group symbol read in MTZ file header has formatting problems ***")
  stop(lista$msg)
 }
 SG_name <- lista$msg

 # Build bl string in two steps
 b <- .crystal_system(sgn)
 if (b == "TRICLINIC") b <- "a"
 if (b == "MONOCLINIC") b <- "m"
 if (b == "ORTHOROMBIC") b <- "o"
 if (b == "TETRAGONAL") b <- "t"
 if (b == "TRIGONAL") b <- "h"
 if (b == "HEXAGONAL") b <- "h"
 if (b == "CUBIC") b <- "c"
 l <- substr(SG_name,1,1)
 bl <- paste(b,l,sep="")

 # Extract cell parameters
 cpar <- lmtz$header$CELL

 # Create object of class UnitCell
 cella <- unitcell(a=cpar[1],b=cpar[2],c=cpar[3],alpha=cpar[4],beta=cpar[5],gamma=cpar[6])

 # Create object of class BravaisType
 bl <- bravaistype(bl)

 # Create object of class Lattice
 object <- lattice(cell=cella,bl=bl)

 return(object)
}
#
# Create a Symmetry object starting from an MTZ file
createSymmetryFromMTZ <- function(file)
{
 # Load mtz file into a named list
 lmtz <- .readMTZ(file,message=FALSE)

 # Extract SG number
 sgn <- lmtz$header$SYMINF[[4]]
 lista <- .translate_SG(sgn)
 if (!lista$ans)
 {
  cat("*** Space Group symbol read in MTZ file header has formatting problems ***")
  stop(lista$msg)
 }
 SG_name <- lista$msg

 # Create object of class Symmetry
 object <- new(Class="Symmetry",sym_xHM=SG_name)

 return(object)
}


### Methods using more than one class. Thus I can't decide for which class they are methods
#
# Check symmetry is compatible with unit cell parameters
setMethod( 
          f="checkSymmetryWithCell",
          signature=c("Symmetry","UnitCell"),
          definition=function(obSymmetry,obUnitCell)
                     {
                      print("Im here!!!")
                      # Space group number
                      sgn <- getSymmetryNumber(obSymmetry)
                      print(sgn)

                      # Extract crystal system associated with symmetry
                      cr_sys <- .crystal_system(sgn[1])
                      print(cr_sys)

                      # Extract centring associated with space group
                      sl <- substr(obSymmetry@sym_xHM,1,1)
                      print(sl)

                      # Build BravaisType object
                      if (cr_sys == "TRICLINIC") fl <- "a"
                      if (cr_sys == "MONOCLINIC") fl <- "m"
                      if (cr_sys == "ORTHOROMBIC") fl <- "o"
                      if (cr_sys == "TETRAGONAL") fl <- "t"
                      if (cr_sys == "CUBIC") fl <- "c"
                      if (cr_sys == "HEXAGONAL") fl <- "h"
                      if (cr_sys == "TRIGONAL") fl <- "h"
                      bl <- bravaistype(paste(fl,sl,sep=""))
                      print(bl)

                      # Now build Lattice object using UnitCell and BravaisType objects. If incompatible, check will fail
                      tmplatt <- lattice(obUnitCell,bl)
                     }
         )
