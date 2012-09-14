## Part of cRy package
#
## J. Foadi & D. G. Waterman - 2011
#
## Functions not immediately asociate with those in other files

#
# This function creates cRy objects from specific files content.
# Valid files are: PDB, MTZ
# Valid classes are: UnitCell, ReciprocalUnitCell, BravaisType, Lattice
#
createFromFile <- function(fileName,fileType=NULL,objectName=NULL)
{
 # Determine type of file from name extension or fileType
 ltmp <- strsplit(fileName,".",fixed=TRUE)
 ext <- tolower(ltmp[[1]][length(ltmp[[1]])])
 if (ext != "pdb" & ext != "mtz" &
     is.null(fileType)
     ) stop("This file cannot be recognised from its extension. Consider using argument [fileType].")
 if (!is.null(fileType))
 {
  ext2 <- tolower(fileType)
  if (ext2 != "pdb" & ext2 != "mtz") stop("fileType not recognised.")
  if (ext2 != ext) stop("Contrasting information between fileName and fileType.")
  ext <- ext2
 }
 # ext correspond to an accepted file extension, now. Can carry on

 # Only following classes accepted:
 # "UnitCell", "ReciprocalUnitCell", "BravaisType", "Lattice"
 if (!is.null(objectName)) 
 {
  if (objectName != "UnitCell" & objectName != "ReciprocalUnitCell" & objectName != "BravaisType" & objectName != "Lattice"
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

  # Here all cases for which code is still not available
  if (objectName != "UnitCell" & objectName != "ReciprocalUnitCell" & objectName != "BravaisType" & objectName != "Lattice"
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

  # Here all cases for which code is still not available
  if (objectName != "UnitCell" & objectName != "ReciprocalUnitCell" & objectName != "BravaisType" & objectName != "Lattice"
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
 # These are cases where xHM in syminfo is different from xHM in the PDB file (add as many as needed)
 if (SG_name == "H 3") SG_name <- "R 3 :H"
 sgn <- .translate_SG(SG_name,SG_in="xHM",SG_out="number")

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
 object <- createBravaisType(bl)
 isCorrect(object,"BravaisType")

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
 # These are cases where xHM in syminfo is different from xHM in the PDB file (add as many as needed)
 if (SG_name == "H 3") SG_name <- "R 3 :H"
 sgn <- .translate_SG(SG_name,SG_in="xHM",SG_out="number")

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

 # Create object of class Lattice
 object <- createLattice(a=cpar[1],b=cpar[2],c=cpar[3],alpha=cpar[4],beta=cpar[5],gamma=cpar[6],bl=bl)

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
 object <- createBravaisType(bl)
 isCorrect(object,"BravaisType")

 return(object)
}
#
# Create a Lattice object starting from an MTZ file
createLatticeFromMTZ <- function(file)
{
 # Load pdb file into a named list
 lmtz <- .readMTZ(file,message=FALSE)

 # Extract SG number
 sgn <- lmtz$header$SYMINF[[4]]
 SG_name <- .translate_SG(sgn)

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

 # Create object of class Lattice
 object <- createLattice(a=cpar[1],b=cpar[2],c=cpar[3],alpha=cpar[4],beta=cpar[5],gamma=cpar[6],bl=bl)

 return(object)
}
