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
# Molecule
# This object has only one slot, atoms, but this is a data frame with 13 columns:
# element  x  y  z  occ  bfac  u11  u22  u33  u12  u13  u23  chain
setClass(
         Class="Molecule",
         representation=representation(atoms="data.frame")
        )
#
# Structure
setClass(
         Class="Structure",
         representation=representation(molecules="Molecule",cell="UnitCell",symmetry="Symmetry")
        )

## Default methods
#
# Print Molecule
setMethod(
          f="print",
          signature="Molecule",
          definition=function(x,...)
                     {
                      # Check Molecule object is correct
                      checkMolecule(x)

                      cat("*** This is an object of class Molecule ***\n")

                      # Count number of molecules and chain id
                      lvs <- levels(as.factor(x@atoms$chain))
                      stringa <- sprintf("    There are %d molecules in this object\n",length(lvs))
                      if (length(lvs) == 1) stringa <- sprintf("    There is one molecule in this object\n")
                      cat(stringa)

                      # Summary description for each molecule
                      for (i in seq(along = lvs))
                      {
                       # Name of chain identifier
                       if (lvs[i] == " ") stringa <- sprintf("    Molecule %2d has no identifier\n",i)
                       if (lvs[i] != " ") stringa <- sprintf("    Molecule %2d has chain identifier %s\n",i,lvs[i])
                       cat(stringa)

                       # Extract data frame of individual molecule
                       tmp <- x@atoms[x@atoms$chain == lvs[i],]

                       # Number of atoms in molecule
                       natoms <- length(tmp$element)
                       stringa <- sprintf("      This molecule is composed of %5d atoms\n",natoms)
                       cat(stringa)

                       # Number of atoms for each chemica element
                       cat("      More specifically, in the molecule there are:\n")
                       lvs_el <- levels(as.factor(tmp$element))
                       for (j in seq(along = lvs_el))
                       {
                        tmp2 <- tmp$element[tmp$element == lvs_el[j]]
                        symb <- .strip_blanks(lvs_el[j])
                        el_name <- .ATOMS_data.frame$name[.ATOMS_data.frame$symbol == symb]
                        nel <- length(tmp2)
                        stringa <- sprintf("        %4d atoms of %s\n",nel,el_name)
                        cat(stringa)
                       }
                      }
                     }
         )
#
# Add two Molecule objects (no steric-clashes check)
setMethod(
          f="+",
          signature=c("Molecule","Molecule"),
          definition=function(e1,e2)
                     {
                      # Check both objects are correct
                      checkMolecule(e1)
                      checkMolecule(e2)

                      # Extract chains information
                      ch1 <- getMoleculeChainSymbols(e1)
                      ch2 <- getMoleculeChainSymbols(e2)
                      union_ch <- unique(c(ch1,ch2))
                      not_used <- .ALPHABET
                      for (i in seq(along = union_ch))
                      {
                       idx <- which(not_used != union_ch[i])
                       if (length(idx) != 0) not_used <- not_used[idx] 
                      }

                      # Change all e2 chain identifiers which are identical with e1's ones
                      inext <- 0
                      for (i in seq(along = ch1))
                      {
                       idx <- which(ch2 == ch1[i])
                       if (length(idx) != 0)
                       {
                        inext <- inext+1
                        tmp <- extractMoleculeArray(e2,"chain")
                        tmp[tmp == ch1[i]] <- not_used[inext]
                        e2 <- changeMolecule(e2,chain=tmp)
                       }
                      }

                      # Bind two data frames together
                      tmp <- rbind(e1@atoms,e2@atoms)

                      # Create new Molecule object
                      new_object <- new(Class="Molecule",atoms=tmp)

                      # Check new object of class Molecule is correct
                      checkMolecule(new_object)
                      
                      return(new_object)
                     }
         )
#
# Print Structure
setMethod(
          f="print",
          signature="Structure",
          definition=function(x,...)
                     {
                      # Check Structure object is correct
                      checkStructure(x)

                      cat("*** This is an object of class Structure ***\n")
                      cat("*** It is composed of 3 different objects, an object of class Molecule, an object of class UnitCell and an object of class Symmetry\n")
                      cat("\n")
                      cat("    --Object of class Molecule--\n")
                      print(x@molecules)
                      cat("\n")
                      cat("    --Object of class UnitCell--\n")
                      print(x@cell)
                      cat("\n")
                      cat("    --Object of class Symmetry--\n")
                      print(x@symmetry)
                     }
         )

## Generic methods
#
# Check Molecule object is correct
setGeneric(
           name="checkMolecule",
           def=function(object,message=FALSE){standardGeneric("checkMolecule")}
          )
#
setMethod(
          f="checkMolecule",
          signature="Molecule",
          definition=function(object,message=FALSE)
                     {
                      # Warning message if object is fully or partially empty
                      if (length(object@atoms) == 0)
                      {
                       warning("This object of class Structure is fully or partly empty")

                       # All tests passed: green light!
                       if (message) print("This object and its slots appear to be correct")

                       return()
                      }

                      # Check number of columns in data frame
                      tmp <- object@atoms
                      namecols <- colnames(tmp)
                      if (length(namecols) != 13) stop("Slot atoms of object Molecule has the wrong number of columns")

                      # Check all columns have correct name
                      if (namecols[1] != "element") stop("Column  1 of atoms slot of Molecule object has the wrong name")
                      if (namecols[2] != "x")       stop("Column  2 of atoms slot of Molecule object has the wrong name")
                      if (namecols[3] != "y")       stop("Column  3 of atoms slot of Molecule object has the wrong name")
                      if (namecols[4] != "z")       stop("Column  4 of atoms slot of Molecule object has the wrong name")
                      if (namecols[5] != "occ")     stop("Column  5 of atoms slot of Molecule object has the wrong name")
                      if (namecols[6] != "bfac")    stop("Column  6 of atoms slot of Molecule object has the wrong name")
                      if (namecols[7] != "u11")     stop("Column  7 of atoms slot of Molecule object has the wrong name")
                      if (namecols[8] != "u22")     stop("Column  8 of atoms slot of Molecule object has the wrong name")
                      if (namecols[9] != "u33")     stop("Column  9 of atoms slot of Molecule object has the wrong name")
                      if (namecols[10] != "u12")    stop("Column 10 of atoms slot of Molecule object has the wrong name")
                      if (namecols[11] != "u13")    stop("Column 11 of atoms slot of Molecule object has the wrong name")
                      if (namecols[12] != "u23")    stop("Column 12 of atoms slot of Molecule object has the wrong name")
                      if (namecols[13] != "chain")  stop("Column 13 of atoms slot of Molecule object has the wrong name")

                      # Check all columns have correct type
                      if (!is.character(tmp[,1]))  stop("Column  1 of atoms slot of Molecule object is supposed to be a character string")
                      if (!is.numeric(tmp[,2]))    stop("Column  2 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,3]))    stop("Column  3 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,4]))    stop("Column  4 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,5]))    stop("Column  5 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,6]))    stop("Column  6 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,7]))    stop("Column  7 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,8]))    stop("Column  8 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,9]))    stop("Column  9 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,10]))   stop("Column 10 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,11]))   stop("Column 11 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.numeric(tmp[,12]))   stop("Column 12 of atoms slot of Molecule object is supposed to be a numeric string")
                      if (!is.character(tmp[,13])) stop("Column 13 of atoms slot of Molecule object is supposed to be a character string")

                      # All tests passed: green light!
                      if (message) print("This object and its slots appear to be correct")
                     }
         )
#
# Check Structure object is correct
setGeneric(
           name="checkStructure",
           def=function(object,message=FALSE){standardGeneric("checkStructure")}
          )
#
setMethod(
          f="checkStructure",
          signature="Structure",
          definition=function(object,message=FALSE)
                     {
                      # Check slots are correct objects of corresponding class
                      slots <- unname(getSlots(class(object)))
                      if (length(slots) != 3) stop("This is not an object of class Structure, as it hasn't got 3 slots")
                      if (slots[1] != "Molecule" | slots[2] != "UnitCell" | slots[3] != "Symmetry")
                             stop("This is not an object of class Structure, as its slots should be of class Molecule, UnitCell and Symmetry")

                      # Check individual slots are objects of correct class
                      checkMolecule(object@molecules)
                      checkUnitCell(object@cell)
                      checkSymmetry(object@symmetry)

                      # Add future checks

                      # All tests passed: green light!
                      if (message) print("This object and its slots appear to be correct")          
                     }
         )

#
# Change Molecule components
setGeneric(
           name="changeMolecule",
           def=function(object,element,x,y,z,occ,bfac,u11,u22,u33,u12,u13,u23,chain){standardGeneric("changeMolecule")}
          )
#
setMethod(
          f="changeMolecule",
          signature="Molecule",
          definition=function(object,element,x,y,z,occ,bfac,u11,u22,u33,u12,u13,u23,chain)
                     {
                      # Check object is correct
                      checkMolecule(object)

                      # Length of all columns in object
                      len_obj <- length(object@atoms$element)

                      # Replace values with new ones if available
                      if (!missing(element))
                      {
                       if (length(element) == len_obj) object@atoms$element <- element
                       if (length(element) != len_obj) stop("Input array element has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(x))
                      {
                       if (length(x) == len_obj) object@atoms$x <- x
                       if (length(x) != len_obj) stop("Input array x has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(y))
                      {
                       if (length(y) == len_obj) object@atoms$y <- y
                       if (length(y) != len_obj) stop("Input array y has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(z))
                      {
                       if (length(z) == len_obj) object@atoms$z <- z
                       if (length(z) != len_obj) stop("Input array z has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(occ))
                      {
                       if (length(occ) == len_obj) object@atoms$occ <- occ
                       if (length(occ) != len_obj) stop("Input array occ has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(bfac))
                      {
                       if (length(bfac) == len_obj) object@atoms$bfac <- bfac
                       if (length(bfac) != len_obj) stop("Input array bfac has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(u11))
                      {
                       if (length(u11) == len_obj) object@atoms$u11 <- u11
                       if (length(u11) != len_obj) stop("Input array u11 has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(u22))
                      {
                       if (length(u22) == len_obj) object@atoms$u22 <- u22
                       if (length(u22) != len_obj) stop("Input array u22 has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(u33))
                      {
                       if (length(u33) == len_obj) object@atoms$u33 <- u33
                       if (length(u33) != len_obj) stop("Input array u33 has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(u12))
                      {
                       if (length(u12) == len_obj) object@atoms$u12 <- u12
                       if (length(u12) != len_obj) stop("Input array u12 has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(u13))
                      {
                       if (length(u13) == len_obj) object@atoms$u13 <- u13
                       if (length(u13) != len_obj) stop("Input array u13 has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(u23))
                      {
                       if (length(u23) == len_obj) object@atoms$u23 <- u23
                       if (length(u23) != len_obj) stop("Input array u23 has length differing from length of corresponding slot in object Molecule")
                      }
                      if (!missing(chain))
                      {
                       if (length(chain) == len_obj) object@atoms$chain <- chain
                       if (length(chain) != len_obj) stop("Input array chain has length differing from length of corresponding slot in object Molecule")
                      }

                      # Check modified object is OK
                      checkMolecule(object)

                      return(object)
                     }
         )
#
# Extract Molecule components
setGeneric(
           name="extractMoleculeArray",
           def=function(object,achoice){standardGeneric("extractMoleculeArray")}
          )
#
setMethod(
          f="extractMoleculeArray",
          signature="Molecule",
          definition=function(object,achoice)
                     {
                      # If achoice is missing, quit
                      if (missing(achoice))
                      {
                       warning("You haven't selected any array from object Molecule")
                       return(NULL)
                      }

                      # Quit if achoice is not an array of object Molecule
                      if (!is.character(achoice)) stop("Input variable achoice of this method for object Molecule can only be a character, the name of one of Molecule's arrays")
                      if (achoice != "element" & achoice != "x" & achoice != "y" & achoice != "z" & achoice != "occ" & achoice != "bfac" &
                          achoice != "u11" & achoice != "u22" & achoice != "u33" & achoice != "u12" & achoice != "u13" & achoice != "u23" & achoice != "chain")
                          stop("Input variable achoice of this method for object Molecule can only be the name of one of Molecule's arrays")
                      
                      # OK. Carry on with extraction
                      idx <- which(colnames(object@atoms) == achoice)
                      extracted_obj <- object@atoms[,idx]

                      return(extracted_obj)
                     }
         )
#
# Display given lines of Molecule object
setGeneric(
           name="displayMoleculeRecords",
           def=function(object,nrec=1:10){standardGeneric("displayMoleculeRecords")}
          )
#
setMethod(
          f="displayMoleculeRecords",
          signature="Molecule",
          definition=function(object,nrec=1:10)
                     {
                      # Check input object is a correct object of class Molecule
                      checkMolecule(object)

                      # Check nrec is a numeric
                      if (!is.numeric(nrec)) stop("Input object nrec has to be an integer or a vector of integers")

                      # Display selected lines
                      tmp <- object@atoms[nrec,]
                      print(tmp)
                     }
         )
#
# Extract name of all chains in a Molecule object
setGeneric(
           name="getMoleculeChainSymbols",
           def=function(object){standardGeneric("getMoleculeChainSymbols")}
          )
#
setMethod(
          f="getMoleculeChainSymbols",
          signature="Molecule",
          definition=function(object)
                     {
                      # Check object of class Molecule is correct
                      checkMolecule(object)

                      # Extract chains info
                      lvs <- levels(as.factor(object@atoms$chain))

                      return(lvs)
                     }
         )
#
# Extract individual molecule (chain) from a Molecule object
setGeneric(
           name="extractMoleculeChain",
           def=function(object,chain){standardGeneric("extractMoleculeChain")}
          )
#
setMethod(
          f="extractMoleculeChain",
          signature="Molecule",
          definition=function(object,chain)
                     {
                      # Check object of class Molecule is correct
                      checkMolecule(object)

                      # Check input chain is a character
                      if (!is.character(chain)) stop("Argument chain of this Molecule method must be a character")

                      # Check input chain is among Molecule existing chains
                      chains <- getMoleculeChainSymbols(object)
                      if (sum(chain == chains) != 1) stop("Your input chain is not included among existing chains in object of class Molecule")

                      # Passed al checks. Extract group of interest
                      tmp <- object@atoms[object@atoms$chain == chain,]
                      new_object <- new(Class="Molecule",atoms=tmp)

                      # Check new object is correct
                      checkMolecule(new_object)

                      return(new_object)
                     }
         )
#
# Rotate Molecule object around its centre of mass
setGeneric(
           name="rotateMolecule",
           def=function(object,t1,t2,t3,chi,psi,fi){standardGeneric("rotateMolecule")}
          )
#
setMethod(
          f="rotateMolecule",
          signature="Molecule",
          definition=function(object,t1,t2,t3,chi,psi,fi)
                     {
                      # Check object of class Molecule is correct
                      checkMolecule(object)

                      # Create rotation matrix using Euler or polar angles
                      boh <- 0
                      if (missing(chi) & missing(psi) & missing(fi))
                      {
                       boh <- 1
                       Rm <- .euler_to_matrix(t1,t2,t3)
                      }
                      if (missing(t1) & missing(t2) & missing(t3))
                      {
                       boh <- 1
                       Rm <- .polar_to_matrix(chi,psi,fi)
                      }

                      # Stop if euler and polar angles have been mixed up
                      if (boh == 0) stop("Object of class Molecule can only be rotated using either Euler or polar angles")

                      # Extract list of atoms position as a nX3 matrix
                      xyz <- unname(as.matrix(object@atoms[,2:4]))

                      # Molecule's centre of mass
                      com <- apply(xyz,2,mean)

                      # Shift molecule to its centre of mass
                      xyz[,1] <- xyz[,1]-com[1]
                      xyz[,2] <- xyz[,2]-com[2]
                      xyz[,3] <- xyz[,3]-com[3]

                      # Rotate molecule
                      new_xyz <- t(Rm%*%t(xyz))

                      # Shift molecule back to initial frame of reference
                      new_xyz[,1] <- new_xyz[,1]+com[1]
                      new_xyz[,2] <- new_xyz[,2]+com[2]
                      new_xyz[,3] <- new_xyz[,3]+com[3]

                      # New molecule object
                      new_object <- changeMolecule(object,x=new_xyz[,1],y=new_xyz[,2],z=new_xyz[,3])

                      return(new_object)
                     }
         )
#
# Translate Molecule object of a fixed amount
setGeneric(
           name="translateMolecule",
           def=function(object,trsl){standardGeneric("translateMolecule")}
          )
#
setMethod(
          f="translateMolecule",
          signature="Molecule",
          definition=function(object,trsl)
                     {
                      # Check object of class Molecule is correct
                      checkMolecule(object)
                      
                      # Check trsl is a 3D vector
                      if (!is.numeric(trsl) | length(trsl) != 3) stop("Input parameter trsl has to be a numeric vector with 3 components")

                      # Extract list of atoms position as a nX3 matrix
                      xyz <- unname(as.matrix(object@atoms[,2:4]))

                      # Translate molecule using vectot trsl
                      xyz[,1] <- xyz[,1]+trsl[1]
                      xyz[,2] <- xyz[,2]+trsl[2]
                      xyz[,3] <- xyz[,3]+trsl[3]

                      # New molecule object
                      new_object <- changeMolecule(object,x=xyz[,1],y=xyz[,2],z=xyz[,3])

                      return(new_object)
                     }
         )
#
# Expand Molecule object according to symmetry. 
# Input is an object of class Molecule, an object of class UnitCell and an object of class Symmetry.
# Output is a list of objects of class Molecule, each molecule being a symmetry equivalent to the input molecule.
setGeneric(
           name="createSymmetryEquivalent",
           def=function(object,cell,SG){standardGeneric("createSymmetryEquivalent")}
          )
#
setMethod(
          f="createSymmetryEquivalent",
          signature="Molecule",
          definition=function(object,cell,SG)
                     {
                      # Check object of class Molecule is correct
                      checkMolecule(object)
                      
                      # Check object of class Symmetry is correct
                      checkSymmetry(SG)

                      # Check object of class UnitCell is correct
                      checkUnitCell(cell)

                      # Extract list of atoms position as a nX3 matrix
                      xyz <- unname(as.matrix(object@atoms[,2:4]))

                      # Convert coordinates to fractional
                      xyzf <- .orth_to_frac(xyz,cell@a,cell@b,cell@c,cell@alpha@ang,cell@beta@ang,cell@gamma@ang)

                      # Extract symmetry operators
                      PG <- SG@PG
                      TR <- SG@TR
                      nsym <- length(PG)
                      if (length(TR) != nsym) stop("NUmber of point group matrices does not match number of translational operators")

                      # Extract centring operators
                      CT <- SG@CT

                      # Create symmetry equivalent
                      symeq <- list()
                      for (i in 1:nsym)
                      {
                       tmp <- t(PG[[i]]%*%t(xyzf))+TR[[i]]
                       symeq <- c(symeq,list(tmp))
                      }

                      # Use the centring, if any
                      if (length(CT) > 1)
                      {
                       for (ic in 2:length(CT))
                       {
                        for (i in seq(along = symeq))
                        {
                         tmp <- symeq[[i]]+CT[[ic]]
                         symeq <- c(symeq,list(tmp))
                        }
                       }
                      }

                      # Convert back to orthogonal coordinates
                      for (i in seq(along = symeq))
                      {
                       tmp <- .frac_to_orth(symeq[[i]],cell@a,cell@b,cell@c,cell@alpha@ang,cell@beta@ang,cell@gamma@ang)
                       symeq[[i]] <- changeMolecule(object,x=tmp[,1],y=tmp[,2],z=tmp[,3])
                      }

                      return(symeq)
                     }
         )

## Friendly constructors
#
# Create a Molecule object starting from a pdb file 
createMoleculeFromPDB <- function(file,chain_sel)
{
 # Load pdb file
 lpdb <- .readPDB(file)
 raw_atoms <- lpdb$Atom

 # If Hetatm exists add it to tmp data frame
 if (sum(names(lpdb) == "Hetatm") == 1) raw_atoms <- rbind(raw_atoms,lpdb$Hetatm)

 # If Anisou exists copy it to a temporary data frame
 if (sum(names(lpdb) == "Anisou") == 1) raw_anisou <- lpdb$Anisou

 # If chain_sel missing fix chain 
 if (missing(chain_sel)) chain_sel <- "@"

 # chain_sel components are characters
 chain_sel <- as.character(chain_sel)
 nchains <- length(chain_sel)
 if (chain_sel[1] == "@") nchains <- 0

 # Check chain_sel characters are included in pdb chainID
 if (nchains != 0)
 {
  lvs <- levels(as.factor(raw_atoms$chainID))
  ntimes <- 0
  for (i in seq(along = chain_sel))
  {
   if (sum(chain_sel[i] == lvs) == 1) ntimes <- ntimes+1
  }
  if (ntimes == 0) stop("The chain or group of chains selected for this molecule or group of molecules are not present in the pdb file")
 }

 # Extract only selected chains
 if (chain_sel[1] != "@")
 {
  idx <- which(raw_atoms$chainID == chain_sel[1])
  tmp <- raw_atoms[idx,]
  if (sum(names(lpdb) == "Anisou") == 1) tmp2 <- raw_anisou[idx,]
  if (nchains > 1)
  {
   for (i in 2:nchains)
   {
    idx <- which(raw_atoms$chainID == chain_sel[i])
    tmp <- rbind(tmp,raw_atoms[idx,])
    if (sum(names(lpdb) == "Anisou") == 1) tmp2 <- rbind(tmp2,raw_anisou[idx,])
   }
  }
 }

 # Extract all chains
 if (chain_sel[1] == "@")
 {
  tmp <- raw_atoms
  if (sum(names(lpdb) == "Anisou") == 1) tmp2 <- raw_anisou
 }

 # Count number of atoms in molecule (or molecules)
 natoms <- length(tmp[,1])

 # If Anisou part is missing fill u11,u22,etc with 0's
 if (sum(names(lpdb) == "Anisou") != 1)
 {
  u11 <- rep(0,times=natoms)
  u22 <- rep(0,times=natoms)
  u33 <- rep(0,times=natoms)
  u12 <- rep(0,times=natoms)
  u13 <- rep(0,times=natoms)
  u23 <- rep(0,times=natoms)
 }

 # If Anisou part is not missing fill u11,u22,etc with values from pdb file
 if (sum(names(lpdb) == "Anisou") == 1)
 {
  u11 <- tmp2$u11
  u22 <- tmp2$u22
  u33 <- tmp2$u33
  u12 <- tmp2$u12
  u13 <- tmp2$u13
  u23 <- tmp2$u23
 }

 # Now build data.frame
 tmp <- data.frame(element=I(tmp$atomSymbol),x=tmp$x,y=tmp$y,z=tmp$z,occ=tmp$occ,bfac=tmp$Bfac,
                   u11=u11,u22=u22,u33=u33,u12=u12,u13=u13,u23=u23,chain=I(tmp$chainID))
 

 # Create new object
 object <- new(Class="Molecule",atoms=tmp)

 # Check new object

 return(object)
}
#
# Create a Molecule object starting from individual arrays
createMoleculeFromScratch <- function(element,x,y,z,occ,bfac,u11,u22,u33,u12,u13,u23,chain)
{
 # Check if some arrays are missing (anisotropy u's and chain are not compulsory)
 if (missing(element) | missing(x) | missing(y) | missing(z) | missing(occ) | missing(bfac)) 
     stop("You are trying to create a Molecule object from individual 'element x y z occ bfac' arrays, but some of them are actually missing")

 # Check array element is of the correct type
 if (!is.character(element)) stop("Input element array is of the wrong type. It needs to be a character array")

 # Check all arrrays have the same length
 n <- length(element)
 if (length(x) != n | length(y) != n | length(z) != n | length(occ) != n) stop("Input arrays have differing lengths")

 # Check all components of array element have a maximm of 2 letters
 if (sum(nchar(element) > 2) > 0) stop("One or more of input array element is not an acepted chemical symbol")

 # Check all components of array element are an accepted chemical-element symbol
 for (i in seq(along=element))
 {
  if (sum(.strip_blanks(element[i]) == .ATOMS_data.frame$symbol) != 1)
  {
   stringa <- sprintf("Element %d of the Molecule object is not an accepted chemical element",i)
   stop(stringa)
  }
 }

 # Check x, y, z, occ and bfac are numerics
 if (!is.numeric(x)) stop("Array x does not appear to be numeric")
 if (!is.numeric(y)) stop("Array y does not appear to be numeric")
 if (!is.numeric(z)) stop("Array z does not appear to be numeric")
 if (!is.numeric(occ)) stop("Array occ does not appear to be numeric")
 if (!is.numeric(bfac)) stop("Array bfac does not appear to be numeric")

 # Check occupancy is between 0 and 1
 if (sum(occ < 0 | occ > 1) > 0) stop("Array occ appears to have some values smaller than 0 or greater than 1")

 # If some u's are not in the input, put corresponding array to 0
 if (missing(u11)) u11 <- rep(0,times=n)
 if (missing(u22)) u22 <- rep(0,times=n)
 if (missing(u33)) u33 <- rep(0,times=n)
 if (missing(u12)) u12 <- rep(0,times=n)
 if (missing(u13)) u13 <- rep(0,times=n)
 if (missing(u23)) u23 <- rep(0,times=n)

 # Stop if u's are not numeric
 if (!is.numeric(u11)) stop("Input array u11 does not appear to be numeric")
 if (!is.numeric(u22)) stop("Input array u22 does not appear to be numeric")
 if (!is.numeric(u33)) stop("Input array u33 does not appear to be numeric")
 if (!is.numeric(u12)) stop("Input array u12 does not appear to be numeric")
 if (!is.numeric(u13)) stop("Input array u13 does not appear to be numeric")
 if (!is.numeric(u23)) stop("Input array u23 does not appear to be numeric")

 # Stop if some or all u's are not of length n
 if (length(u11) != n) stop("Length of array u11 does not match length of other input arrays")
 if (length(u22) != n) stop("Length of array u22 does not match length of other input arrays")
 if (length(u33) != n) stop("Length of array u33 does not match length of other input arrays")
 if (length(u12) != n) stop("Length of array u12 does not match length of other input arrays")
 if (length(u13) != n) stop("Length of array u13 does not match length of other input arrays")
 if (length(u23) != n) stop("Length of array u23 does not match length of other input arrays")

 # If chain is not included in input fix it to " "
 if (missing(chain)) chain <- rep(" ",times=n)

 # Check existing chain has correct length
 if (length(chain) != n) stop("Array chain has length different from other input arrays")

 # If existing chain is numeric, turn it into corresponding character (can't have mixed number-characters vectors, they are not allowed in R)
 if (is.numeric(chain)) chain <- as.character(chain)

 # All checks passed. Now build data.frame
 tmp <- data.frame(element=element,x=x,y=y,z=z,occ=occ,bfac=bfac,u11=u11,u22=u22,u33=u33,u12=u12,u13=u13,u23=u23,chain=I(chain))

 # Create new object
 object <- new(Class="Molecule",atoms=tmp)

 # Check new object
 checkMolecule(object)

 return(object)
}
#
# Create Structure object from objects of class Molecule, UnitCell and Symmetry
createStructureFromMUS <- function(molecules,cell,symmetry)
{
 # Check molecules is of class Molecule
 checkMolecule(molecules)

 # Check cell is of class UnitCell
 checkUnitCell(cell)

 # Check symmetry is of class Symmetry
 checkSymmetry(symmetry)

 # Create new object
 object <- new(Class="Structure",molecules=molecules,cell=cell,symmetry=symmetry)

 # Check object of class Structure is correct
 checkStructure(object)

 return(object)
}
#
# Create Structure object from PDB
createStructureFromPDB <- function(file)
{
 # Create Molecule object from pdb
 mol <- createMoleculeFromPDB(file)

 # Create UnitCell object from pdb
 cell <- createUnitCellFromPDB(file)
 
 # Create Symmetry object from pdb
 SG <- createSymmetryFromPDB(file)

 # Create new object
 object <- createStructureFromMUS(mol,cell,SG)

 # Check object of class Structure is correct
 checkStructure(object)

 return(object)
}
