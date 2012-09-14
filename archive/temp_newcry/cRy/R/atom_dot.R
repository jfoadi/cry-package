# For use in conjunction with S4-classes based modules in cRy.
# The word "dot" has been adopted because all functions here start with a ".", so
# to be invisible.

.readPDB <- function(filename,message=FALSE)
{
 # First read PDB ascii file and store it
 pdb0 <- scan(filename,what="character",sep="\n",quiet=TRUE)

 # Output object is a named list
 lpdb <- list()

 # Extract header
 value <- "HEADER"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$header <- c(substr(pdb0[bsg],11,50),substr(pdb0[bsg],51,59),substr(pdb0[bsg],63,66))
 if (length(bsg) != 1) lpdb$header <- NULL 
 
 # Extract cell parameters and space group
 value <- "CRYST1"
 bsg <- grep(value,pdb0,fixed=TRUE)
 aa <- as.numeric(substr(pdb0[bsg],7,15))
 bb <- as.numeric(substr(pdb0[bsg],16,24))
 cc <- as.numeric(substr(pdb0[bsg],25,33))
 alpha <- as.numeric(substr(pdb0[bsg],34,40))
 beta  <- as.numeric(substr(pdb0[bsg],41,47))
 gamma <- as.numeric(substr(pdb0[bsg],48,54))
 cell_par <- c(aa,bb,cc,alpha,beta,gamma)
 ltmp <- substr(pdb0[bsg],56,66)

 # Trim empty space at the end of SG string
 nempty <- 0
 while (substr(ltmp,nchar(ltmp)-nempty,nchar(ltmp)-nempty) == " ") nempty <- nempty+1
 lSGroup <- substr(ltmp,1,nchar(ltmp)-nempty)
 z <- as.numeric(substr(pdb0[bsg],67,70))
 if (length(bsg) == 1) lpdb$cryst1 <- list(cell_par=cell_par,SG=lSGroup,z=z)
 if (length(bsg) != 1) lpdb$cryst1 <- NULL

 # Extract scales
 value <- "SCALE1"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$scale1 <- c(as.numeric(substr(pdb0[bsg],11,20)),as.numeric(substr(pdb0[bsg],21,30)),as.numeric(substr(pdb0[bsg],31,40)),
                                        as.numeric(substr(pdb0[bsg],46,55)))
 if (length(bsg) != 1) lpdb$scale1 <- NULL
 value <- "SCALE2"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$scale2 <- c(as.numeric(substr(pdb0[bsg],11,20)),as.numeric(substr(pdb0[bsg],21,30)),as.numeric(substr(pdb0[bsg],31,40)),
                                        as.numeric(substr(pdb0[bsg],46,55)))
 if (length(bsg) != 1) lpdb$scale2 <- NULL
 value <- "SCALE3"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$scale3 <- c(as.numeric(substr(pdb0[bsg],11,20)),as.numeric(substr(pdb0[bsg],21,30)),as.numeric(substr(pdb0[bsg],31,40)),
                                        as.numeric(substr(pdb0[bsg],46,55)))
 if (length(bsg) != 1) lpdb$scale3 <- NULL

 # Extract origin shifts
 value <- "ORIGX1"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$origx1 <- c(as.numeric(substr(pdb0[bsg],11,20)),as.numeric(substr(pdb0[bsg],21,30)),as.numeric(substr(pdb0[bsg],31,40)),
                                        as.numeric(substr(pdb0[bsg],46,55)))
 if (length(bsg) != 1) lpdb$origx1 <- NULL
 value <- "ORIGX2"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$origx2 <- c(as.numeric(substr(pdb0[bsg],11,20)),as.numeric(substr(pdb0[bsg],21,30)),as.numeric(substr(pdb0[bsg],31,40)),
                                        as.numeric(substr(pdb0[bsg],46,55)))
 if (length(bsg) != 1) lpdb$origx2 <- NULL
 value <- "ORIGX3"
 bsg <- grep(value,pdb0,fixed=TRUE)
 if (length(bsg) == 1) lpdb$origx3 <- c(as.numeric(substr(pdb0[bsg],11,20)),as.numeric(substr(pdb0[bsg],21,30)),as.numeric(substr(pdb0[bsg],31,40)),
                                        as.numeric(substr(pdb0[bsg],46,55)))
 if (length(bsg) != 1) lpdb$origx3 <- NULL

 # Atoms: atom names, residue name, atomic coordinates, etc.
 # Raw data are at present organized in a data frame
 value <- "ATOM  "
 bsg <- which(substr(pdb0,1,6) == value)                    # Done slightly different way as there are other "ATOM  " parts of string in the middle of lines
 atom_name <- substr(pdb0[bsg],13,16)
 residue_name <- substr(pdb0[bsg],18,20)
 structure <- data.frame(          n = as.integer(substr(pdb0[bsg],7,11)),
                            atomName = I(atom_name),
                              altLoc = I(substr(pdb0[bsg],17,17)),
                         residueName = I(residue_name),
                             chainID = I(substr(pdb0[bsg],22,22)),
                                resN = as.integer(substr(pdb0[bsg],23,26)),
                               iCode = I(substr(pdb0[bsg],27,27)),
                                   x = as.numeric(substr(pdb0[bsg],31,38)),
                                   y = as.numeric(substr(pdb0[bsg],39,46)),
                                   z = as.numeric(substr(pdb0[bsg],47,54)),
                                 occ = as.numeric(substr(pdb0[bsg],55,60)),
                                Bfac = as.numeric(substr(pdb0[bsg],61,66)),
                          atomSymbol = I(substr(pdb0[bsg],77,78)),
                          atomCharge = I(substr(pdb0[bsg],79,80)))
 if (length(bsg) > 0) lpdb$Atom <- structure
 if (length(bsg) == 0) lpdb$Atom <- NULL

 # Hetatm
 # Raw data are at present organized in a data frame
 value <- "HETATM"
 bsg <- which(substr(pdb0,1,6) == value)
 atom_name <- substr(pdb0[bsg],13,16)
 residue_name <- substr(pdb0[bsg],18,20)
 structure <- data.frame(          n = as.integer(substr(pdb0[bsg],7,11)),
                            atomName = I(atom_name),
                              altLoc = I(substr(pdb0[bsg],17,17)),
                         residueName = I(residue_name),
                             chainID = I(substr(pdb0[bsg],22,22)),
                                resN = as.integer(substr(pdb0[bsg],23,26)),
                               iCode = I(substr(pdb0[bsg],27,27)),
                                   x = as.numeric(substr(pdb0[bsg],31,38)),
                                   y = as.numeric(substr(pdb0[bsg],39,46)),
                                   z = as.numeric(substr(pdb0[bsg],47,54)),
                                 occ = as.numeric(substr(pdb0[bsg],55,60)),
                                Bfac = as.numeric(substr(pdb0[bsg],61,66)),
                          atomSymbol = I(substr(pdb0[bsg],77,78)),
                          atomCharge = I(substr(pdb0[bsg],79,80)))
 if (length(bsg) > 0) lpdb$Hetatm <- structure
 if (length(bsg) == 0) lpdb$Hetatm <- NULL

 # Anisotropic values line (if it exists)
 # Raw data are at present organized in a data frame
 value <- "ANISOU"
 bsg <- which(substr(pdb0,1,6) == value)
 atom_name <- substr(pdb0[bsg],13,16)
 residue_name <- substr(pdb0[bsg],18,20)
 structure <- data.frame(          n = as.integer(substr(pdb0[bsg],7,11)),
                            atomName = I(atom_name),
                              altLoc = I(substr(pdb0[bsg],17,17)),
                         residueName = I(residue_name),
                             chainID = I(substr(pdb0[bsg],22,22)),
                                resN = as.integer(substr(pdb0[bsg],23,26)),
                               iCode = I(substr(pdb0[bsg],27,27)),
                                 u11 = as.integer(substr(pdb0[bsg],29,35)),
                                 u22 = as.integer(substr(pdb0[bsg],36,42)),
                                 u33 = as.integer(substr(pdb0[bsg],43,49)),
                                 u12 = as.integer(substr(pdb0[bsg],50,56)),
                                 u13 = as.integer(substr(pdb0[bsg],57,63)),
                                 u23 = as.integer(substr(pdb0[bsg],64,70)),
                          atomSymbol = I(substr(pdb0[bsg],77,78)),
                          atomCharge = I(substr(pdb0[bsg],79,80)))
 if (length(bsg) > 0) lpdb$Anisou <- structure
 if (length(bsg) == 0) lpdb$Anisou <- NULL

 return(lpdb)
}

#
.writePDB <- function(lpdb,filename)
{
 # Remove pdb file, if it does exist, because data will be appended to it
 if (file.exists(filename)) file.remove(filename)

 # Start adding lines to PDB file in correct order

 # *** TITLE SECTION ***
 # HEADER
 idx <- which(names(lpdb) == "header")
 if (length(idx) != 0)
 {
  # Fill description if needed
  if (nchar(lpdb$header[1]) > 40) lpdb$header[1] <- substr(lpdb$header[1],1,40)
  if (nchar(lpdb$header[1]) < 40) for (i in 1:(40-nchar(lpdb$header[1]))) lpdb$header[1] <- paste(lpdb$header[1]," ",sep="")
  if (nchar(lpdb$header[2]) != 9)
  {
   warning("Date part of header appears to be badly formatted. Amending it as feasibly as possible")
   if (nchar(lpdb$header[2]) > 9) lpdb$header[2] <- substr(lpdb$header[2],1,9)
   if (nchar(lpdb$header[2]) < 9) for (i in 1:(9-nchar(lpdb$header[2]))) lpdb$header[2] <- paste(lpdb$header[2]," ",sep="")
  }
  if (nchar(lpdb$header[3]) != 4) stop("Wrong or badly formatted pdb code. Please amend it and start again")
  linea <- paste("HEADER    ",lpdb$header[1],lpdb$header[2],"   ",lpdb$header[3],"              \n",sep="")
  cat(linea,file=filename)
 }

 # OBSLTE (to be completed)
 idx <- which(names(lpdb) == "obslte")
 #if (length(idx) != 0)
 #{
 #}

 # TITLE (to be completed)
 idx <- which(names(lpdb) == "title")
 #if (length(idx) != 0)
 #{
 #}

 # CAVEAT (to be completed)
 idx <- which(names(lpdb) == "caveat")
 #if (length(idx) != 0)
 #{
 #}

 # COMPND (to be completed)
 idx <- which(names(lpdb) == "compnd")
 #if (length(idx) != 0)
 #{
 #}

 # SOURCE (to be completed)
 idx <- which(names(lpdb) == "source")
 #if (length(idx) != 0)
 #{
 #}

 # KEYWDS (to be completed)
 idx <- which(names(lpdb) == "keywds")
 #if (length(idx) != 0)
 #{
 #}

 # EXPDTA (to be completed)
 idx <- which(names(lpdb) == "expdta")
 #if (length(idx) != 0)
 #{
 #}

 # AUTHOR (to be completed)
 idx <- which(names(lpdb) == "author")
 #if (length(idx) != 0)
 #{
 #}

 # REVDAT (to be completed)
 idx <- which(names(lpdb) == "revdat")
 #if (length(idx) != 0)
 #{
 #}

 # SPRSDE (to be completed)
 idx <- which(names(lpdb) == "sprsde")
 #if (length(idx) != 0)
 #{
 #}

 # JRNL (to be completed)
 idx <- which(names(lpdb) == "jrnl")
 #if (length(idx) != 0)
 #{
 #}

 # REMARK (to be completed)
 idx <- which(names(lpdb) == "remark")
 #if (length(idx) != 0)
 #{
 #}

 # CRYST1
 idx <- which(names(lpdb) == "cryst1")
 if (length(idx) != 0)
 {
  ntmp <- 11-nchar(lpdb$cryst1$SG)
  lsgroup <- lpdb$cryst1$SG
  for (i in 1:ntmp) lsgroup <- paste(lsgroup," ",sep="")
  linea <- sprintf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s%4d          \n",lpdb$cryst1[[1]][1],lpdb$cryst1[[1]][2],lpdb$cryst1[[1]][3],
                                                                            lpdb$cryst1[[1]][4],lpdb$cryst1[[1]][5],lpdb$cryst1[[1]][6],lsgroup,lpdb$cryst1$z)
  cat(linea,file=filename,append=TRUE)
 }

 # ORIGX1
 idx <- which(names(lpdb) == "origx1")
 if (length(idx) != 0)
 {
  linea <- sprintf("ORIGX1    %10.6f%10.6f%10.6f     %10.5f                         \n",lpdb$origx1[1],lpdb$origx1[2],lpdb$origx1[3],lpdb$origx1[4])
  cat(linea,file=filename,append=TRUE)
 }

 # ORIGX2
 idx <- which(names(lpdb) == "origx3")
 if (length(idx) != 0)
 {
  linea <- sprintf("ORIGX2    %10.6f%10.6f%10.6f     %10.5f                         \n",lpdb$origx2[1],lpdb$origx2[2],lpdb$origx2[3],lpdb$origx2[4])
  cat(linea,file=filename,append=TRUE)
 }

 # ORIGX3
 idx <- which(names(lpdb) == "origx3")
 if (length(idx) != 0)
 {
  linea <- sprintf("ORIGX3    %10.6f%10.6f%10.6f     %10.5f                         \n",lpdb$origx3[1],lpdb$origx3[2],lpdb$origx3[3],lpdb$origx3[4])
  cat(linea,file=filename,append=TRUE)
 }

 # SCALE1
 idx <- which(names(lpdb) == "scale1")
 if (length(idx) != 0)
 {
  linea <- sprintf("SCALE1    %10.6f%10.6f%10.6f     %10.5f                         \n",lpdb$scale1[1],lpdb$scale1[2],lpdb$scale1[3],lpdb$scale1[4])
  cat(linea,file=filename,append=TRUE)
 }

 # SCALE2
 idx <- which(names(lpdb) == "scale2")
 if (length(idx) != 0)
 {
  linea <- sprintf("SCALE2    %10.6f%10.6f%10.6f     %10.5f                         \n",lpdb$scale2[1],lpdb$scale2[2],lpdb$scale2[3],lpdb$scale2[4])
  cat(linea,file=filename,append=TRUE)
 }

 # SCALE3
 idx <- which(names(lpdb) == "scale3")
 if (length(idx) != 0)
 {
  linea <- sprintf("SCALE3    %10.6f%10.6f%10.6f     %10.5f                         \n",lpdb$scale3[1],lpdb$scale3[2],lpdb$scale3[3],lpdb$scale3[4])
  cat(linea,file=filename,append=TRUE)
 }

 # For ATOM, ANISOU and HETATM (a few if conditions needed here as ANISOU is interspersed between ATOM and HETATM lines)

 # Indices to check what is available
 idxato <- which(names(lpdb) == "Atom")
 idxhet <- which(names(lpdb) == "Hetatm")
 idxani <- which(names(lpdb) == "Anisou")
 
 # ATOM only or ATOM and HETATM only
 if (length(idxato) != 0 & length(idxani) == 0)
 {
  stc <- lpdb$Atom
  natoms <- length(stc$n)
  tmp <- as.factor(stc$chainID)
  nter <- 0
  for (i in 1:natoms)
  {
   linea <- sprintf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %s%s\n",i+nter,stc$atomName[i],stc$altLoc[i],stc$residueName[i],
                                                                                           stc$chainID[i],stc$resN[i],stc$iCode[i],stc$x[i],stc$y[i],stc$z[i],
                                                                                           stc$occ[i],stc$Bfac[i],stc$atomSymbol[i],stc$atomCharge[i])
   cat(linea,file=filename,append=TRUE)
   
   # Add TER if needed
   if (i != natoms)
   {
    lev1 <- tmp[i]
    lev2 <- tmp[i+1]
    if (lev1 != lev2)
    {
     nter <- nter+1
     linea <- sprintf("TER   %5d      %s %s%4d%s                                                     \n",i+nter,stc$residueName[i],stc$chainID[i],stc$resN[i],stc$iCode[i])
     cat(linea,file=filename,append=TRUE)
    }
   }
   if (i == natoms)
   {
     nter <- nter+1
     linea <- sprintf("TER   %5d      %s %s%4d%s                                                     \n",i+nter,stc$residueName[i],stc$chainID[i],stc$resN[i],stc$iCode[i])
     cat(linea,file=filename,append=TRUE)
   }
  }

  # HETATM
  if (length(idxhet) > 0)
  {
   stc <- lpdb$Hetatm
   nadd <- natoms+nter
   natoms <- length(stc$n)
   for (i in 1:natoms)
   {
    linea <- sprintf("HETATM%5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %s%s\n",i+nadd,stc$atomName[i],stc$altLoc[i],stc$residueName[i],
                                                                                            stc$chainID[i],stc$resN[i],stc$iCode[i],stc$x[i],stc$y[i],stc$z[i],
                                                                                            stc$occ[i],stc$Bfac[i],stc$atomSymbol[i],stc$atomCharge[i])
    cat(linea,file=filename,append=TRUE)
   }
  }
 }

 # ANISOU
 if (length(idxato) == 1 & length(idxani) == 1)
 {
  stcato <- lpdb$Atom
  if (length(idxhet) > 0) stchet <- lpdb$Hetatm
  stcani <- lpdb$Anisou
  natom <- length(stcato$n)
  if (length(idxhet) > 0) nhetatm <- length(stchet$n)
  if (length(idxhet) == 0) nhetatm <- 0
  nter <- 0
  tmp <- as.factor(lpdb$Atom$chainID)

  # Alternate lines: ATOM - ANISOU or HETATM - ANISOU (sometimes TER is at the end of chains)
  for (i in 1:natom)
  {
   j <- nter+i
   linea <- sprintf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %s%s\n",j,stcato$atomName[i],stcato$altLoc[i],stcato$residueName[i],
                                                                                           stcato$chainID[i],stcato$resN[i],stcato$iCode[i],stcato$x[i],stcato$y[i],stcato$z[i],
                                                                                           stcato$occ[i],stcato$Bfac[i],stcato$atomSymbol[i],stcato$atomCharge[i])
   cat(linea,file=filename,append=TRUE)
   linea <- sprintf("ANISOU%5d %s%s%s %s%4d%s %7d%7d%7d%7d%7d%7d      %s%s\n",j,stcani$atomName[i],stcani$altLoc[i],stcani$residueName[i],
                                                                              stcani$chainID[i],stcani$resN[i],stcani$iCode[i],stcani$u11[i],stcani$u22[i],stcani$u33[i],
                                                                              stcani$u12[i],stcani$u13[i],stcani$u23[i],stcani$atomSymbol[i],stcani$atomCharge[i])
   cat(linea,file=filename,append=TRUE)
   
   # Add TER if needed
   if (i != natom)
   {
    lev1 <- tmp[i]
    lev2 <- tmp[i+1]
    if (lev1 != lev2)
    {
     nter <- nter+1
     linea <- sprintf("TER   %5d      %s %s%4d%s                                                     \n",i+nter,stcato$residueName[i],stcato$chainID[i],
                                                                                                         stcato$resN[i],stcato$iCode[i])
     cat(linea,file=filename,append=TRUE)
    }
   }
   if (i == natom)
   {
     nter <- nter+1
     linea <- sprintf("TER   %5d      %s %s%4d%s                                                     \n",i+nter,stcato$residueName[i],stcato$chainID[i],
                                                                                                         stcato$resN[i],stcato$iCode[i])
     cat(linea,file=filename,append=TRUE)
   }
  }
  for (i in 1:nhetatm)
  {
   j <- natom+nter+i
   jani <- natom+i
   linea <- sprintf("HETATM%5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %s%s\n",j,stchet$atomName[i],stchet$altLoc[i],stchet$residueName[i],
                                                                                           stchet$chainID[i],stchet$resN[i],stchet$iCode[i],stchet$x[i],stchet$y[i],stchet$z[i],
                                                                                           stchet$occ[i],stchet$Bfac[i],stchet$atomSymbol[i],stchet$atomCharge[i])
   cat(linea,file=filename,append=TRUE)
   linea <- sprintf("ANISOU%5d %s%s%s %s%4d%s %7d%7d%7d%7d%7d%7d      %s%s\n",j,stcani$atomName[jani],stcani$altLoc[jani],stcani$residueName[jani],
                                                                              stcani$chainID[jani],stcani$resN[jani],stcani$iCode[jani],
                                                                              stcani$u11[jani],stcani$u22[jani],stcani$u33[jani],
                                                                              stcani$u12[jani],stcani$u13[jani],stcani$u23[jani],stcani$atomSymbol[jani],
                                                                              stcani$atomCharge[jani])
   cat(linea,file=filename,append=TRUE)
  }
 }
}
