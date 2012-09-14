# R code to implement the crystallographic ideas related to reflections.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## Classes
#
# MergedReflection
setClass(
         Class="MergedReflection",
         representation=representation(rcell="ReciprocalUnitCell",symmetry="Symmetry",records="data.frame",dtypes="character")
        )

## Default methods
#
# Print
#setMethod(
#          f="print",
#          signature="Lattice",
#          definition=function(x,...)
#                     {
#                     }
#         )

## Generic methods
#
# Check MergedReflection object is correct
setGeneric(
           name="checkMergedReflection",
           def=function(object){standardGeneric("checkMergedReflection")}
          )
#
setMethod(
          f="checkMergedReflection",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Check records slot has more than 3 columns
                      if (length(object@records) <= 3) stop("Slot records in MergedReflection object needs to have more than 3 columns")

                      # Check first three columns of records slot data.frame are H,K,L
                      tmp <- names(object@records[,1:3])
                      if (tmp[1] != "H" | tmp[2] != "K" | tmp[3] != "L") stop("First three columns of records slot in MergedReflection object have to be 'H', 'K' and 'L'")

                      # Check first three data types are H H H
                      tmp <- object@dtypes
                      if (tmp[1] != "H" | tmp[2] != "H" | tmp[3] != "H") stop("First three columns of records slot in MrgedReflection object have to be Miller indices, of data types 'H', 'H', 'H'")

                      # Check number of data types is equal to number of columns in records slot
                      if (length(object@records) != length(object@dtypes)) stop("Number of data types does not coincide with number of records columns in MergedReflectio object")

                      # Check compatibility of space group and cell parameters
                      cell <- computeUnitCell(object@rcell)
                      tmplatt <- checkSymmetryWithCell(object@symmetry,cell)

                      # Check systematic absences
                      sym_number <- object@symmetry@sym_number
                      setting <- object@symmetry@setting
                      hkl <- object@records[,1:3]
                      fans <- .check_systematic_absences(hkl,sym_number,setting)
                      if (!fans) stop("MergedReflection object includes systematic absences")

                      # Everything OK. Green light.
                     }
         )
#
# Extract ReciprocalCell object
setGeneric(
           name="extractReciprocalCell",
           def=function(object){standardGeneric("extractReciprocalCell")}
          )
#
setMethod(
          f="extractReciprocalCell",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Check object is a MergedReflection object

                      # Extract rcell slot
                      tmp <- object@rcell

                      return(tmp)
                     }
         )
#
# Extract UnitCell object
setGeneric(
           name="extractUnitCell",
           def=function(object){standardGeneric("extractUnitCell")}
          )
#
setMethod(
          f="extractUnitCell",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Check object is a MergedReflection object

                      # Extract rcell slot
                      tmp <- object@rcell

                      # Compute cell from rcell
                      tmp <- computeUnitCell(tmp)

                      return(tmp)
                     }
         )
#
# Extract Symmetry object
setGeneric(
           name="extractSymmetry",
           def=function(object){standardGeneric("extractSymmetry")}
          )
#
setMethod(
          f="extractSymmetry",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Check object is a MergedReflection object

                      # Extract symmetry slot
                      tmp <- object@symmetry

                      return(tmp)
                     }
         )
#
# Extract reflections records (a data frame)
setGeneric(
           name="extractData",
           def=function(object){standardGeneric("extractData")}
          )
#
setMethod(
          f="extractData",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Check object is a MergedReflection object

                      # Extract records slot 
                      tmp <- object@records

                      return(tmp)
                     }
         )
#
# Compute resolutions
setGeneric(
           name="computeResolutions",
           def=function(object){standardGeneric("computeResolutions")}
          )
#
setMethod(
          f="computeResolutions",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Check object is a MergedReflection object

                      # Extract records slot 
                      #tmp <- na.omit(object@records)
                      tmp <- object@records

                      # Extract unit cell
                      cell <- extractUnitCell(object)

                      # Compute resolutions for all reflections
                      resos <- .d_hkl(tmp$H,tmp$K,tmp$L,cell@a,cell@b,cell@c,cell@alpha@ang,cell@beta@ang,cell@gamma@ang)

                      # Vector containing min and max resolution
                      #resos <- sort(range(resos),decreasing=TRUE)

                      return(resos)
                     }
         )
#
# Expand to full reciprocal space
setGeneric(
           name="expandMergedReflection",
           def=function(object){standardGeneric("expandMergedReflection")}
          )
#
setMethod(
          f="expandMergedReflection",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Extract Symmetry object from MergedReflections object
                      symmetry <- extractSymmetry(object)  
                    
                      # Point group operators
                      lpg <- extractSymmetryOperators(symmetry)$PG

                      # Extract actual data
                      refs <- extractData(object)

                      # Expansion
                      ncol <- length(refs)
                      newrefs <- refs
                      for (i in 1:length(lpg))
                      {
                       # Non-Friedel
                       if (i != 1)
                       {
                        tmp <- t(mapply((function(h,k,l,m) {r <- c(h,k,l)%*%m}),refs$H,refs$K,refs$L,MoreArgs=list(m=lpg[[i]])))
                        tmp2 <- data.frame(H=tmp[,1],K=tmp[,2],L=tmp[,3])
                        tmp2 <- cbind(tmp2,refs[,4:ncol])
                        colnames(tmp2) <- colnames(refs)
                        newrefs <- rbind(newrefs,tmp2)
                       }

                       # Friedel
                       tmp <- t(mapply((function(h,k,l,m) {r <- c(h,k,l)%*%m}),refs$H,refs$K,refs$L,MoreArgs=list(m=-lpg[[i]])))
                       tmp2 <- data.frame(H=tmp[,1],K=tmp[,2],L=tmp[,3])
                       tmp2 <- cbind(tmp2,refs[,4:ncol])
                       colnames(tmp2) <- colnames(refs)
                       newrefs <- rbind(newrefs,tmp2)
                      }

                      # Get rid of duplicates
                      andare <- duplicated(newrefs[,1:3])
                      idx <- which(!andare)
                      newrefs <- newrefs[idx,]

                      # Copy expanded reflections into records slot of MergedReflection object
                      object@records <- newrefs

                      return(object)
                     }
         )
#
# Compute cartesian coordinates for reciprocal nodes in a MergedReflection object
setGeneric(
           name="computeCartesianCoordinates",
           def=function(object,...){standardGeneric("computeCartesianCoordinates")}
          )
#
setMethod(
          f="computeCartesianCoordinates",
          signature="MergedReflection",
          definition=function(object,convention=2)
                     {
                      # Extract reflection data
                      refs <- extractData(object)

                      # Extract UnitCell object
                      cella <- extractUnitCell(object)

                      # Extract cell parameters
                      cpar <- extractCellParameters(cella)

                      # Crystal axes
                      if (convention == 1) m <- .triclinic_to_orthogonal_01(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])
                      if (convention == 2) m <- .triclinic_to_orthogonal_02(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])
                      av <- m[1,]
                      bv <- m[2,]
                      cv <- m[3,]
                      V <- computeCellVolume(cella)

                      # Reciprocal axes
                      avr <- bv%X%cv/V
                      bvr <- cv%X%av/V
                      cvr <- av%X%bv/V
                      mr <- matrix(c(avr,bvr,cvr),nrow=3,ncol=3,byrow=T)

                      # calculate cartesian coordinates for all reciprocal points
                      xyz <- t(mapply((function(h,k,l,m) {r <- c(h,k,l)%*%m}),refs$H,refs$K,refs$L,MoreArgs=list(m=mr)))

                      # Append coordinates to original data frame
                      refs <- cbind(refs,data.frame(X=xyz[,1],Y=xyz[,2],Z=xyz[,3]))

                      # Copy enriched reflections into records slot of MergedReflection object
                      object@records <- refs

                      return(object)
                     }
         )
#
# Compute spherical polar coordinates for reciprocal nodes in a MergedReflection object
setGeneric(
           name="computeSphericalCoordinates",
           def=function(object,...){standardGeneric("computeSphericalCoordinates")}
          )
#
setMethod(
          f="computeSphericalCoordinates",
          signature="MergedReflection",
          definition=function(object,convention=2)
                     {
                      # Extract reflection data
                      refs <- extractData(object)

                      # Extract UnitCell object
                      cell <- extractUnitCell(object)

                      # Extract cell parameters
                      cpar <- extractCellParameters(cella)

                      # Crystal axes
                      if (convention == 1) m <- .triclinic_to_orthogonal_01(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])
                      if (convention == 2) m <- .triclinic_to_orthogonal_02(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])
                      av <- m[1,]
                      bv <- m[2,]
                      cv <- m[3,]
                      V <- computeCellVolume(cella)

                      # Reciprocal axes
                      avr <- bv%X%cv/V
                      bvr <- cv%X%av/V
                      cvr <- av%X%bv/V
                      mr <- matrix(c(avr,bvr,cvr),nrow=3,ncol=3,byrow=T)

                      # calculate cartesian coordinates for all reciprocal points
                      xyz <- t(mapply((function(h,k,l,m) {r <- c(h,k,l)%*%m}),refs$H,refs$K,refs$L,MoreArgs=list(m=mr)))

                      # Calculate spherical-polar coordinates
                      r <- sqrt(xyz[,1]^2+xyz[,2]^2+xyz[,3]^2)
                      t <- acos(xyz[,3]/r)
                      p <- atan2(xyz[,2]/(r*sin(t)),xyz[,1]/(r*sin(t)))

                      # For all theta == 0 p is undefined; let us fix it to 0
                      idx <- which(t > -0.000001 & t < 0.000001)
                      p[idx] <- 0

                      # Also all p have to be between -pi (excluded) and +pi 
                      idx <- which(p < -pi)
                      p[idx] <- p[idx]+pi
                      idx <- which(p > pi)
                      p[idx] <- p[idx]-pi
                      idx <- which(p < -pi+0.00001)
                      p[idx] <- pi

                      # Append coordinates to original data frame
                      refs <- cbind(refs,data.frame(R=r,THETA=t,PHI=p))

                      # Copy enriched reflections into records slot of MergedReflection object
                      object@records <- refs

                      return(object)
                     }
         )
#
# Compute completeness
setGeneric(
           name="computeCompleteness",
           def=function(object){standardGeneric("computeCompleteness")}
          )
#
setMethod(
          f="computeCompleteness",
          signature="MergedReflection",
          definition=function(object)
                     {
                      # Compute max resolution (this will also check object is of class MergedReflection)
                      max_reso <- range(computeResolutions(object))[2]

                      # Extract unit cell
                      cell <- extractUnitCell(object)

                      # Extract symmetry
                      symmetry <- extractSymmetry(object)

                      # Compute MergedReflection object from scratch to generate all reflections up to max_reso resolution
                      tmp <- createMergedReflection(cell,symmetry,max(sort(c(cell@a,cell@b,cell@c)))+100,max_reso)

                      # Extract data for original object and for tmp
                      d1 <- extractData(object)
                      d2 <- extractData(tmp)

                      # Get rid of NA's in original object
                      d1 <- na.omit(d1)

                      ctness <- (length(d1$H)/length(d2$H))*100

                      # Delete unwanted objects
                      rm(tmp,d1,d2)

                      return(ctness)
                     }
         )
#
# Write reflection object into a file with MTZ format
setGeneric(
           name="exportMergedReflectionToMTZ",
           def=function(object,...){standardGeneric("exportMergedReflectionToMTZ")}
          )
#
setMethod(
          f="exportMergedReflectionToMTZ",
          signature="MergedReflection",
          definition=function(object,filename="export.mtz",title="Untitled",sort=c(1,2,3,0,0),
                              records_id=NA,project=NA,crystal=NA,dataset=NA,dcell=NA,dwavel=NA,history=NA)
                     {
                      # Check filename is a character string
                      if (!is.character(filename)) stop("MTZ file name needs to be a character string")

                      # Check title is a character string
                      if (!is.character(title)) stop("Title needs to be a character string")

                      # Check sort is a numeric vector
                      if (!is.numeric(sort)) stop("Sort needs to be a vector of integer numbers")
 
                      # First part of mtz are data
                      dt1 <- extractData(object)  # This method checks also for object correctness

                      # Compute various quantities to build header
                      ncol <- c(as.integer(length(dt1)),as.integer(length(dt1[,1])),as.integer(0))    # Number of columns, number of reflections 
                                                                                                      # and number of batch header; 0 because this is a merged file
                      cell <- extractCellParameters(extractUnitCell(object)) # Cell parameters
                      idx <- c()                                                                                                 #
                      idx <- c(idx,which(sort == 1))                                                                             #
                      idx <- c(idx,which(sort == 2))                                                                             # To sort data according to user choice
                      idx <- c(idx,which(sort == 3))                                                                             #
                      idx <- c(idx,which(sort == 4))                                                                             #
                      idx <- c(idx,which(sort == 5))                                                                             #
                      sort <- as.integer(sort)
                      if (length(idx) == 5) dt1 <- dt1[order(dt1[,idx[1]],dt1[idx[2]],dt1[,idx[3]],dt1[,idx[4]],dt1[,idx[5]]),]  #
                      if (length(idx) == 4) dt1 <- dt1[order(dt1[,idx[1]],dt1[idx[2]],dt1[,idx[3]],dt1[,idx[4]]),]               #
                      if (length(idx) == 3) dt1 <- dt1[order(dt1[,idx[1]],dt1[idx[2]],dt1[,idx[3]]),]                            #
                      if (length(idx) == 2) dt1 <- dt1[order(dt1[,idx[1]],dt1[idx[2]]),]                                         #
                      if (length(idx) == 1) dt1 <- dt1[order(dt1[,idx[1]]),]                                                     #
                      sym_old <- paste("'",.translate_SG(object@symmetry@sym_xHM,SG_in="xHM",SG_out="old"),"'",sep="")
                      stmp <- .extract_symmetry_info(object@symmetry@sym_xHM)$PGRP                                               # Point group info from syminfo.lib file
                      first <- 1
                      last <- nchar(stmp)
                      while (substr(stmp,last-first,last-first) != "'") first <- first+1
                      pgrp <- substr(stmp,last-first,last)
                      syminf <- list(length(object@symmetry@PG)*length(object@symmetry@CT),                               # Number of symmetry operators
                                     length(object@symmetry@PG),                                                          # Number of primitive symmetry ops
                                     substring(.translate_SG(object@symmetry@sym_xHM,SG_in="xHM",SG_out="old"),1,1),      # Lattice type
                                     object@symmetry@sym_number,                                                          # Space group number
                                     sym_old,                                                                             # Space group name (old)
                                     pgrp)                                                                                # Point group name
                      symm <- .full_symm_strings(object@symmetry@sym_xHM)                                                 #
                      reso <- 1/range(computeResolutions(object))^2                                                              # Resolution limits (1/d^2)
                      labels <- names(object@records)                                                                     #
                      types <- object@dtypes                                                                              #
                      rmin <- sprintf("%18.4f",apply(object@records,2,min,na.rm=TRUE))
                      rmax <- sprintf("%18.4f",apply(object@records,2,max,na.rm=TRUE))
                      if (is.na(records_id[1]))                                                                           # Columns: names, types, range and id
                      {                                                                                                   #
                       rid <- rep(1,times=length(object@records))                                                         #
                       rid[1:3] <- c(0,0,0)                                                                               #
                      }                                                                                                   #
                      if (!is.na(records_id[1]))                                                                          #
                      {                                                                                                   #
                       if (length(records_id) != length(object@records)) stop("Argument records_id is not of correct length")
                       rid <- records_id                                                                                  #
                      }                                                                                                   #
                      rid <- as.character(rid)
                      column <- data.frame(labels=I(labels),types=I(types),min=I(rmin),max=I(rmax),id=I(rid))                      #
                      rownames(column) <- 1:length(labels)
                      ndif <- length(unique(rid))                                                                         # Number of datasets
                      if (!is.data.frame(project))                                                                        #
                      {                                                                                                   #
                       rid <- c(0,1)                                                                                      #
                       pname <- c("HKL_base","Unnamed project")                                                           # Project id and name
                       tproj <- data.frame(id=rid,pname=I(pname))                                                            #
                      }                                                                                                   #
                      if (is.data.frame(project))                                                                         #
                      {                                                                                                   #
                       if (length(project[,1]) != ndif) stop("Number of datasets incompatible with project assignement")  #
                       tproj <- project                                                                                   #
                      }                                                                                                   #
                      if (!is.data.frame(crystal))                                                                        #
                      {                                                                                                   #
                       if (length(unique(column$id)) != 2) stop("Number of datasets incompatible with crystal assignement")
                       rid <- c(0,1)                                                                                      
                       cname <- c("HKL_base","Unnamed crystal") 
                       tcryst <- data.frame(id=rid,cname=I(cname))                                                           #
                      }                                                                                                   #
                      if (is.data.frame(crystal))                                                                         #
                      {                                                                                                   #
                       if (length(crystal[,1]) != ndif) stop("Number of datasets incompatible with crystal assignement")  #
                       tcryst <- crystal                                                                                  #
                      }                                                                                                   #
                      if (!is.data.frame(dataset))                                                                        #
                      {                                                                                                   #
                       if (length(unique(column$id)) != 2) stop("Number of datasets incompatible with dataset assignement")
                       rid <- c(0,1)                                                                                      #
                       dname <- c("HKL_base","Unnamed dataset")                                                           # Dataset id and name
                       tdata <- data.frame(id=rid,dname=I(dname))                                                            #
                      }                                                                                                   #
                      if (is.data.frame(dataset))                                                                         #
                      {                                                                                                   #
                       if (length(dataset[,1]) != ndif) stop("Number of datasets incompatible with dataset assignement")  #
                       tdata <- dataset                                                                                   #
                      }                                                                                                   #
                      if (!is.data.frame(dcell))                                                                          #
                      {                                                                                                   #
                       if (length(unique(column$id)) != 2) stop("Number of datasets incompatible with dcell assignement") #
                       rid <- c(0,1)                                                                                      #
                       a <- c(cell[1],cell[1])                                                                            # Dcell id and parameters
                       b <- c(cell[2],cell[2])                                                                            #
                       c <- c(cell[3],cell[3])                                                                            #
                       alpha <- c(cell[4],cell[4])                                                                        #
                       beta <- c(cell[5],cell[5])                                                                         #
                       gamma <- c(cell[6],cell[6])                                                                        #
                       tcell <- data.frame(id=rid,a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)                          #
                      }                                                                                                   #
                      if (is.data.frame(dcell))                                                                           #
                      {                                                                                                   #
                       if (length(dcell[,1]) != ndif) stop("Number of datasets incompatible with dcell assignement")      #
                       for (i in 1:ndif)
                       {
                        if (abs(dcell[i,2]-cell[1]) > 1 | abs(dcell[i,3]-cell[2]) > 1 | abs(dcell[i,4]-cell[3]) > 1 |
                            abs(dcell[i,5]-cell[4]) > 1 | abs(dcell[i,6]-cell[5]) > 1 | abs(dcell[i,7]-cell[6]) > 1)
                        {
                         warning(paste("Cell parameters in dataset",i,"in DCELL dataframe are different from those of reflection object"))
                        }
                       }
                       tcell <- dcell                                                                                     #
                      }                                                                                                   #
                      if (!is.data.frame(dwavel))                                                                          #
                      {                                                                                                    #
                       if (length(unique(column$id)) != 2) stop("Number of datasets incompatible with dwavel assignement") #
                       rid <- c(0,1)                                                                                       #
                       lambda <- c(0,0)                                                                                    # Dwavel id and values
                       twavel <- data.frame(id=rid,lambda=lambda)                                                          #
                      }                                                                                                    #
                      if (is.data.frame(dwavel))                                                                           #
                      {                                                                                                    #
                       if (length(dwavel[,1]) != ndif) stop("Number of datasets incompatible with dwavel assignement")     #
                       twavel <- dwavel                                                                                    #
                      }                                                                                                    #
                      if (is.na(history[1]))
                      {
                       thistory <- "NO HISTORY FOR THIS MTZ FILE"
                       while (nchar(thistory) != 80) thistory <- paste(thistory," ",sep="") 
                      }
                      if (!is.na(history[1]))
                      {
                       if (!is.character(history)) stop("history field in this function needs to be a character vector")
                       thistory <- history
                       for (i in 1:length(thistory))
                       {
                        if (nchar(thistory[i]) != 80) while (nchar(thistory[i]) != 80) thistory[i] <- paste(thistory[i]," ",sep="")
                       }
                      }
                      ndif <- as.numeric(ndif)       # Just for being able to write mtz binary format

                      # Named list to return
                      dt2 <- list(NCOL=ncol,CELL=cell,SORT=sort,SYMINF=syminf,RESO=reso,NDIF=ndif,SYMM=symm,
                                  PROJECT=tproj,CRYSTAL=tcryst,DATASET=tdata,DCELL=tcell,DWAVEL=twavel,COLUMN=column,HISTORY=thistory)
                      dt3 <- NULL

                      # Assemble all parts in a named list
                      dta <- list(reflections=dt1,header=dt2,batch_header=dt3)

                      # Write out to an mtz file
                      .writeMTZ(dta,file=filename,title=title)

                      #return(dta)
                     }
         )
#
# Modify slots of MergedReflection object
setGeneric(
           name="changeMergedReflection",
           def=function(object,...){standardGeneric("changeMergedReflection")}
          )
#
setMethod(
          f="changeMergedReflection",
          signature="MergedReflection",
          definition=function(object,rcell=NA,symmetry=NA,records=NA,dtypes=NA)
                     {
                      options(warn=-1)  # Suppress warnings
                      if (!is.na(rcell)) object@rcell <- rcell
                      if (!is.na(symmetry)) object@symmetry <- symmetry
                      if (!is.na(records)) object@records <- records
                      if (!is.na(dtypes)) object@dtypes <- dtypes
                      options(warn=0)   # Re-establish warnings

                      # Check new MergedReflection object is OK
                      #checkMergedReflection(object)

                      return(object)
                     }
         )


## Friendly constructors
#
# Create a MergedReflection object with reflections slot having the 5 columns H K L FP SIGFP, and both
# FP and SIGFP set to 0. All reflections between given resolutions will be produced.
createMergedReflection <- function(cell,symmetry,reso_min,reso_max)
{
 # Check cell is a proper UnitCell object
 checkUnitCell(cell)

 # Compute reciprocal unit cell
 rcell <- computeReciprocalUnitCell(cell)

 # Check symmetry is a proper Symmetry object
 checkSymmetry(symmetry)

 # Generate all Miller indices up to reso_max resolution. tmp_hkl has the following columns: h k l s
 tmp_hkl <- .generate_reflections_up_to_resolution(cell@a,cell@b,cell@c,cell@alpha@ang,cell@beta@ang,cell@gamma@ang,reso_max)

 # Get rid of reflections with resolution lower than reso_min
 tmp_hkl <- tmp_hkl[tmp_hkl$s > (1/reso_min),]

 # Reduce to hkl asymmetric unit (here we also get rid of systematic absences)
 tmp_hkl <- .apply_hklasu(symmetry@sym_number,tmp_hkl,symmetry@setting)
 
 # Final data frame: H K L FP SIGFP
 hkl <- data.frame(H=tmp_hkl$h,K=tmp_hkl$k,L=tmp_hkl$l,FP=rep(0,times=length(tmp_hkl[,1])),SIGFP=rep(0,times=length(tmp_hkl[,1])))

 # Data types
 dtypes <- c("H","H","H","F","Q")

 # Assemble MergedReflection object
 mr <- new(Class="MergedReflection",rcell=rcell,symmetry=symmetry,records=hkl,dtypes=dtypes)

 # Check object created is a valid MergedReflection object
 checkMergedReflection(mr)

 return(mr)
}
#
# Create a MergedReflection object starting from an mtz file
createMergedReflectionFromMTZ <- function(file)
{
 # Load mtz file in named list
 lmtz <- .readMTZ(file,messages=FALSE)

 # If batch_header is not NULL can't proceed, as this is an unmerged file
 if (!is.null(lmtz$batch_header[1])) stop("This mtz file contains unmerged data")

 # Extract cell and create UnitCell object
 cpar <- lmtz$header$CELL
 rcell <- computeReciprocalUnitCell(createUnitCell(a=cpar[1],b=cpar[2],c=cpar[3],alpha=cpar[4],beta=cpar[5],gamma=cpar[6]))
 checkReciprocalUnitCell(rcell)

 # Extract symmetry number and create Symmetry object
 sym_number <- lmtz$header$SYMINF[[4]]
 SG <- createSymmetry(sym_number)
 checkSymmetry(SG)

 # Check columns H, K, L are included in lmtz$reflections
 refcolnames <- colnames(lmtz$reflections)
 if (sum("H" == refcolnames) != 1) stop("mtz file contains no or more than one 'H' columns")
 if (sum("K" == refcolnames) != 1) stop("mtz file contains no or more than one 'K' columns")
 if (sum("L" == refcolnames) != 1) stop("mtz file contains no or more than one 'L' columns")

 # Extract columns data types
 dtypes <- as.character(lmtz$header$COLUMN$types)

 # Create MergedReflection object
 object <- new(Class="MergedReflection",rcell=rcell,symmetry=SG,records=lmtz$reflections,dtypes=dtypes)

 # Check object created is a valid MergedReflection object
 #checkMergedReflection(object)

 return(object)
}
