## Module containing stuff for crystallography with CCP4 
## Suggested names:  cRy -   


# Generic functions ### GENERIC ### ### GENERIC ### ### GENERIC ### ### GENERIC ### ### GENERIC ### ### GENERIC ### START
.tagRead <- function(stringa)
{
 # Given a string whose words are separated by blanks, returns
 # a list whose first element is all that matter. This is a vector
 # containing the individual words composing the string.
 blanks <- "[[:blank:]]+"  
 ltmp <- strsplit(stringa,blanks)

 return(ltmp[[1]])
}

# Generic functions ### GENERIC ### ### GENERIC ### ### GENERIC ### ### GENERIC ### ### GENERIC ### ### GENERIC ###   END


# Functions to read MTZ ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### START
.readIR <- function(f)
{
 # Reads identification record of an MTZ binary file.
 # 1) "MTZ "
 # 2) Machine stamp

 # Hex version of ascii "MTZ "
 rawBlank <- charToRaw(" ")
 rawMTZ <- c(as.raw(77),as.raw(84),as.raw(90),rawBlank)

 # errF is zero if identification record is OK; otherwise is 1
 errF <- 0

 # Check initial word is "MTZ "
 mtzStamp <- readBin(f,"raw",n=4)
 for (i in 1:4)
 {
  if (mtzStamp[i] != rawMTZ[i])
  {
   errF <- 1
   return(list(errF,NULL))
  }
 }
 
 # Location of start for data records. If subsequently-read machineStamp is not "DA", then use big endian
 edn <- 1
 headerLoc <- readBin(f,"integer",n=1,size=4)

 # Machine stamp
 machineStamp <- readChar(f,4) #number formats of the architecture f was written on

 # If machineStamp != "DA" re-read headerLoc with big endian and change edn into -1
 if (machineStamp != "DA")
 {
  seek(f,4)
  headerLoc <- readBin(f,"integer",n=1,size=4,endian="big")
  edn <- -1
 }

 # Prepare for data records reading: pPosition connection at beginning of reflection data (21st byte=4*20+1)
 seek(f,80)

 return(list(errF,headerLoc,edn))
}
 
.readH <- function(f,messages)
{
 # Given a connection f corresponding to the binary MTZ file immediately
 # after end of reflection records, this function keeps reading the file to extract
 # information included in the header.

 col_symm <- c()
 col_labels <- c()
 col_types <- c()
 col_min <- c()
 col_max <- c()
 col_id <- c()
 d_project <- data.frame(id=1,pname=NA)
 d_crystal <- data.frame(id=1,cname=NA)
 d_dataset <- data.frame(id=1,dname=NA)
 d_dcell <- data.frame(id=1,a=NA,b=NA,c=NA,alpha=NA,beta=NA,gamma=NA)
 d_dwavel <- data.frame(id=1,lambda=NA)
 col_batch <- c()
 l_data <- list()
 endTag=""
 while (endTag != "END")
 {
  hdata <- readChar(f,80)
  #print(hdata)
  fields <- .tagRead(hdata)
  endTag <- fields[1]

  # Distribute values according to keyword
  if (endTag == "NCOL") l_data$NCOL <- c(as.integer(fields[2]),           # Number of columns
                                         as.integer(fields[3]),           # Number of reflection records
                                         as.integer(fields[4]))           # Number of batches
  if (endTag == "CELL") l_data$CELL <- c(as.numeric(fields[2]),           # Cell parameter a
                                         as.numeric(fields[3]),           # Cell parameter b
                                         as.numeric(fields[4]),           # Cell parameter c
                                         as.numeric(fields[5]),           # Cell parameter alpha
                                         as.numeric(fields[6]),           # Cell parameter beta
                                         as.numeric(fields[7]))           # Cell parameter gamma
  if (endTag == "SORT") l_data$SORT <- c(as.integer(fields[2]),           # Sort order 1
                                         as.integer(fields[3]),           # Sort order 2
                                         as.integer(fields[4]),           # Sort order 3
                                         as.integer(fields[5]),           # Sort order 4
                                         as.integer(fields[6]))           # Sort order 5
  if (endTag == "SYMINF")
  {
   # Need to treat this case differently from the others
   cplace <- c()
   for (i in 1:nchar(hdata))
   {
    if (substr(hdata,i,i) == "'") cplace <- c(cplace,i)
   }
   if (length(cplace) != 0 & length(cplace) != 2 & length(cplace) != 4) stop("Wrongly formatted SYMINF line in header")
   if (length(cplace) == 0)
   {
    sgname <- fields[length(fields)-1]
    pgname <- fields[length(fields)]
   }
   if (length(cplace) == 2)
   {
    sgname <- substr(hdata,cplace[1],cplace[2])
    pgname <- fields[length(fields)]
   }
   if (length(cplace) == 4)
   {
    sgname <- substr(hdata,cplace[1],cplace[2])
    pgname <- substr(hdata,cplace[3],cplace[4])
   }
   #sgname <- paste("'",strsplit(hdata,"'",fixed=TRUE)[[1]][2],"'",sep="")
   #pgname <- paste("'",strsplit(hdata,"'",fixed=TRUE)[[1]][4],"'",sep="")
   l_data$SYMINF <- list(as.integer(fields[2]),    # Number of symmetry operations
                         as.integer(fields[3]),    # Number of primitive symmetry operations
                                    fields[4] ,    # Lattice type
                         as.integer(fields[5]),    # Space group number
                                       sgname ,    # Space group name
                                    #fields[7])     # Point group name
                                       pgname)     # Point group name
  }
  #if (endTag == "SYMINF") l_data$SYMINF <- list(as.integer(fields[2]),    # Number of symmetry operations
  #                                              as.integer(fields[3]),    # Number of primitive symmetry operations
  #                                              fields[4],                # Lattice type
  #                                              as.integer(fields[5]),    # Space group number
  #                                              fields[6],                # Space group name
  #                                              fields[7])                # Point group name
  if (endTag == "SYMM") col_symm <- c(col_symm,hdata)
  if (endTag == "RESO") l_data$RESO <- c(as.numeric(fields[2]),           # Minimum resolution (stored as 1/d-squared)
                                         as.numeric(fields[3]))           # Maximum resolution (stored as 1/d-squared)
  if (endTag == "NDIF") l_data$NDIF <- c(as.numeric(fields[2]))           # Number of datasets included
  if (endTag == "PROJECT")
  {
   id <- as.integer(fields[2])                                             # ID of dataset
   nome <- paste(fields[3:length(fields)],collapse=" ")                    # Project name
   d_project <- rbind(d_project,data.frame(id=id,pname=nome))
  }
  if (endTag == "CRYSTAL")
  {
   id <- as.integer(fields[2])                                             # ID of dataset
   nome <- paste(fields[3:length(fields)],collapse=" ")                    # Crystal name
   d_crystal <- rbind(d_crystal,data.frame(id=id,cname=nome))
  }
  if (endTag == "DATASET")
  {
   id <- as.integer(fields[2])                                             # ID of dataset
   nome <- paste(fields[3:length(fields)],collapse=" ")                    # Dataset name
   d_dataset <- rbind(d_dataset,data.frame(id=id,dname=nome))
  }
  if (endTag == "DCELL")
  {
   id <- as.integer(fields[2])                                             # ID of dataset
   d_dcell <- rbind(d_dcell,data.frame(id=id,a=as.numeric(fields[3]),
                                             b=as.numeric(fields[4]),
                                             c=as.numeric(fields[5]),
                                             alpha=as.numeric(fields[6]),
                                             beta=as.numeric(fields[7]),
                                             gamma=as.numeric(fields[8]))) # Cell parameters (one set for each dataset)
  }
  if (endTag == "DWAVEL")
  {
   id <- as.integer(fields[2])                                             # ID of dataset
   lambda <- as.numeric(fields[3])
   if (is.na(lambda)) lambda <- -1
   d_dwavel <- rbind(d_dwavel,
   #                  data.frame(id=id,lambda=as.numeric(fields[3])))       # Wavelength (one for each dataset)
                     data.frame(id=id,lambda=lambda))       # Wavelength (one for each dataset)
  }
  if (endTag == "COLUMN")
  {
   col_labels <- c(col_labels,fields[2])
   col_types <- c(col_types,fields[3])
   col_min <- c(col_min,fields[4])
   col_max <- c(col_max,fields[5])
   col_id <- c(col_id,fields[6])
   if (messages) print(hdata)
  }
  if (endTag == "BATCH")
  {
   col_batch <- c(col_batch,as.integer(fields[2:length(fields)]))
  }
 }
 l_data$SYMM <- col_symm
 if (length(d_project[,1]) != 1) d_project <- na.omit(d_project)
 rownames(d_project) <- 1:length(d_project[,1])
 l_data$PROJECT <- d_project
 if (length(d_crystal[,1]) != 1) d_crystal <- na.omit(d_crystal)
 rownames(d_crystal) <- 1:length(d_crystal[,1])
 l_data$CRYSTAL <- d_crystal
 if (length(d_dataset[,1]) != 1) d_dataset <- na.omit(d_dataset)
 rownames(d_dataset) <- 1:length(d_dataset[,1])
 l_data$DATASET <- d_dataset
 if (length(d_dcell[,1]) != 1) d_dcell <- na.omit(d_dcell)
 rownames(d_dcell) <- 1:length(d_dcell[,1])
 l_data$DCELL <- d_dcell
 if (length(d_dwavel[,1]) != 1) d_dwavel <- na.omit(d_dwavel)
 for (i in 1:length(d_dwavel[,1]))
 {
  if (!is.na(d_dwavel[i,2])) if (d_dwavel[i,2] == -1) d_dwavel[i,2] <- NA
 }
 rownames(d_dwavel) <- 1:length(d_dwavel[,1])
 l_data$DWAVEL <- d_dwavel
 l_data$COLUMN <- data.frame(labels=I(col_labels),types=I(col_types),min=I(col_min),max=I(col_max),id=I(col_id))  # The I() avoids values being turned into factor levels
 l_data$BATCH <- col_batch

 # History
 col_history <- c()
 bheaderTag=" "
 while (bheaderTag != "MTZBATS" & bheaderTag != "MTZENDOFHEADERS")
 {
  hdata <- readChar(f,80)
  col_history <- c(col_history,hdata)
  bheaderTag <- .tagRead(hdata)[1]
 }
 l_data$HISTORY <- col_history[1:(length(col_history)-1)]

 if (bheaderTag == "MTZENDOFHEADERS")
 {
  hF <- 1
 }
 else
 {
  hF <- 0
 }

 data <- list(l_data,hF)
 return(data)
}

.readMTZHeader <- function(filename,messages=TRUE)
{
 # Reads and output MTZ header (version running independently of readMTZ)

 # Create a connection to binary file
 f <- file(filename, open="rb")
 
 # Reads initial record
 irdata <- .readIR(f)
 errF <- irdata[[1]]
 if (errF != 0) stop("readMTZ: MTZ binary file does not contain initial \"MTZ \" tag. It is either a badly formatted or a wrong MTZ file.")

 # Work out how many reflection records are contained in this MTZ file (= nref X ncols)
 numDataItems <- irdata[[2]] - 20 - 1  

 # Load all reflection records (only use here is for getting to the right point of binary file)
 if (irdata[[3]] == 1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4)
 if (irdata[[3]] == -1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4,endian="big")
 
 # Load header information
 hData <- .readH(f,messages)

 # Close connection
 close(f)

 # Return header stuff as list
 return(hData[[1]])
}

.readBH <- function(f)
{
 # Given a connection f corresponding to the binary MTZ file immediately
 # after "MTZBATS", this function keeps reading the file to extract the
 # batch header content. 
 # It returns a list whose elements include information on each batch header.
 # This is the same as the one documented in the MTZLIB library of CCP4
 # There are still 19 integers which I can't recognize and, therefore, do not output.
 # From the documentation there seem to be included also a character X 3, which I can't
 # find in the binary.

 # Hex version of ascii NULL and "B"
 rawNULL <- as.raw(0) 
 rawB <- as.raw(66)

 # Flag to stop reading batch headers
 bflag <- 0

 # List containing all batch headers
 l_data <- list()

 # Main cycle through all batch headers
 while (bflag != 1)
 {
  # Each batch header is temporarily stored in a list on its own.
  # This is initialized after each cycle
  single_batch <- list()

  batchData <- readChar(f,80)
  batchData <- .tagRead(batchData)              # Initial string starting by "BH" and containing batch number, number of records, number of integers and number of reals
  idx <- as.integer(batchData[2])
  if (batchData[1] == "BH")
  {
   batchData <- readChar(f,80)
   batchData <- .tagRead(batchData)              # Title of batch header
   if (length(batchData) == 1)
   {
    single_batch$TITLE <- ""
   }
   if (length(batchData) > 1)
   {
    parola <- ""
    for (i in 2:length(batchData)) parola <- paste(parola,batchData[i])
    single_batch$TITLE <- parola
   }
  
   # Numbers for batch header, 29 integers and 156 reals
   batchData <- readBin(f,"integer",n=29,size=4)
   single_batch$NWORDS <- batchData[1]
   single_batch$NINTGR <- batchData[2]
   single_batch$NREALS <- batchData[3]
   single_batch$IORTYP <- batchData[4]
   single_batch$LBCELL <- batchData[5:10]
   single_batch$MISFLG <- batchData[11]
   single_batch$JUMPAX <- batchData[12]
   single_batch$NCRYST <- batchData[13]
   single_batch$LCRFLG <- batchData[14]
   single_batch$LDTYPE <- batchData[15]
   single_batch$JSCAX <- batchData[16]
   single_batch$NBSCAL <- batchData[17]
   single_batch$NGONAX <- batchData[18]
   single_batch$LBMFLG <- batchData[19]
   single_batch$NDET <- batchData[20]
   single_batch$LBSETID <- batchData[21]
   single_batch$INTPAD <- batchData[22:29]
   batchData <- readBin(f,"numeric",n=156,size=4)
   single_batch$CELL <- batchData[1:6]
   single_batch$UMAT <- matrix(data=batchData[7:15],nrow=3,ncol=3)
   single_batch$PHIXYZ <- matrix(data=batchData[16:21],nrow=3,ncol=2)
   single_batch$CRYDAT$ETAD <- batchData[22]
   single_batch$CRYDAT$ETADH <- batchData[23]
   single_batch$CRYDAT$ETADV <- batchData[24]
   single_batch$CRYDAT$GENERIC <- matrix(data=batchData[25:33],nrow=3,ncol=3)
   single_batch$DATUM <- batchData[34:36]
   single_batch$PHISTT <- batchData[37]
   single_batch$PHIEND <- batchData[38]
   single_batch$SCANAX <- batchData[39:41]
   single_batch$TIME1 <- batchData[42]
   single_batch$TIME2 <- batchData[43]
   single_batch$BSCALE <- batchData[44]
   single_batch$BBFAC <-  batchData[45]
   single_batch$SDBSCL <- batchData[46]
   single_batch$SDBFAC <- batchData[47]
   single_batch$PHIRANGE <- batchData[48]
   single_batch$BATPAD <- batchData[49:59]
   single_batch$E1 <- batchData[60:62]
   single_batch$E2 <- batchData[63:65]
   single_batch$E3 <- batchData[66:68]
   single_batch$GONPAD <- batchData[69:80]
   single_batch$SOURCE <- batchData[81:83]
   single_batch$S0 <- batchData[84:86]
   single_batch$BEMDAT$ALAMBD <- batchData[87]
   single_batch$BEMDAT$DELAMB <- batchData[88]
   single_batch$BEMDAT$DELCOR <- batchData[89]
   single_batch$BEMDAT$DIVHD  <- batchData[90]
   single_batch$BEMDAT$DIVVD  <- batchData[91]
   single_batch$BEMDAT$rest   <- batchData[92:111]
                                                 # First detector
   single_batch$DX1 <- batchData[112]
   single_batch$THETA1 <- batchData[113] 
   single_batch$DETLM1 <- matrix(data=batchData[114:117],nrow=2,ncol=2)
                                                 # Second detector
   single_batch$DX2 <- batchData[118]
   single_batch$THETA2 <- batchData[119] 
   single_batch$DETLM2 <- matrix(data=batchData[120:123],nrow=2,ncol=2)

   single_batch$DETPAD <- batchData[124:156]

   # Last line
   batchData <- readChar(f,80)
   fields <- .tagRead(batchData)
   if (length(fields) > 1) single_batch$GONLAB <- fields[2:length(fields)]
   if (length(fields) == 1) single_batch$GONLAB <- NULL
  
   # Store this batch header
   l_data[[idx]] <- single_batch
  }
  else
  {
   bflag <- 1
  }
 }
 
 return(l_data)
}

.readMTZ <- function(filename, messages = TRUE){

 ################################################################################
 ## a function for reading an MTZ file                                         ##
 ## N.B. assumes 4 bytes per data item, not sure if this is always true        ##
 ##                                                                            ##
 ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
 ## Started off by David Waterman in June 2009.                                ##
 ## Extended by James Foadi in July 2009                                       ##
 ################################################################################
 
 # Create a connection to binary file
 f <- file(filename, open="rb")
 
 # Reads initial record
 irdata <- .readIR(f)
 errF <- irdata[[1]]
 if (errF != 0) stop("readMTZ: MTZ binary file does not contain initial \"MTZ \" tag. It is either a badly formatted or a wrong MTZ file.")

 # Work out how many reflection records are contained in this MTZ file (= nref X ncols)
 numDataItems <- irdata[[2]] - 20 - 1  

 # Load all reflection records
 if (irdata[[3]] == 1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4)
 if (irdata[[3]] == -1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4,endian="big")
 
 # Load header information
 hData <- .readH(f,messages)
 hF <- hData[[2]]                                    # hF == 1 means no batch headers
 # Turn reflnData into a dataframe where columns are named
 tmpmatrix <- matrix(data=reflnData,nrow=hData[[1]]$NCOL[2],ncol=hData[[1]]$NCOL[1],byrow=TRUE,dimnames=list(NULL,hData[[1]]$COLUMN$labels))
 tmpdataframe <- as.data.frame(tmpmatrix)
 reflnData <- tmpdataframe
 rm(tmpmatrix,tmpdataframe)

 # If no batch headers are contained, program ends here
 if (hF == 1)
 {
  close(f)                  # Close connection
  if (hData[[1]]$NCOL[3] > 0) warning("MTZ appears to be a multi-record file, but it does not include batch headers.")

  # Output data are packed in a list
  data <- list(reflections=reflnData,header=hData[[1]],batch_header=NULL)

  # Create "mtz" class
  #class(data) <- "mtz"

  return(data)
 }

 # Batch headers loop
 bhData <- .readBH(f)

 # Close connection
 close(f)

 # Output data are packed in a list
 data <- list(reflections=reflnData,header=hData[[1]],batch_header=bhData)

 # Create "mtz" class
 #class(data) <- "mtz"

 return(data)
}
# Functions to read MTZ ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ###   END


# Functions to write MTZ ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### START
.writeMTZ <- function(data,filename,title="Untitled")
{

 ################################################################################
 ## a function for writing an MTZ file                                         ##
 ## N.B. assumes 4 bytes per data item, not sure if this is always true        ##
 ##                                                                            ##
 ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
 ################################################################################
 
 # Create a connection to binary file
 f <- file(filename, open="wb")

 # Hex version of ascii "MTZ "
 rawNULL <- as.raw(0)
 rawBlank <- charToRaw(" ")
 rawMTZ <- c(as.raw(77),as.raw(84),as.raw(90),rawBlank)

 # Write MTZ file signature
 writeBin(rawMTZ,f)

 # Determine data length
 dl <- length(data[[1]][,1])
 cl <- length(data[[1]])
 data_length <- dl*cl

 # Location of header
 #headerLoc <- 80+data_length*4+1
 headerLoc <- data_length+20+1
 storage.mode(headerLoc) <- "integer"  # Turn headerLoc into integer for correct binary storage
 writeBin(headerLoc,f,size=4)

 # Machine stamp (to be improved)
 machineStamp <- "DA"
 writeChar(machineStamp,f,nchars=2)

 # Fill with NULLs for 69 bytes (2 to complete machine stamp, plus 68 to get to a total of 80. We write one less as this is automatically added when connection is closed)
 for (i in 1:69) writeBin(rawNULL,f)

 # Data (data frame -> matrix -> transpose(matrix) -> vector -> write to binary file)
 #numDataItems <- (headerLoc-80-1)/4
 numDataItems <- data_length
 linedata <- as.vector(t(as.matrix(data[[1]])))
 writeBin(linedata,f,size=4)
 rm(linedata)

 # Build header 80-characters lines and write to binary file
 fullhdata <- ""
 hdata <- "VERS MTZ:V1.1                                                                   "
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # MTZ stamp
 hdata <- paste("TITLE",title)
 nblanks <- 80-nchar(hdata)
 for (i in 1:nblanks) hdata <- paste(hdata," ",sep="")
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # TITLE
 hdata <- sprintf("NCOL%9d%13d%9d                                             ",data[[2]]$NCOL[1],data[[2]]$NCOL[2],data[[2]]$NCOL[3])
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # NCOL
 hdata <- sprintf("CELL %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f               ",data[[2]]$CELL[1],data[[2]]$CELL[2],data[[2]]$CELL[3],
                                                                           data[[2]]$CELL[4],data[[2]]$CELL[5],data[[2]]$CELL[6])
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # CELL
 hdata <- sprintf("SORT %4d%4d%4d%4d%4d                                                       ",data[[2]]$SORT[1],data[[2]]$SORT[2],data[[2]]$SORT[3],
                                                                                                data[[2]]$SORT[4],data[[2]]$SORT[5])
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # SORT
 hdata <- sprintf("SYMINF %3d%3d %1s%6d",data[[2]]$SYMINF[[1]],data[[2]]$SYMINF[[2]],data[[2]]$SYMINF[[3]],data[[2]]$SYMINF[[4]])
 nblanks <- 23-nchar(data[[2]]$SYMINF[[5]])
 stmp <- ""
 for (i in 1:nblanks) stmp <- paste(stmp," ",sep="")
 stmp <- paste(stmp,data[[2]]$SYMINF[[5]],sep="")
 hdata <- paste(hdata,stmp,sprintf("%6s                              ",data[[2]]$SYMINF[[6]]),sep="")
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # SYMINF
 for (i in 1:length(data[[2]]$SYMM))
 {
  hdata <- data[[2]]$SYMM[i]
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # SYMM
 }
 hdata <- sprintf("RESO %8.6f   %8.6f                                                        ",data[[2]]$RESO[1],data[[2]]$RESO[2])
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # RESO
 hdata <- "VALM NAN                                                                        "
 fullhdata <- paste(fullhdata,hdata,sep="")                                                  # VALM
 for (i in 1:length(data[[2]]$COLUMN[,1]))
 {
  linea <- c()
  linea <- c(linea,data[[2]]$COLUMN$labels[i],data[[2]]$COLUMN$types[i],data[[2]]$COLUMN$min[i],data[[2]]$COLUMN$max[i],data[[2]]$COLUMN$id[i])
  hdata <- sprintf("COLUMN %-31s%1s%18s%18s  %3s",linea[1],linea[2],linea[3],linea[4],linea[5])
  fullhdata <- paste(fullhdata,hdata,sep="")                                                 # COLUMN
 }
 hdata <- sprintf("NDIF %8d                                                                   ",data[[2]]$NDIF)
 fullhdata <- paste(fullhdata,hdata,sep="")                                                    # NDIF
 for (i in 1:data[[2]]$NDIF)
 {
  hdata <- sprintf("PROJECT%8d %-64s",data[[2]]$PROJECT$id[i],data[[2]]$PROJECT$pname[i])
  fullhdata <- paste(fullhdata,hdata,sep="")
  hdata <- sprintf("CRYSTAL%8d %-64s",data[[2]]$CRYSTAL$id[i],data[[2]]$CRYSTAL$cname[i])
  fullhdata <- paste(fullhdata,hdata,sep="")
  hdata <- sprintf("DATASET%8d %-64s",data[[2]]$DATASET$id[i],data[[2]]$DATASET$dname[i])
  fullhdata <- paste(fullhdata,hdata,sep="")
  hdata <- sprintf("DCELL%10d %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f    ",data[[2]]$DCELL$id[i],data[[2]]$DCELL$a[i],
                                                                         data[[2]]$DCELL$b[i],data[[2]]$DCELL$c[i],
                                                                         data[[2]]$DCELL$alpha[i],data[[2]]$DCELL$beta[i],
                                                                         data[[2]]$DCELL$gamma[i])
  fullhdata <- paste(fullhdata,hdata,sep="")
  hdata <- sprintf("DWAVEL%9d %10.5f                                                      ",data[[2]]$DWAVEL$id[i],data[[2]]$DWAVEL$lambda[i])
  fullhdata <- paste(fullhdata,hdata,sep="")
 }
 plength <- length(data[[2]]$BATCH)
 dnd <- 0
 while (plength > 0)
 {
  dnd <- dnd+1
  rem <- plength%%12
  istart <- (dnd-1)*12+1
  iend <- istart+11
  hdata <- sprintf("BATCH ")
  if (plength >= 12)
  {
   for (i in istart:iend) hdata <- paste(hdata,sprintf("%6d",data[[2]]$BATCH[i]),sep="")
   hdata <- paste(hdata,sprintf("  "),sep="")
   fullhdata <- paste(fullhdata,hdata,sep="")
  }
  if (plength < 12)
  {
   for (i in istart:length(data[[2]]$BATCH)) hdata <- paste(hdata,sprintf("%6d",data[[2]]$BATCH[i]),sep="")
   iend <- 12-rem
   for (i in 1:iend) hdata <- paste(hdata,sprintf("      "),sep="")
   hdata <- paste(hdata,sprintf("  "),sep="")
   fullhdata <- paste(fullhdata,hdata,sep="")
  } 
  plength <- plength-12
 }
 hdata <- sprintf("END                                                                             ")
 fullhdata <- paste(fullhdata,hdata,sep="")
 #hdata <- sprintf("MTZHIST%4d                                                                     ",(length(data[[2]]$HISTORY)-1))
 #fullhdata <- paste(fullhdata,hdata,sep="")
 for (i in 1:length(data[[2]]$HISTORY))
 {
  hdata <- data[[2]]$HISTORY[i]
  fullhdata <- paste(fullhdata,hdata,sep="")                                                # HISTORY
 }

 # End of MTZ or beginning of batch headers
 if (data[[2]]$NCOL[3] == 0)                                                                 # MTZENDOFHEADERS
 {
  hdata <- sprintf("MTZENDOFHEADERS                                                                 ")
  fullhdata <- paste(fullhdata,hdata,sep="")
  writeChar(fullhdata,f,eos=NULL)

  # Close connection
  close(f)
 }

 if (data[[2]]$NCOL[3] != 0)
 {
  # Carry on with batch header if it is an unmerged file
  hdata <- sprintf("MTZBATS                                                                         ")
  fullhdata <- paste(fullhdata,hdata,sep="")
  writeChar(fullhdata,f,eos=NULL)

  # All BATCH HEADER stuff is written from here onward 
  nbatches <- length(data[[2]]$BATCH)
  for (num in data[[2]]$BATCH[1]:data[[2]]$BATCH[nbatches])
  {
   hdata <- sprintf("BH %8d%8d%8d%8d                                             ",num,data[[3]][[num]]$NWORDS,data[[3]][[num]]$NINTGR,data[[3]][[num]]$NREALS)
   writeChar(hdata,f,eos=NULL)
   hdata <- paste("TITLE",data[[3]][[num]]$TITLE)                                                             # TITLE
   nblanks <- 80-nchar(hdata)
   for (i in 1:nblanks) hdata <- paste(hdata," ",sep="")
   writeChar(hdata,f,eos=NULL)
   ## INTEGERS
   writeBin(data[[3]][[num]]$NWORDS,f,size=4)                                                                 # NWORDS
   writeBin(data[[3]][[num]]$NINTGR,f,size=4)                                                                 # NINTGR
   writeBin(data[[3]][[num]]$NREALS,f,size=4)                                                                 # NREALS
   writeBin(data[[3]][[num]]$IORTYP,f,size=4)                                                                 # IORTYP
   writeBin(data[[3]][[num]]$LBCELL,f,size=4)                                                                 # LBCELL (6)
   writeBin(data[[3]][[num]]$MISFLG,f,size=4)                                                                 # MISFLG
   writeBin(data[[3]][[num]]$JUMPAX,f,size=4)                                                                 # JUMPAX
   writeBin(data[[3]][[num]]$NCRYST,f,size=4)                                                                 # NCRYST
   writeBin(data[[3]][[num]]$LCRFLG,f,size=4)                                                                 # LCRFLG
   writeBin(data[[3]][[num]]$LDTYPE,f,size=4)                                                                 # LDTYPE
   writeBin(data[[3]][[num]]$JSCAX,f,size=4)                                                                  # JSCAX
   writeBin(data[[3]][[num]]$NBSCAL,f,size=4)                                                                 # NBSCAL
   writeBin(data[[3]][[num]]$NGONAX,f,size=4)                                                                 # NGONAX
   writeBin(data[[3]][[num]]$LBMFLG,f,size=4)                                                                 # LBMFLG
   writeBin(data[[3]][[num]]$NDET,f,size=4)                                                                   # NDET
   writeBin(data[[3]][[num]]$LBSETID,f,size=4)                                                                # LBSETID
   writeBin(data[[3]][[num]]$INTPAD,f,size=4)                                                                 # INTPAD (8)
   ## REALS
   writeBin(data[[3]][[num]]$CELL,f,size=4)                                                                   # CELL (6)
   writeBin(as.vector(data[[3]][[num]]$UMAT),f,size=4)                                                        # UMAT (9)
   writeBin(as.vector(data[[3]][[num]]$PHIXYZ),f,size=4)                                                      # PHIXYZ (6)
   writeBin(data[[3]][[num]]$CRYDAT$ETAD,f,size=4)                                                            # CRYDAT$ETAD
   writeBin(data[[3]][[num]]$CRYDAT$ETADH,f,size=4)                                                           # CRYDAT$ETADH
   writeBin(data[[3]][[num]]$CRYDAT$ETADV,f,size=4)                                                           # CRYDAT$ETADH
   writeBin(as.vector(data[[3]][[num]]$CRYDAT$GENERIC),f,size=4)                                              # CRYDAT$GENERIC (9)
   writeBin(data[[3]][[num]]$DATUM,f,size=4)                                                                  # DATUM (3)
   writeBin(data[[3]][[num]]$PHISTT,f,size=4)                                                                 # PHISTT
   writeBin(data[[3]][[num]]$PHIEND,f,size=4)                                                                 # PHIEND
   writeBin(data[[3]][[num]]$SCANAX,f,size=4)                                                                 # SCANAX (3)
   writeBin(data[[3]][[num]]$TIME1,f,size=4)                                                                  # TIME1
   writeBin(data[[3]][[num]]$TIME2,f,size=4)                                                                  # TIME2
   writeBin(data[[3]][[num]]$BSCALE,f,size=4)                                                                 # BSCALE
   writeBin(data[[3]][[num]]$BBFAC,f,size=4)                                                                  # BBFAC
   writeBin(data[[3]][[num]]$SDBSCL,f,size=4)                                                                 # SDBSCL
   writeBin(data[[3]][[num]]$SDBFAC,f,size=4)                                                                 # SDBFAC
   writeBin(data[[3]][[num]]$PHIRANGE,f,size=4)                                                               # PHIRANGE
   writeBin(data[[3]][[num]]$BATPAD,f,size=4)                                                                 # BATPAD (11)
   writeBin(data[[3]][[num]]$E1,f,size=4)                                                                     # E1 (3)
   writeBin(data[[3]][[num]]$E2,f,size=4)                                                                     # E2 (3)
   writeBin(data[[3]][[num]]$E3,f,size=4)                                                                     # E3 (3)
   writeBin(data[[3]][[num]]$GONPAD,f,size=4)                                                                 # GONPAD (12)
   writeBin(data[[3]][[num]]$SOURCE,f,size=4)                                                                 # SOURCE (3)
   writeBin(data[[3]][[num]]$S0,f,size=4)                                                                     # S0 (3)
   writeBin(data[[3]][[num]]$BEMDAT$ALAMBD,f,size=4)                                                          # BEMDAT$ALAMBD
   writeBin(data[[3]][[num]]$BEMDAT$DELAMB,f,size=4)                                                          # BEMDAT$DELAMB
   writeBin(data[[3]][[num]]$BEMDAT$DELCOR,f,size=4)                                                          # BEMDAT$DELCOR
   writeBin(data[[3]][[num]]$BEMDAT$DIVHD,f,size=4)                                                           # BEMDAT$DIVHD
   writeBin(data[[3]][[num]]$BEMDAT$DIVVD,f,size=4)                                                           # BEMDAT$DIVVD
   writeBin(data[[3]][[num]]$BEMDAT$rest,f,size=4)                                                            # rest of BEMDAT (20)
   writeBin(data[[3]][[num]]$DX1,f,size=4)                                                                    # DX1
   writeBin(data[[3]][[num]]$THETA1,f,size=4)                                                                 # THETA1
   writeBin(as.vector(data[[3]][[num]]$DETLM1),f,size=4)                                                      # DETLM1 (4)
   writeBin(data[[3]][[num]]$DX2,f,size=4)                                                                    # DX2
   writeBin(data[[3]][[num]]$THETA2,f,size=4)                                                                 # THETA2
   writeBin(as.vector(data[[3]][[num]]$DETLM2),f,size=4)                                                      # DETLM2 (4)
   writeBin(data[[3]][[num]]$DETPAD,f,size=4)                                                                 # DETPAD (33)
 
   # Last line "BHCH"
   if (is.null(data[[3]][[num]]$GONLAB)) batchData <- "BHCH                                                                            "
   if (!is.null(data[[3]][[num]]$GONLAB))
   {
    lengthGONLAB <- length(data[[3]][[num]]$GONLAB)
    batchData <- "BHCH "
    for (i in 1:lengthGONLAB) batchData <- paste(batchData,sprintf("%-9s%-9s%-9s",
                                                 data[[3]][[num]]$GONLAB[1],data[[3]][[num]]$GONLAB[2],data[[3]][[num]]$GONLAB[3]),sep="")
    nblanks <- 80-nchar(hdata)
    for (i in 1:nblanks) batchData <- paste(batchData," ",sep="")
   }
   writeChar(batchData,f,eos=NULL)
  }
  # End of MTZ file
  batchData <- "MTZENDOFHEADERS                                                                 "
  writeChar(batchData,f,eos=NULL)

  # Close connection
  close(f)
 }
}
# Functions to write MTZ ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ### ### MTZ ###   END

#
.squared_resolution_coeffs <- function(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
{
 # Given cell parameters return coefficients for the expression of squared resolution:
 #
 #  s^2 = a*h^2+b*k^2+c*l^2+2*d*h*k+2*e*h*l+2*f*k*l

 ctmp <- unname(.dcl(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma))
 den <- 1-ctmp[4]^2-ctmp[5]^2-ctmp[6]^2+2*ctmp[4]*ctmp[5]*ctmp[6]
 a <- ctmp[1]^2/(cell_a^2*den)
 b <- ctmp[2]^2/(cell_b^2*den)
 c <- ctmp[3]^2/(cell_c^2*den)
 d <- (ctmp[4]*ctmp[5]-ctmp[6])/(cell_a*cell_b*den)
 e <- (ctmp[4]*ctmp[6]-ctmp[5])/(cell_a*cell_c*den)
 f <- (ctmp[5]*ctmp[6]-ctmp[4])/(cell_b*cell_c*den)

 return(c(a,b,c,d,e,f))
}

#
.generate_reflections_up_to_resolution <- function(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,reso)
{
 # Generate all reflections (Miller indices) with resolutions lower or equal
 # to reso.

 # Inverse, squared resolution
 ss <- (1/reso)^2

 # Coefficients of squared resolution function
 ctmp <- .squared_resolution_coeffs(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
 a <- ctmp[1]
 b <- ctmp[2]
 c <- ctmp[3]
 d <- ctmp[4]
 e <- ctmp[5]
 f <- ctmp[6]

 # Find max h, k, l to generate all (h,k,l) within minimal box

 # Coefficient determinant (common to both h, and k and l)
 D <- a*(b*c-f^2)-d*(c*d-e*f)+e*(d*f-b*e)

 # Max h
 A <- (b*c-f^2)/D 
 B <- (e*f-c*d)/D
 C <- (d*f-b*e)/D
 t <- sqrt(ss/(a*A^2+b*B^2+c*C^2+2*d*A*B+2*e*A*C+2*f*B*C))
 max_h <- ceiling(A*t)

 # Max k
 A <- (e*f-c*d)/D 
 B <- (a*c-e^2)/D
 C <- (d*e-a*f)/D
 t <- sqrt(ss/(a*A^2+b*B^2+c*C^2+2*d*A*B+2*e*A*C+2*f*B*C))
 max_k <- ceiling(B*t)

 # Max l
 A <- (d*f-b*e)/D 
 B <- (d*e-a*f)/D
 C <- (a*b-d^2)/D
 t <- sqrt(ss/(a*A^2+b*B^2+c*C^2+2*d*A*B+2*e*A*C+2*f*B*C))
 max_l <- ceiling(C*t)

 # Dataframe with all h, k, l
 hkl <- expand.grid(h=-max_h:max_h,k=-max_k:max_k,l=-max_l:max_l)

 # Strictly reflections with resolution lower than reso
 new_hkl <- hkl[1/.d_hkl(hkl$h,hkl$k,hkl$l,cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma) <= 1/reso,]
 
 # Add inverse of resolution (s) as last column of data frame
 s <- 1/.d_hkl(new_hkl$h,new_hkl$k,new_hkl$l,cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
 new_hkl <- cbind(new_hkl,data.frame(s=s))

 return(new_hkl)
}

#
.apply_hklasu <- function(sym_number,hkl,setting=1)
{
 # This auxiliary function starts from all reflections, up to a given resolution, and delete all
 # reflectios outside asymmetric reciprocal unit. It also deletes systematic absences, if present.
 # hkl is a data frame or matrix whose first 3 columns are h, k, l. It has to have at least 4 columns.

 # There are only 9 asymmetric units for Miller indices
 if (sym_number == 1 | sym_number ==2)
 {
  tmp <- hkl[hkl$l >= 0,]
  tmp[tmp$h == 0 & tmp$l == 0 & tmp$k < 0,4] <- NA
  tmp <- na.omit(tmp)
 }
 if (sym_number >= 3 & sym_number <= 15)
 {
  tmp <- hkl[hkl$k >= 0 & hkl$l >=0,]
  tmp[tmp$l == 0 & tmp$h < 0,4] <- NA
  tmp <- na.omit(tmp)
 }
 if (sym_number >= 16 & sym_number <= 74)
 {
  tmp <- hkl[hkl$h >= 0 & hkl$k >= 0 & hkl$l >= 0,]
 }
 if ((sym_number >= 75 & sym_number <= 88) | (sym_number >= 168 & sym_number <= 176))
 {
  tmp <- hkl[hkl$l >= 0 & hkl$h >= 0 & hkl$k >= 0,]
  tmp[tmp$h == 0 & tmp$k > 0,4] <- NA
  tmp <- na.omit(tmp)
 }
 if ((sym_number >= 89 & sym_number <= 142) | (sym_number >= 177 & sym_number <= 194))
 {
  tmp <- hkl[hkl$h >= 0 & hkl$k >= hkl$h & hkl$l >= 0,]
  tmp[tmp$h == 0 & tmp$k > 0,4] <- NA
  tmp <- na.omit(tmp)
 }
 if (sym_number >= 143 & sym_number <= 148)
 {
  tmp <- hkl[hkl$h >= 0 & hkl$k >= 0,]
  tmp[tmp$k == 0 & tmp$l < 0,4] <- NA
  tmp[tmp$h == 0 & tmp$l <= 0,4] <- NA
  tmp <- na.omit(tmp)
 }
 if (sym_number >= 149 & sym_number <= 167)
 {
  tmp <- hkl[hkl$h >= 0 & hkl$k >= hkl$h,]
  tmp[tmp$h == 0 & tmp$l < 0,4] <- NA
  tmp <- na.omit(tmp)
 }
 if (sym_number >= 195 & sym_number <= 206)
 {
  tmp <- hkl[hkl$h >= 0 & hkl$k >= hkl$h & hkl$l >= hkl$h,]
  tmp[tmp$h < tmp$k,4] <- NA
  tmp <- na.omit(tmp)
 }
 if (sym_number >= 207 & sym_number <= 230)
 {
  tmp <- hkl[hkl$h >= 0 & hkl$k >= hkl$h & hkl$l >= hkl$k,4] <- NA
  tmp <- na.omit(tmp)
 }

 # Now delete systematic absences, if present in arrays
 tmp <- .find_systematic_absences(tmp,sym_number,setting)

 return(tmp)
}

#
.find_systematic_absences <- function(hkl,sym_number,setting)
{
 # Given a data frame hkl with first three columns being H, K, L (Miller indices)
 # this function finds out if the space group of number sym_number includes any
 # systematic absences and, if any, deletes them from hkl, returning the resulting
 # data frame.

 # Make a m X 3 matrix of original Miller indices
 hkl2 <- as.matrix(hkl[,1:3])

 # Delete incorrect Miller indices if they are systematic absences
 idx <- sysabs(hkl2,sym_number,setting)
 hkl <- hkl[idx,]

 return(hkl)
}

#
.check_systematic_absences <- function(hkl,sym_number,setting=1)
{
 # Given a data frame hkl with first three columns being H, K, L (Miller indices)
 # this function finds out if the space group of number sym_number includes any
 # systematic absences. It next check whether hkl obeys any found systematic absence.
 # It returns TRUE if no systematic absences are included by the space group, or if
 # hkl does not contain systematic absences; it returns FALSE if hkl includes systematic
 # absences which shuld not be there.
 # Careful!!! This function does not check if setting is correct, as this should be done 
 # somewhere else.

 # Final answer
 fans <- TRUE

 # Make a m X 3 matrix of original Miller indices
 hkl2 <- as.matrix(hkl[,1:3])

 # Delete incorrect Miller indices if they are systematic absences
 idx <- sysabs(hkl2,sym_number,setting)
 if (length(idx) != length(hkl2[,1])) fans <- FALSE

 return(fans)
}
