###################################################################
###################################################################
### This module is part of the cRy package for crystallography. ###
### Authors: J. Foadi & D. G. Waterman                          ###
### MPL - Imperial College London / Diamond Light Source Ltd    ###
### CCP4 - Research Complex at Harwell                          ###
###                                                             ###
### Auxiliary functions used throughout the code. They are the  ### 
### first to be collated in the package.                        ###
###################################################################
###################################################################



#################################### MOSTLY RELATED TO SYMMETRY (OLD symmetry_dot.R) #############################################

## Translate space group number of symbol into other formats
#
.translate_SG <- function(value,SG_in="number",SG_out="xHM",setting=1)
{
 # Input: space group in one of known formats; output: space group in another
 # known format. Possible formats are:
 # 1) Space group number
 # 2) Hall symbol (e.g. ' P 2yb (z,x,y)')
 # 3) Extended Hermann-Maguin symbol (e.g. 'P 1 1 21')
 # If more than one setting is implied in an unbiguous way in the input value,
 # then the first setting will be selected by default for the output value, unless
 # argument "setting" is set to another value
 # Input: SG_in is a string representing the space group format. Possible values are:
 # 1) "number"
 # 2) "ccp4" 
 # 3) "Hall" 
 # 4) "xHM" 
 # 5) "old" 
 # Input: SG_out is a string representing the space group format needed as output. Possible
 #        values are the same as those for SG_in
 # Input: setting is a number like 1,2,...If for a same symbol there are more settings, use this. 
 
 # A few checks
 if (SG_in != "number" & SG_in != "ccp4" & SG_in != "Hall" & SG_in != "xHM" & SG_in != "old") stop("Wrong SG_in string. Valid strings are: number ccp4 Hall xHM old")
 if (SG_out != "number" & SG_out != "ccp4" & SG_out != "Hall" & SG_out != "xHM" & SG_out != "old") stop("Wrong SG_out string. Valid strings are: number ccp4 Hall xHM old")

 # If number is not in the range 1:230 stop
 if (SG_in == "number")
 {
  idx <- which(value == 1:230)
  if (length(idx) == 0)
  {
   msg <- "There is not a space group associated with the input number."
   return(list(msg=msg,ans=FALSE))
  }
 }

 # If input is "number", turn into a string
 if (SG_in == "number") value <- paste("number ",value)
 if (SG_in == "ccp4") value <- paste("symbol ccp4",value)

 # Complete string for non-number cases
 if (SG_in == "Hall") value <- paste("symbol Hall '",value,"'",sep="")
 if (SG_in == "xHM") value <- paste("symbol xHM  '",value,"'",sep="")
 if (SG_in == "old") value <- paste("symbol old  '",value,"'",sep="")

 # Content of syminfo.lib is automatically loaded with cRy package
 bsg <- grep(value,syminfo,fixed=TRUE)
 

 # Select correct one in the "number" case
 if (SG_in == "number")
 {
  bsg2 <- c()
  for (i in 1:length(bsg))
  {
   numero1 <- strsplit(syminfo[bsg[i]],"  ")[[1]][2]
   numero2 <- strsplit(value,"  ")[[1]][2]
   if (numero1 == numero2) bsg2 <- c(bsg2,bsg[i])
  }
  bsg <- bsg2
  rm(bsg2)
 }
 #print(paste("There are ",length(bsg)," settings for this space group"))

 # Select case based on setting
 if (length(bsg) < setting)
 {
  msg <- "Something wrong in your input:"
  msg <- paste(msg,"   1) the symbol or number input for this space group does not exist",sep="\n")
  msg <- paste(msg,"   2) if your inpur was a number, perhaps for this space group there are not that many settings",sep="\n")
  return(list(msg=msg,ans=FALSE))
 }
 bsg <- bsg[setting]

 bini <- grep("begin_spacegroup",syminfo)
 prima <- bini[length(bini[bini < bsg])]
 dopo  <- bini[bini > bsg][1]-1   # Add 1 for those cases when prima = dopo
 if (is.na(dopo)) dopo <- length(syminfo)
 if (SG_out == "number")
 {
  key <- "number"
 }
 if (SG_out != "number")
 {
  key <- paste("symbol",SG_out)
 }
 tmp <- syminfo[prima:dopo][grep(key,syminfo[prima:dopo])]
 if (key == "number") translated_value <- strsplit(tmp,"  ")[[1]][2] 
 if (key != "number")
 {
  tmp2 <- strsplit(tmp," ")[[1]]
  if (tmp2[2] == "ccp4") translated_value <- tmp2[3]
  if (tmp2[2] == "Hall") translated_value <- strsplit(tmp,"'")[[1]][2]
  if (tmp2[2] == "xHM") translated_value <- strsplit(tmp,"'")[[1]][2]
  if (tmp2[2] == "old") translated_value <- strsplit(tmp,"'")[[1]][2]
  rm(tmp2)
 }
 rm(tmp)

 # If output requires number, turn character into numeric
 if (SG_out == "number" | SG_out == "ccp4") translated_value <- as.integer(translated_value)

 return(list(msg=translated_value,ans=TRUE))
}

## Given the extended Hermann-Maguin symbol it returns list of symmetry operators in string form
#
.syminfo_to_op_xyz_list <- function(SG)
{
 # Input is spacegroup symbol in xHM format.
 # Output is a list of 2 vector of strings,
 # The first string describes point group and translation; the
 # second string describes cell centering 

 # Extract full symmetry information 
 data <- .extract_symmetry_info(SG)

 # Extract "symop" bit
 tmp2 <- strsplit(data$SYMOP," ")
 symop_xyz <- c()
 for (i in 1:length(tmp2))
 {
  symop_xyz <- c(symop_xyz,tmp2[[i]][2])
 }
 
 # Extract "cenop" bit
 tmp2 <- strsplit(data$CENOP," ")
 cenop_xyz  <- c()
 for (i in 1:length(tmp2))
 {
  cenop_xyz <- c(cenop_xyz,tmp2[[i]][2])
 }
 op_xyz_list <- list(symop_xyz,cenop_xyz)

 return(op_xyz_list)
}

## Given the list of symmetry operators in string form, it returns it in matrix form
#
.op_xyz_list_to_matrix_list <- function(op_xyz_list)
{
 # Input: "syminfo_to_op_xyz_list" output, i.e. a list of 2 character vectors.
 # The first one contains symmetry operators in x,y,z format; the second one centering operators
 # in x,y,z format.
 # Returns a list consisting of 3 lists. The first list contains 3X3 point group matrices; 
 # the second list contains the same number of 3X1 translation vectors. First matrix is always the
 # identity matrix, the first translation vector is always the null vector. the third
 # list consists of centering vectors; the first centering vector is always the null
 # vector. To summarize, the output looks like the following:
 # [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
 # where:
 # I            = identity    3X3 matrix
 # 0            = null        3X1 vector
 # M2,M3,...,Mn = point group 3X3 matrices
 # V2,V3,...,Cn = translation 3X1 vectors
 # C2,C3,...,Cm = centering   3X1 vectors

 # Create empty lists
 matrix_list <- list() 
 vector_list <- list()
 centering_list <- list()

 # Point group matrices and translation vectors
 for (i in 1:length(op_xyz_list[[1]]))
 {
  data <- .op_xyz_to_matrix(op_xyz_list[[1]][i])
  matrix_list[i] <- list(data[[1]])
  vector_list[i] <- list(data[[2]])
 }

 # Centering vectors
 for (i in 1:length(op_xyz_list[[2]]))
 {
  data <- .op_xyz_to_matrix(op_xyz_list[[2]][i])
  centering_list[i] <- list(data[[2]])
 }

 overall_list <- list(PG=matrix_list,T=vector_list,C=centering_list)

 return(overall_list)
}

## Given a symmetry operator in string form, it returns it in matrix form
#
.op_xyz_to_matrix <- function(op_xyz)
{
 # Reads in a symmetry or centering operator in character form and output it in
 # matrix or vector form. If input is a symmetry operator, output is a list of
 # a 3X3 matrix and a 3X1 vector; if input is a centering operator, output is
 # still a matrix and a vector, but the matrix is always the identity matrix.

 ltmp <- strsplit(op_xyz,",")
 if (substr(ltmp[[1]][1],1,1) != "-" & substr(ltmp[[1]][1],1,1) != "+")
 {
  stmp <- paste("+",ltmp[[1]][1],sep="")
  ltmp[[1]][1] <- stmp
 }
 if (substr(ltmp[[1]][2],1,1) != "-" & substr(ltmp[[1]][2],1,1) != "+")
 {
  stmp <- paste("+",ltmp[[1]][2],sep="")
  ltmp[[1]][2] <- stmp
 }
 if (substr(ltmp[[1]][3],1,1) != "-" & substr(ltmp[[1]][3],1,1) != "+")
 {
  stmp <- paste("+",ltmp[[1]][3],sep="")
  ltmp[[1]][3] <- stmp
 }

 # Get 3X3 point group matrix
 m <- matrix(rep(0,times=9),nrow=3,ncol=3)
 for (j in 1:length(ltmp[[1]]))
 {
  if (grepl("\\+x",ltmp[[1]][j])) m[j,1] <-  1
  if (grepl("\\-x",ltmp[[1]][j])) m[j,1] <- -1
  if (grepl("\\+y",ltmp[[1]][j])) m[j,2] <-  1
  if (grepl("\\-y",ltmp[[1]][j])) m[j,2] <- -1
  if (grepl("\\+z",ltmp[[1]][j])) m[j,3] <-  1
  if (grepl("\\-z",ltmp[[1]][j])) m[j,3] <- -1
 }

 v <- rep(0,times=3)
 for (j in 1:length(ltmp[[1]]))
 {
  bb <- strsplit(ltmp[[1]][j],"/")
  if (length(bb[[1]]) > 1)
  {
   v[j] <- as.integer(substr(bb[[1]][1],nchar(bb[[1]][1])-1,nchar(bb[[1]][1])))/
           as.integer(bb[[1]][2])
  }
  if (length(bb[[1]]) <= 1)
  {
   v[j] <- 0
  }
 }

 # Final list
 final_list <- list(m,v) 

 return(final_list)
}

## Given space group symbol in Hermann-Maguin form, it returns the list of symmetry operations in matrix form
#
.syminfo_to_matrix_list <- function(SG)
{
 # Input: space group name in xHM format.
 # Returns a list consisting of 3 lists. The first list contains 3X3 point group matrices; 
 # the second list contains the same number of 3X1 translation vectors. First matrix is always the
 # identity matrix, the first translation vector is always the null vector. the third
 # list consists of centering vectors; the first centering vector is always the null
 # vector. To summarize, the output looks like the following:
 # [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
 # where:
 # I            = identity    3X3 matrix
 # 0            = null        3X1 vector
 # M2,M3,...,Mn = point group 3X3 matrices
 # V2,V3,...,Cn = translation 3X1 vectors
 # C2,C3,...,Cm = centering   3X1 vectors
 #
 # This function is simply a wrapper for .op_xyz_list_to_matrix_list.

 # 3 functions called in one line (cool, isn't it?)
 return(.op_xyz_list_to_matrix_list(.syminfo_to_op_xyz_list(SG)))
}

## Extracts symmetry information corresponding to space group SG (Hermann-Maguin form), and returns it in a list
#
.extract_symmetry_info <- function(SG)
{
 # Input is spacegroup symbol in xHM format.
 # Output is a list with all symmetry information
 # for the specific SG group, as contained in syminfo.lib.

 # Content of syminfo.lib is automatically loaded with cRy package
 bsg <- grep(SG,syminfo)
 bini <- grep("begin_spacegroup",syminfo)
 bend <- grep("end_spacegroup",syminfo)
 prima <- bini[length(bini[bini < bsg[1]])]
 dopo  <- bend[bend > bsg[1]][1]

 # Extract info in string format
 infostring <- list()
 infostring$NUMBER       <- syminfo[prima+ 1] 
 infostring$BASISOP      <- syminfo[prima+ 2] 
 infostring$CCP4         <- syminfo[prima+ 3] 
 infostring$HALL         <- syminfo[prima+ 4] 
 infostring$XHM          <- syminfo[prima+ 5] 
 infostring$OLD          <- syminfo[prima+ 6] 
 infostring$LAUE         <- syminfo[prima+ 7] 
 infostring$PATT         <- syminfo[prima+ 8] 
 infostring$PGRP         <- syminfo[prima+ 9] 
 infostring$HKLASU       <- syminfo[prima+10] 
 infostring$MAPASU_CCP4  <- syminfo[prima+11] 
 infostring$MAPASU_ZERO  <- syminfo[prima+12] 
 infostring$MAPASU_NONZ  <- syminfo[prima+13] 
 infostring$CHESHIRE     <- syminfo[prima+14] 
 infostring$SYMOP        <- syminfo[prima:dopo][grep("symop",syminfo[prima:dopo])]
 infostring$CENOP        <- syminfo[prima:dopo][grep("cenop",syminfo[prima:dopo])]
 
 return(infostring)
}

## Returns all symmetry operators in string format (the kind SYMM)
# 
.full_symm_strings <- function(SG)
{
 # Extract all symmetry information
 sinfo <- .extract_symmetry_info(SG)

 # 1. Point symmetry
 pgstring <- sinfo$SYMOP
 npgops <- length(pgstring)

 # 2. Centring
 cstring <- sinfo$CENOP
 ncops <- length(cstring)

 # Produce npgops X ncops vectors
 final_symm <- c()

 # Double loop
 for (ssym in pgstring)
 {
  # 3 character vectors for non-centring part
  pgvec <- c("","","")
  i <- 1
  while (substr(ssym,i,i) != ",") i <- i+1
  pgvec[1]<- paste(pgvec[1],substr(ssym,7,i-1),sep="")
  ssym <- substr(ssym,i+1,nchar(ssym))
  i <- 1
  while (substr(ssym,i,i) != ",") i <- i+1
  pgvec[2] <- paste(pgvec[2],substr(ssym,1,i-1),sep="")
  ssym <- substr(ssym,i+1,nchar(ssym))
  pgvec[3] <- paste(pgvec[3],ssym,sep="")
  for (iscen in 1:length(cstring))
  {
   cnvec <- c("","","")
   scen <- cstring[iscen]
   if (iscen == 1)
   {
    final_symm <- c(final_symm,paste("symm ",pgvec[1],",  ",pgvec[2],",  ",pgvec[3],sep=""))
   }
   if (iscen > 1)
   {
    i <- 1
    while (substr(scen,i,i) != ",") i <- i+1
    sa <- substr(scen,7,i-1)
    scen <- substr(scen,i+1,nchar(scen))
    i <- 1
    while (substr(scen,i,i) != ",") i <- i+1
    sb <- substr(scen,1,i-1)
    scen <- substr(scen,i+1,nchar(scen))
    sc <- scen
    cnvec <- c("","","")
    if (nchar(sa) > 2 & substr(sa,2,2) == "+" | substr(sa,2,2) == "-") cnvec[1] <- substr(sa,2,nchar(sa))
    if (nchar(sa) > 2 & substr(sa,3,3) == "+" | substr(sa,3,3) == "-") cnvec[1] <- substr(sa,3,nchar(sa))
    if (nchar(sb) > 2 & substr(sb,2,2) == "+" | substr(sb,2,2) == "-") cnvec[2] <- substr(sb,2,nchar(sb))
    if (nchar(sb) > 2 & substr(sb,3,3) == "+" | substr(sb,3,3) == "-") cnvec[2] <- substr(sb,3,nchar(sb))
    if (nchar(sc) > 2 & substr(sc,2,2) == "+" | substr(sc,2,2) == "-") cnvec[3] <- substr(sc,2,nchar(sc))
    if (nchar(sc) > 2 & substr(sc,3,3) == "+" | substr(sc,3,3) == "-") cnvec[3] <- substr(sc,3,nchar(sc))
    final_symm <- c(final_symm,paste("symm ",pgvec[1],cnvec[1],",  ",pgvec[2],cnvec[2],",  ",pgvec[3],cnvec[3],sep=""))
   }
  }
 } 
 final_symm <- toupper(final_symm)
 for (i in 1:length(final_symm)) while(nchar(final_symm[i]) != 80) final_symm[i] <- paste(final_symm[i]," ",sep="")

 return(final_symm)
}

## Function to find correct Herman-Maguin spelling starting from commonly used strings
#
.findHM <- function(sym_xHM)
{
 # Possible inputs  
 if (sym_xHM == "P 2") sym_xHM <- "P 1 2 1"
 if (sym_xHM == "P 21") sym_xHM <- "P 1 21 1"
 if (sym_xHM == "C 2") sym_xHM <- "C 1 2 1"
 if (sym_xHM == "P m") sym_xHM <- "P 1 m 1"
 if (sym_xHM == "P c") sym_xHM <- "P 1 c 1"
 if (sym_xHM == "C m") sym_xHM <- "C 1 m 1"
 if (sym_xHM == "C c") sym_xHM <- "C 1 c 1"
 if (sym_xHM == "P 2/m") sym_xHM <- "P 1 2/m 1"
 if (sym_xHM == "P 21/m") sym_xHM <- "P 1 21/m 1"
 if (sym_xHM == "C 2/m") sym_xHM <- "C 1 2/m 1"
 if (sym_xHM == "P 2/c") sym_xHM <- "P 1 2/c 1"
 if (sym_xHM == "P 21/c") sym_xHM <- "P 1 21/c 1"
 if (sym_xHM == "C 2/c") sym_xHM <- "C 1 2/c 1"
 if (sym_xHM == "P n n n") sym_xHM <- "P n n n :1"
 if (sym_xHM == "P b a n") sym_xHM <- "P b a n :1"
 if (sym_xHM == "P m m n") sym_xHM <- "P m m n :1"
 if (sym_xHM == "C c c a") sym_xHM <- "C c c a :1"
 if (sym_xHM == "F d d d") sym_xHM <- "F d d d :1"
 if (sym_xHM == "P 4/n") sym_xHM <- "P 4/n :1"
 if (sym_xHM == "P 42/n") sym_xHM <- "P 42/n :1"
 if (sym_xHM == "P 41/a") sym_xHM <- "P 41/a :1"
 if (sym_xHM == "P 4/n b m") sym_xHM <- "P 4/n b m :1"
 if (sym_xHM == "P 4/n n c") sym_xHM <- "P 4/n n c :1"
 if (sym_xHM == "P 4/n m m") sym_xHM <- "P 4/n m m :1"
 if (sym_xHM == "P 4/n c c") sym_xHM <- "P 4/n c c :1"
 if (sym_xHM == "P 42/n b c") sym_xHM <- "P 42/n b c :1"
 if (sym_xHM == "P 42/n n m") sym_xHM <- "P 42/n n m :1"
 if (sym_xHM == "P 42/n m c") sym_xHM <- "P 42/n m c :1"
 if (sym_xHM == "P 42/n c m") sym_xHM <- "P 42/n c m :1"
 if (sym_xHM == "I 41/a m d") sym_xHM <- "I 41/a m d :1"
 if (sym_xHM == "I 41/a c d") sym_xHM <- "I 41/a c d :1"
 if (sym_xHM == "H 3") sym_xHM <- "R 3 :H"
 if (sym_xHM == "H -3") sym_xHM <- "R -3 :H"
 if (sym_xHM == "H 3 2") sym_xHM <- "R 3 2 :H"
 if (sym_xHM == "H 3 m") sym_xHM <- "R 3 m :H"
 if (sym_xHM == "H 3 c") sym_xHM <- "R 3 c :H"
 if (sym_xHM == "H -3 m") sym_xHM <- "R -3 m :H"
 if (sym_xHM == "H -3 c") sym_xHM <- "R -3 c :H"
 if (sym_xHM == "R 3") sym_xHM <- "R 3 :R"
 if (sym_xHM == "R -3") sym_xHM <- "R -3 :R"
 if (sym_xHM == "R 3 2") sym_xHM <- "R 3 2 :R"
 if (sym_xHM == "R 3 m") sym_xHM <- "R 3 m :R"
 if (sym_xHM == "R 3 c") sym_xHM <- "R 3 c :R"
 if (sym_xHM == "R -3 m") sym_xHM <- "R -3 m :R"
 if (sym_xHM == "R -3 c") sym_xHM <- "R -3 c :R"
 if (sym_xHM == "P n -3") sym_xHM <- "P n -3 :1"
 if (sym_xHM == "F d -3") sym_xHM <- "F d -3 :1"
 if (sym_xHM == "P n -3 n") sym_xHM <- "P n -3 n :1"
 if (sym_xHM == "F d -3 m") sym_xHM <- "F d -3 m :1"
 if (sym_xHM == "F d -3 c") sym_xHM <- "F d -3 c :1"

 # Returned, corrected xHM symbol
 return(sym_xHM)
}

## Get symmetry number and setting given object of class Symmetry
#
.getSymmetryNumber <- function(object)
{
 # Nothing to return if Symmetry is empty
 if (length(object@sym_xHM) == 0) return(NULL)
 
 # Get symmetry number first
 tmp <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
 if (!tmp$ans) return(tmp$msg)
 sgn <- tmp$msg

 # Count how many settings
 nsett <- 0
 while (tmp$ans)
 {
  nsett <- nsett+1
  tmp <- .translate_SG(value=sgn,SG_in="number",SG_out="xHM",setting=nsett)
 }
 nsett <- nsett-1

 # To finish, check which xHM symbol coincides with the one related to number and setting
 setting <- NULL
 for (i in 1:nsett)
 {
  tmp <- .translate_SG(value=sgn,SG_in="number",SG_out="xHM",setting=i)
  if (tmp$msg == object@sym_xHM) setting <- i                    
 }

 return(c(sgn,setting))
}

## Using space group number and setting returns code related to unit cell parameters constraints
## It is used only for monoclinic groups
#
.getMonoclinicConstraints <- function(sgn)
{
 # sgn is a vector of length 2: c(space group number,setting)

 # Monoclinic constraints
 if (sgn[1] %in% 3:15)
 {
  if (sgn[1] == 3)
  {
   if (sgn[2] == 1) idx <- 5
   if (sgn[2] == 2) idx <- 6
  }
  if (sgn[1] == 4)
  {
   if (sgn[2] == 1) idx <- 5
   if (sgn[2] == 2) idx <- 6
  }
  if (sgn[1] == 5)
  {
   if (sgn[2] %in% c(1:3,10,11)) idx <- 5
   if (sgn[2] %in% 4:6) idx <- 6
   if (sgn[2] %in% 7:9) idx <- 4
  }
  if (sgn[1] == 6)
  {
   if (sgn[2] == 1) idx <- 5
   if (sgn[2] == 2) idx <- 6
   if (sgn[2] == 3) idx <- 4
  }
  if (sgn[1] == 7)
  {
   if (sgn[2] %in% 1:3) idx <- 5
   if (sgn[2] %in% 4:6) idx <- 6
   if (sgn[2] %in% 7:9) idx <- 4
  }
  if (sgn[1] == 8)
  {
   if (sgn[2] %in% 1:3) idx <- 5
   if (sgn[2] %in% c(4:6,10)) idx <- 6
   if (sgn[2] %in% 7:9) idx <- 4
  }
  if (sgn[1] == 9)
  {
   if (sgn[2] %in% 1:6) idx <- 5
   if (sgn[2] %in% 7:12) idx <- 6
   if (sgn[2] %in% 13:18) idx <- 4
  } 
  if (sgn[1] == 10)
  {
   if (sgn[2] == 1) idx <- 5
   if (sgn[2] == 2) idx <- 6
   if (sgn[2] == 3) idx <- 4
  }
  if (sgn[1] == 11)
  {
   if (sgn[2] == 1) idx <- 5
   if (sgn[2] == 2) idx <- 6
   if (sgn[2] == 3) idx <- 4
  }
  if (sgn[1] == 12)
  {
   if (sgn[2] %in% 1:3) idx <- 5
   if (sgn[2] %in% 4:6) idx <- 6
   if (sgn[2] %in% 7:9) idx <- 4
  }
  if (sgn[1] == 13)
  {
   if (sgn[2] %in% 1:3) idx <- 5
   if (sgn[2] %in% 4:6) idx <- 6
   if (sgn[2] %in% 7:9) idx <- 4
  }
  if (sgn[1] == 14)
  {
   if (sgn[2] %in% 1:3) idx <- 5
   if (sgn[2] %in% 4:6) idx <- 6
   if (sgn[2] %in% 7:9) idx <- 4
  }
  if (sgn[1] == 15)
  {
   if (sgn[2] %in% 1:6) idx <- 5
   if (sgn[2] %in% 7:12) idx <- 6
   if (sgn[2] %in% 13:18) idx <- 4
  }
 }

 # Trigonal constraints

 return(idx)
}


############# MOSTLY RELATED TO UNIT CELL, LATTICE AND FRACTIONAL/ORTHOGONAL COORDINATES (OLD lattice_dot.R) #############################

## Input: cell parameters. Output: values useful to all crystallographic calculations
#
.dcl <- function(a,b,c,aa,bb,gg)
{
# Input: cell parameters. Output: values useful to all crystallographic
# calculations. These are: 
# 1) sa = sin(alpha), sb = sin(beta), sc = sin(gamma)
# 2) ca = cos(alpha), cb = cos(beta). cc = cos(gamma)
# 3) sides of reciprocal cell: ar = a*, br = b*, cr = c*
# 4) sines of angles of reciprocal cell: sar = sin(alpha*), sbr = sin(beta*), scr = sin(gamma*)
# 5) cosines of angles of reciprocal cell: car = cos(alpha*), cbr = cos(beta*), ccr = cos(gamma*)
# 6) Volume of unit cell: V
 aa <- aa*pi/180
 bb <- bb*pi/180
 gg <- gg*pi/180
 sa <- sin(aa)
 sb <- sin(bb)
 sc <- sin(gg)
 ca <- cos(aa)
 cb <- cos(bb)
 cc <- cos(gg)

 # To avoid NaN generated by rounding off errors, use for cell-derived quantities formulas
 # derived previously by computationa crystallographers
 sang <- 0.5*(aa+bb+gg)
 V2 <- sqrt(sin(sang-aa)*sin(sang-bb)*sin(sang-gg)*sin(sang))
 V <- 2*a*b*c*V2
 ar <- b*c*sa/V
 br <- a*c*sb/V
 cr <- a*b*sc/V
 car <- (cb*cc-ca)/(sb*sc)
 cbr <- (ca*cc-cb)/(sa*sc)
 ccr <- (ca*cb-cc)/(sa*sb)
 sar <- sqrt(1-car*car)
 sbr <- sqrt(1-cbr*cbr)
 scr <- sqrt(1-ccr*ccr)
 l <- c(sa,sb,sc,ca,cb,cc,ar,br,cr,sar,sbr,scr,car,cbr,ccr,V)
 names(l) <- c("SIN_ALPHA","SIN_BETA","SIN_GAMMA","COS_ALPHA","COS_BETA","COS_GAMMA","A*","B*","C*","SIN_ALPHA*","SIN_BETA*","SIN_GAMMA*",
               "COS_ALPHA*","COS_BETA*","COS_GAMMA*","V")

 return(l)
}

## Input: space group number. Output: crystal system
#
.crystal_system <- function(gn)
{
 # Given the space group number (gn), returns the crystal system
 # bs = 1 TRICLINIC
 # bs = 2 MONOCLINIC
 # bs = 3 ORTHOROMBIC
 # bs = 4 TETRAGONAL
 # bs = 5 CUBIC
 # bs = 6 HEXAGONAL
 # bs = 7 TRIGONAL
 if (gn >=   1 & gn <=   2) bs <- 1
 if (gn >=   3 & gn <=  15) bs <- 2
 if (gn >=  16 & gn <=  74) bs <- 3
 if (gn >=  75 & gn <= 142) bs <- 4
 if (gn >= 143 & gn <= 167) bs <- 7
 if (gn >= 168 & gn <= 194) bs <- 6
 if (gn >= 195 & gn <= 230) bs <- 5
 bs_name <- c("TRICLINIC","MONOCLINIC","ORTHOROMBIC","TETRAGONAL","CUBIC","HEXAGONAL","TRIGONAL")

 return(bs_name[bs])
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (First choice in Giacovazzo's book)
.triclinic_to_orthogonal_01 <- function(a,b,c,aa,bb,cc)
{
 lp <- .dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(a,0,0,b*lp[6],b*lp[3],0,c*lp[5],-c*lp[2]*lp[13],1/lp[9]),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (MOSFLM choice, second choice in Giacovazzo's book)
.triclinic_to_orthogonal_02 <- function(a,b,c,aa,bb,cc)
{
 lp <- .dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(1/lp[7],-lp[15]/(lp[7]*lp[12]),a*lp[5],0,1/(lp[8]*lp[12]),b*lp[4],0,0,c),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

#
# Transform fractional into orthonormal coordinates
.frac_to_orth <- function(xyzf,a,b,c,aa,bb,cc,ochoice=1)
{
 # a,b,c = cell sides in angstroms; aa,bb,cc = cell angles in degrees
 # The parameter ochoice controls which convention is being used to
 # collocate the cell in an orthonormal cartesian frame.
 # ochoice = 1: X axis along a; Y axis normal to a, in the (a,b) plane;
 #              Z axis normal to X and Y (and therefore parallel to
 #              c*).
 # ochoice = 2: this is also called "Cambridge setting". The X axis is
 #              along a*; the Y axis lies in the (a*,b*) plane; the Z
 #              axis is, consequently, along c. 
 # xyzf is a vector, matrix or data frame of fractional crystal coordinates.
 # If matrix or data frame it needs to have 3 columns. 
 # This function returns a data frame with 3 columns, the cartesian orthonormal coordinates.

 # Check xyzf object is either vector, matrix or data.frame
 if (!is.vector(xyzf) & !is.matrix(xyzf) & !is.data.frame(xyzf)) stop("Input object is not of a valid type (vector, matrix or data.frame)")

 # If xyzf is a vector, turn it into a matrix
 if (is.vector(xyzf)) dim(xyzf) <- c(1,length(xyzf))

 # If xyzf has less or more than 3 columns stop
 if (dim(xyzf)[2] != 3) stop("Input object has more or less than 3 columns")

 # If xyzf is a data.frame turn it into a matrix
 if (is.data.frame(xyzf)) xyzf <- as.matrix(xyzf)

 # Orthogonalization matrix
 if (ochoice == 1) M_1 <- .triclinic_to_orthogonal_01(a,b,c,aa,bb,cc)
 if (ochoice == 2) M_1 <- .triclinic_to_orthogonal_02(a,b,c,aa,bb,cc)

 # Transformed coordinates
 #x <- av[1]*xf+bv[1]*yf+cv[1]*zf
 #y <- av[2]*xf+bv[2]*yf+cv[2]*zf
 #z <- av[3]*xf+bv[3]*yf+cv[3]*zf
 xyz <- xyzf%*%M_1 

 # Turn matrix back into a data frame
 xyz <- as.data.frame(xyz)
 colnames(xyz) <- c("x","y","z")

 return(xyz)
}

#
# Transforms orthogonal into fractional coordinates
.orth_to_frac <- function(xyz,a,b,c,aa,bb,cc,ochoice=1)
{
 # Does the inverse job of "frac_to_orth"

 # Check xyz object is either vector, matrix or data.frame
 if (!is.vector(xyz) & !is.matrix(xyz) & !is.data.frame(xyz)) stop("Input object is not of a valid type (vector, matrix or data.frame)")

 # If xyz is a vector, turn it into a matrix
 if (is.vector(xyz)) dim(xyz) <- c(1,length(xyz))

 # If xyz has less or more than 3 columns stop
 if (dim(xyz)[2] != 3) stop("Input object has more or less than 3 columns")

 # If xyz is a data.frame turn it into a matrix
 if (is.data.frame(xyz)) xyz <- as.matrix(xyz)

  # Orthogonalization matrix
 if (ochoice == 1) M_1 <- .triclinic_to_orthogonal_01(a,b,c,aa,bb,cc)
 if (ochoice == 2) M_1 <- .triclinic_to_orthogonal_02(a,b,c,aa,bb,cc)

 # Inverse matrix
 M <- solve(M_1)

 # Transform coordinates
 xyzf <- xyz%*%M

 # Turn matrix back into a data frame
 xyzf <- as.data.frame(xyzf)
 colnames(xyzf) <- c("xf","yf","zf")

 return(xyzf)
}

.d_hkl <- function(h,k,l,a,b,c,aa,bb,cc)
{
 # Given Miller indices and cell parameters, this function returns
 # resolution corresponding to the specific Miller indices.

 aa <- aa*pi/180
 bb <- bb*pi/180
 cc <- cc*pi/180
 top <- 1-(cos(aa))^2-(cos(bb))^2-(cos(cc))^2+2*cos(aa)*cos(bb)*cos(cc)
 b1 <- h^2*(sin(aa))^2/a^2
 b2 <- k^2*(sin(bb))^2/b^2
 b3 <- l^2*(sin(cc))^2/c^2
 b4 <- 2*h*k*(cos(aa)*cos(bb)-cos(cc))/(a*b)
 b5 <- 2*h*l*(cos(aa)*cos(cc)-cos(bb))/(a*c)
 b6 <- 2*k*l*(cos(bb)*cos(cc)-cos(aa))/(b*c)
 d2 <- top/(b1+b2+b3+b4+b5+b6)
 return(sqrt(d2))
}

## Extract information from Lattice object. Warning: as this is not meant for public use, the function
## does not test for object to be of class Lattice.
#
.extractLatticeStuff <- function(object)
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


############################# MOSTLY RELATED MTZ AND REFLECTIONS (OLD reflections_dot.R) ##########################################

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

## Used to check beginning of an MTZ file
.try_read_MTZ <- function(filename)
{
 # Create a connection to binary file
 f <- file(filename, open="rb")

 # Reads initial record
 irdata <- .readIR(f)
 close(f)
 errF <- irdata[[1]]

 # Return FALSE if errF != 0
 rval <- TRUE
 if (errF != 0) rval <- FALSE

 return(rval)
}

############################# MOSTLY RELATED PDB AND ATOMS (OLD atom_dot.R) ##########################################

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

## Used to verify a given file is a PDB file
.try_read_PDB <- function(filename)
{
 pdb0 <- scan(filename,what="character",sep="\n",quiet=TRUE)
 value <- "CRYST1"
 suppressWarnings(bsg <- grep(value,pdb0,fixed=TRUE))
 if (length(bsg) == 0)
 {
  rval <- FALSE
 }
 if (length(bsg) != 0)
 {
  keychar <- substr(pdb0[bsg[1]],1,1)
  if (keychar == "C") rval <- TRUE
  if (keychar != "C") rval <- FALSE
 }

 return(rval)
}


############################# MOSTLY RELATED TO GENERIC STUFF (OLD generic_dot.R) ##########################################
# For use in conjunction with S4-classes based modules in cRy.
# The word "dot" has been adopted because all functions here start with a ".", so
# to be invisible.

#
# Input word "wordin" is stripped of blanks before, after or before and after the word itself.
# side can be "left","right", "both".
.strip_blanks <- function(wordin,side="both")
{
 if (side != "left" & side != "right" & side != "both") stop("Input variable side can only be 'left', 'right' or 'both'")

 if (side == "left" | side == "both")
 {
  nb <- 1
  while(substr(wordin,nb,nb) == " ") nb <- nb+1
  wordin <- substr(wordin,nb,nchar(wordin))
 }
  
 if (side == "right" | side == "both")
 {
  nb <- nchar(wordin)
  while(substr(wordin,nb,nb) == " ") nb <- nb-1
  wordin <- substr(wordin,1,nb)
 }

 return(wordin)
}
#
# Given 3 Euler angles returns 3X3 rotation matrix
.euler_to_matrix <- function(alpha,beta,gamma)
{
 # Input angles are read in degrees
 aa <- alpha*pi/180
 bb <- beta*pi/180
 cc <- gamma*pi/180
 R1 <- matrix(c(cos(cc),sin(cc),0,-sin(cc),cos(cc),0,0,0,1),nrow=3,ncol=3)
 R2 <- matrix(c(1,0,0,0,cos(bb),sin(bb),0,-sin(bb),cos(bb)),nrow=3,ncol=3)
 R3 <- matrix(c(cos(aa),sin(aa),0,-sin(aa),cos(aa),0,0,0,1),nrow=3,ncol=3)
 Rfinal <- R3%*%R2%*%R1

 return(Rfinal)
}
#
# Given 3 polar angles returns 3X3 rotation matrix
.polar_to_matrix <- function(chi,psi,fi)
{
 # Input angles are in degrees
 chi <- chi*pi/180
 psi <- psi*pi/180
 fi <- fi*pi/180

 # Matrix components
 o11 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(cos(fi))^2
 o12 <- -sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o13 <- -sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o21 <- sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o22 <- cos(chi)+(1-cos(chi))*(cos(psi))^2
 o23 <- sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o31 <- sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o32 <- -sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o33 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(sin(fi))^2

 # Matrix
 oline <- c(o11,o12,o13,o21,o22,o23,o31,o32,o33)
 Om <- matrix(oline,nrow=3,ncol=3,byrow=TRUE)

 return(Om)
}
#
# Cross product u X v in 3D. Returns a 3 components vector
"%X%" <- function(u,v)
{
 x <- u[2]*v[3]-u[3]*v[2]
 y <- u[3]*v[1]-u[1]*v[3]
 z <- u[1]*v[2]-u[2]*v[1]
 return(c(x,y,z))
}

## Find file type from extension (so far: MTZ, PDB)
