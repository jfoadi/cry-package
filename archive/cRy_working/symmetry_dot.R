# For use in conjunction with S4-classes based modules in cRy.
# The word "dot" has been adopted because all functions here start with a ".", so
# to be invisible.

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
  if (length(idx) == 0) stop("There is not a space group associated with the input number",call.=FALSE)
 }

 # If input is "number", turn into a string
 if (SG_in == "number") value <- paste("number ",value)
 if (SG_in == "ccp4") value <- paste("symbol ccp4",value)

 # Complete string for non-number cases
 if (SG_in == "Hall") value <- paste("symbol Hall '",value,"'",sep="")
 if (SG_in == "xHM") value <- paste("symbol xHM  '",value,"'",sep="")
 if (SG_in == "old") value <- paste("symbol old  '",value,"'",sep="")

 # Read whole content of "syminfo.lib"
 syminfo <- scan(.SYMMETRY_file,what="character",sep="\n",quiet=TRUE)
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
  print("Something wrong in your input:")
  print("   1) the symbol or number input for this space group does not exist")
  print("   2) if your inpur was a number, perhaps for this space group there are not that many settings")
  #print("   3) wrong SG_in/SG_out combination with value")
  stop(call.=FALSE)
 }
 bsg <- bsg[setting]

 bini <- grep("begin_spacegroup",syminfo)
 prima <- bini[length(bini[bini < bsg])]
 dopo  <- bini[bini > bsg][1]-1   # Add 1 for those cases when prima = dopo
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

 return(translated_value)
}

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

#
.extract_symmetry_info <- function(SG)
{
 # Input is spacegroup symbol in xHM format.
 # Output is a list with all symmetry information
 # for the specific SG group, as contained in syminfo.lib.

 # Read whole content of "syminfo.lib"
 syminfo <- scan(.SYMMETRY_file,what="character",sep="\n",quiet=TRUE)
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

# 
.full_symm_strings <- function(SG)
{
 # Returns all symmetry operators in string format (the kind SYMM)

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
