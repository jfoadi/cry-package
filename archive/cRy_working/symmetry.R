# Functions to handle symmetry ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### START
#
#####################################################################
#
# James Foadi and David Waterman
#
# Diamond Light Source & Imperial College London
#
# September 2009
#
#####################################################################

# All symmetry starts with CCP4 file "syminfo.lib"
symmetry_file <- SYMMETRY_file

# Functions
source(paste(CRY_HOME,"lattice.R",sep="/"))

#
syminfo_to_op_xyz_list <- function(SG)
{
 # Input is spacegroup symbol in xHM format.
 # Output is a list of 2 vector of strings,
 # The first string describes point group and translation; the
 # second string describes cell centering 

 # Extract full symmetry information 
 data <- extract_symmetry_info(SG)

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
translate_SG <- function(value,SG_in="number",SG_out="xHM",setting=1)
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

 # If input is "number", turn into a string
 if (SG_in == "number") value <- paste("number ",value)
 if (SG_in == "ccp4") value <- paste("symbol ccp4",value)

 # Complete string for non-number cases
 if (SG_in == "Hall") value <- paste("symbol Hall '",value,"'",sep="")
 if (SG_in == "xHM") value <- paste("symbol xHM  '",value,"'",sep="")
 if (SG_in == "old") value <- paste("symbol old  '",value,"'",sep="")

 # Read whole content of "syminfo.lib"
 syminfo <- scan(symmetry_file,what="character",sep="\n",quiet=TRUE)
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
  print("   2) for this space group there are not that many settings")
  print("   3) wrong SG_in/SG_out combination with value")
  stop()
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
  data <- op_xyz_to_matrix(op_xyz_list[[1]][i])
  matrix_list[i] <- list(data[[1]])
  vector_list[i] <- list(data[[2]])
 }

 # Centering vectors
 for (i in 1:length(op_xyz_list[[2]]))
 {
  data <- op_xyz_to_matrix(op_xyz_list[[2]][i])
  centering_list[i] <- list(data[[2]])
 }

 overall_list <- list(PG=matrix_list,T=vector_list,C=centering_list)

 return(overall_list)
}

#
op_xyz_to_matrix <- function(op_xyz)
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
.prepare_op_xyz <- function(SYMMstring)
{
 # Turn a string of the "SYMM    " format into one interpretable for function op_xyz_to_matrix 

 vstmp <- strsplit(SYMMstring,",",fixed=TRUE)[[1]]
 vstmp <- unname(sapply(vstmp,noBlanks))
 op_xyz <- tolower(paste(vstmp[1],vstmp[2],vstmp[3],sep=","))

 return(op_xyz)
}

#
syminfo_to_matrix_list <- function(SG)
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
 return(.op_xyz_list_to_matrix_list(syminfo_to_op_xyz_list(SG)))
}

#
extract_symmetry_info <- function(SG)
{
 # Input is spacegroup symbol in xHM format.
 # Output is a list with all symmetry information
 # for the specific SG group, as contained in syminfo.lib.

 # Read whole content of "syminfo.lib"
 syminfo <- scan(symmetry_file,what="character",sep="\n",quiet=TRUE)
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
equiv_hkl <- function(hkl,SG)
{
 # Given a point in reciprocal space with Miller indices h, k and l
 # and a space group with xHM symbol SG, produces a data frame of equivalent
 # points with Miller indices grouped in columns named "h","k","l","pg_sym", "orig_hkl".
 # The last columns indicate which symmetry operator was responsible for the specific
 # equivalent reflection and to which set of original Miller indices is the specific
 # index related.
 # Input: hkl can be a vector, matrix or data.frame; if a matrix or a data.frame,
 #        it needs to have 3 columns
 # Input: SG is the space group symbol in xHM format
 # Output: a data frame with 5 columns: h k l pg_sym orig_hkl

 # Check hkl object is either vector, matrix or data.frame
 if (!is.vector(hkl) & !is.matrix(hkl) & !is.data.frame(hkl)) stop("Input object is not of a valid type (vector, matrix or data.frame)")

 # If hkl is a vector, turn it into a matrix
 if (is.vector(hkl)) dim(hkl) <- c(1,length(hkl))

 # If hkl has less or more than 3 columns stop
 if (dim(hkl)[2] != 3) stop("Input object has more or less than 3 columns")

 # If hkl is a data.frame turn it into a matrix
 if (is.data.frame(hkl)) hkl <- as.matrix(hkl)

 # Name row in all cases
 rownames(hkl) <- 1:dim(hkl)[1]

 # Get list of point group matrices
 PG <- syminfo_to_matrix_list(SG)[[1]]

 # Initialise data frame of equivalent reciprocal-lattice points
 hkl_out <- data.frame(h=NA,k=NA,l=NA,pg_sym=NA,orig_hkl=NA)

 # Cycle through point group matrices
 for (i in 1:length(PG))
 {
  # Transpose vector of equivalent Miller indices
  tmp_hkl <- t(mapply((function(h,k,l,PG,i) {r <- c(h,k,l)%*%PG[[i]]}),hkl[,1],hkl[,2],hkl[,3],MoreArgs=list(PG=PG,i=i)))

  # Fill temporary data frame
  tmp_data <- data.frame(h=tmp_hkl[,1],k=tmp_hkl[,2],l=tmp_hkl[,3],pg_sym=i,orig_hkl=rownames(hkl))

  # Update main data frame
  hkl_out <- rbind(hkl_out,tmp_data)
 }
 rm(tmp_hkl,tmp_data)

 # Get rid of first NULL row
 hkl_out <- na.omit(hkl_out)

 # Get rid of duplicated reflections
 tmp <- unique(hkl_out[,1:3])
 hkl_out <- hkl_out[rownames(tmp),]
 rm(tmp)

 return(hkl_out)
}

#
equiv_xyzf <- function(xyzf,SG,in_cell=FALSE)
{
 # Given fractional coordinates for a 3D point, returns all equivalent points
 # according to symmetry operators of space group with xHM symbol SG.
 # If "in_cell" is TRUE, all points outside cell are translated back into cell.
 # Input: xyzf can be a vector, matrix or data.frame; if a matrix or a data.frame,
 #        it needs to have 3 columns
 # Input: SG xHM format for space group
 # Input: in_cell, TRUE or FALSE
 # Output: data frame with 6 columns: xf yf zf sym c_sym orig_xyzf. "sym" indicates which
 #         operator produced the specific "xf", "yf", "zf". "c_sym" specifies the centering.
 #         "orig_xyzf" indicates to which initial point the produced point is related.

 # Check xyzf object is either vector, matrix or data.frame
 if (!is.vector(xyzf) & !is.matrix(xyzf) & !is.data.frame(xyzf)) stop("Input object is not of a valid type (vector, matrix or data.frame)")

 # If xyzf is a vector, turn it into a matrix
 if (is.vector(xyzf)) dim(xyzf) <- c(1,length(xyzf))

 # If xyzf has less or more than 3 columns stop
 if (dim(xyzf)[2] != 3) stop("Input object has more or less than 3 columns")

 # If xyzf is a data.frame turn it into a matrix
 if (is.data.frame(xyzf)) xyzf <- as.matrix(xyzf)

 # Name row in all cases
 rownames(xyzf) <- 1:dim(xyzf)[1]

 # Get list of space group matrices
 MT <- syminfo_to_matrix_list(SG)

 # Number of non-cell-centering operators
 nsym <- length(MT$PG)
 ncent <- length(MT$C)
 
 # Create dataframe
 eq_xyzf <- data.frame(xf=NA,yf=NA,zf=NA,sym=NA,c_sym=NA,orig_xyzf=NA)

 # Loop over elements of the list. External loop for centering, internal for symmetry
 for (j in 1:ncent)
 {
  for (i in 1:nsym)
  {
   new_xyzf <- t(mapply((function(xf,yf,zf,MT,i,j) {r <- c(xf,yf,zf)%*%t(MT$PG[[i]])+MT$T[[i]]+MT$C[[j]]}),xyzf[,1],xyzf[,2],xyzf[,3],MoreArgs=list(MT=MT,i=i,j=j)))
   #new_xyzf <- xyzf%*%t(MT$PG[[i]])+MT$T[[i]]+MT$C[[j]]
   eq_xyzf <- rbind(eq_xyzf,data.frame(xf=new_xyzf[,1],yf=new_xyzf[,2],zf=new_xyzf[,3],sym=i,c_sym=j,orig_xyzf=rownames(xyzf)))
  }
 }

 # Get rid of NA
 eq_xyzf <- na.omit(eq_xyzf)

 # Put points inside unit cell, if required
 if (in_cell) eq_xyzf[,1:3] <- eq_xyzf[,1:3]%%1

 # Get rid of duplicated reflections
 tmp <- unique(eq_xyzf[,1:3])
 eq_xyzf <- eq_xyzf[rownames(tmp),]
 rm(tmp)

 return(eq_xyzf)
}

#
asu_ccp4_limits <- function(SG)
{
 # Given space group xHM symbol, returns asymmetric unit limits along
 # x, y and z (according to ccp4 convention).

 # Extract string containing "mapasu_ccp4"
 mapasu_ccp4 <- extract_symmetry_info(SG)$MAPASU_CCP4
 tmp <- substr(mapasu_ccp4,13,nchar(mapasu_ccp4))
 tmp <- unlist(strsplit(tmp," "))
 tmp <- unlist(strsplit(tmp,fixed=TRUE,";"))
 
 # x limits (xmin is always 0)
 xmin <- 0
 if (substr(tmp[1],nchar(tmp[1])-1,nchar(tmp[1])-1) == "/")
 {
  xtop <- as.numeric(substr(tmp[1],nchar(tmp[1])-2,nchar(tmp[1])-2))
  xbot <- as.numeric(substr(tmp[1],nchar(tmp[1]),nchar(tmp[1])))
  xmax <- xtop/xbot
 }
 if (substr(tmp[1],nchar(tmp[1])-1,nchar(tmp[1])-1) != "/")
 {
  xmax <- as.numeric(substr(tmp[1],nchar(tmp[1]),nchar(tmp[1])))
 }

 # y limits (ymin is always 0)
 ymin <- 0
 if (substr(tmp[2],nchar(tmp[2])-1,nchar(tmp[2])-1) == "/")
 {
  ytop <- as.numeric(substr(tmp[2],nchar(tmp[2])-2,nchar(tmp[2])-2))
  ybot <- as.numeric(substr(tmp[2],nchar(tmp[2]),nchar(tmp[2])))
  ymax <- ytop/ybot
 }
 if (substr(tmp[2],nchar(tmp[2])-1,nchar(tmp[2])-1) != "/")
 {
  ymax <- as.numeric(substr(tmp[2],nchar(tmp[2]),nchar(tmp[2])))
 }

 # z limits (zmin is always 0)
 zmin <- 0
 if (substr(tmp[3],nchar(tmp[3])-1,nchar(tmp[3])-1) == "/")
 {
  ztop <- as.numeric(substr(tmp[3],nchar(tmp[3])-2,nchar(tmp[3])-2))
  zbot <- as.numeric(substr(tmp[3],nchar(tmp[3]),nchar(tmp[3])))
  zmax <- ztop/zbot
 }
 if (substr(tmp[3],nchar(tmp[3])-1,nchar(tmp[3])-1) != "/")
 {
  zmax <- as.numeric(substr(tmp[3],nchar(tmp[3]),nchar(tmp[3])))
 }
 limits <- matrix(c(xmin,xmax,ymin,ymax,zmin,zmax),ncol=2,byrow=TRUE)

 return(limits)
}
# Functions to handle symmetry ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ### ### SYMMETRY ###   END
