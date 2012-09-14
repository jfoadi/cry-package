

#
create_SYMM_line <- function(SG)
{
 # Given space group number calls some symmetry functions to return operators in string forms and
 # prepare them in capital SYMM format

 # Load symmetry module
 source(paste(CRY_HOME,"symmetry.R"))

 # Find xHM symbol
 xHM <- translate_SG(SG)

 # Create space after ","
 lista <- syminfo_to_op_xyz_list(xHM)[[1]]
 new_lista <- c()
 for (i in 1:length(lista))
 {
  svector <- strsplit(lista[i],",")[[1]]
  stmp <- ""
  for (j in 1:length(svector))
  {
   if (j != length(svector)) stmp <- paste(stmp,svector[j],",  ",sep="")
   if (j == length(svector)) stmp <- paste(stmp,svector[j],sep="")
 }
  new_lista[i] <- stmp
 }

 # Find list of operators in lowercase
 op_list <- toupper(new_lista)

 return(op_list)
}
