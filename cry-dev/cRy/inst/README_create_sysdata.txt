!!! This information is for cRy developers only !!!

The file "sysdata.rda" contained under the "R/" directory, includes vital
information for the whole cRy. This information is hidden to the user.

If at any time the file "sysdata.rda" is lost it can be recreated starting
from the 3 files,

atomsf.lib
syminfo.lib
elements_list.dat
amino_list.dat

simply by starting a pristine R session under "inst/" and typing the following:

 syminfo <- scan("syminfo.lib",what="character",sep="\n",quiet=TRUE)
 .ATOMS_data.frame <- read.table("elements_list.dat",header=TRUE,as.is=1:2)
 .ALPHABET <- c(LETTERS," ")
 .AMINO_ACIDS_data.frame <- read.table("amino_list.dat")
 save.image("sysdata.rda")

Information from atomsf.lib is temporarily unloaded. We will add it to "sysdata.rda" once
it is included in the source code.
