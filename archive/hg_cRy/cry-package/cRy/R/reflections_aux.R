sysabs  <- function(hkl,sym_number,setting)
{
 # Input: hkl is a m X 3 matrix whose lines are Miller indices
 #        sym_number is the symmetry number in the International Tables for Crystallography
 #        setting is the setting as in the "syminfo.lib" file
 # Output is a vector containing indices of Miller indices obeying systematic absences


 # First add a fourth column with 0 to matrix
 nrefs <- nrow(hkl)
 hkl <- cbind(hkl,matrix(rep(0,times=nrefs),nrow=nrefs))
 colnames(hkl)[4] <- "FLAG"

 # Long list for all space groups and settings

 #########################################################################################################################
 #########
 #########            TRICLINIC: 1 to 2
 #########
 #########################################################################################################################
 # Nothing to do for sym_number 1
 # Nothing to do for sym_number 2
 #########################################################################################################################
 #########
 #########            MONOCLINIC: 3 to 15 
 #########
 #########################################################################################################################
 # Nothing to do for sym_number 3
 if (sym_number == 4)
 {
  if (setting == 1) hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 2) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 3) hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 5)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 10) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 11) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 6
 if (sym_number == 7)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 2,3] <- NA
  if (setting == 2) hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 4) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 5) hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 7) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 8) hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 9) hkl[hkl[,1] == 0 & hkl[,3]%%2,3] <- NA
 }
 if (sym_number == 8)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA 
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 10) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 9)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 7)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 8)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 9)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 10)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 11)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 12)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 13)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 14)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 15)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 16)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 17)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 18)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 # Nothing to do for sym_number 10, 11
 if (sym_number == 12)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 13)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 2) hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 4) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 5) hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 7) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 8) hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 9) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 14)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 7)
  {
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 8)
  {
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 9)
  {
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 15)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 7)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 8)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 9)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 10)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 11)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 12)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 13)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 14)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 15)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 16)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 17)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 18)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 #########################################################################################################################
 #########
 #########            ORTHOROMBIC: 16 to 74 
 #########
 #########################################################################################################################
 # Nothing to do for sym_number 16
 if (sym_number == 17)
 {
  if (setting == 1) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 2) hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 3) hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 18
 if (sym_number == 19)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 if (sym_number == 20)
 { 
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 21)
 { 
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- 3
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- 3
  if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- 3
 }
 if (sym_number == 22)
 {
  hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 if (sym_number == 23) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 24)
 {
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 25
 if (sym_number == 26)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 2) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 3) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 27)
 {
  if (setting == 1) 
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2) 
  {
   hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3) 
  {
   hkl[hkl[,3] == 0 & hkl[,1] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 28)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 2 | setting == 3) hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 4) hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 5) hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 29)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,3] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 30)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 31)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 32)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 33)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 34)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 35)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 36)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 37)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 38)
 {
  if (setting == 1) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 2) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 4) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 6) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 39)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 40)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 41)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 42)
 {
  hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 if (sym_number == 43)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%4 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%4 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%4 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%4 != 0,4] <- NA
  }
 }
 if (sym_number == 44) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 45)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 46)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 # Nothing to do for sym_number 47
 if (sym_number == 48)
 {
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 49)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 50)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 51)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 52)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 53)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 54)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 55)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 56)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 57)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 58)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 59)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 60)
 {
  if (setting == 1)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {  
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {  
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 61)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 62)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 63)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 2) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 3) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 4) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
  if (setting == 5) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  if (setting == 1 | setting == 2) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 3 | setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 5 | setting == 6) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 64)
 {
  if (setting == 1) 
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 2) 
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 3) 
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 4) 
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 5) 
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (setting == 6) 
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (setting == 1 | setting == 2) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 3 | setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting == 5 | setting == 6) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 65)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,l])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,l])%%2 != 0,4] <- NA
 }
 if (sym_number == 66)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,l])%%2 != 0,4] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,l])%%2 != 0,4] <- NA
 }
 if (sym_number == 67)
 {
  if (setting == 1 | setting == 2 | setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 4 | setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
 }
 if (sym_number == 68)
 {
  if (setting == 1 | setting == 2 | setting == 3 | setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 5 | setting == 6 | setting == 7)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 8)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 9 | setting == 10 | setting == 11)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 12)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting >= 1 & setting <= 4) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (setting >= 5 & setting <= 8) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (setting >= 9 & setting <= 12) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 69)
 {
  hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 if (sym_number == 70)
 {
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (3*(hkl[,1]+hkl[,2]))%%4 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (3*(hkl[,2]+hkl[,3]))%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (3*(hkl[,1]+hkl[,3]))%%4 != 0,4] <- NA
  hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 if (sym_number == 71)
 {
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 72)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 73)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 74)
 {
  if (setting == 1 | setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 #########################################################################################################################
 #########
 #########            TETRAGONAL: 75 to 142 
 #########
 #########################################################################################################################
 # Nothing to do for sym_number 75
 if (sym_number == 76) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
 if (sym_number == 77) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 if (sym_number == 78) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
 if (sym_number == 79) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 80) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 81
 if (sym_number == 82) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 # Nothing to do for sym_number 83
 if (sym_number == 84) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 if (sym_number == 85) hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 if (sym_number == 86) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 if (sym_number == 87) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 88)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 89
 if (sym_number == 90)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 if (sym_number == 91) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 if (sym_number == 92)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 if (sym_number == 93) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 if (sym_number == 94)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 if (sym_number == 95) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
 if (sym_number == 96)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 if (sym_number == 97) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 98)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 99
 if (sym_number == 100)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 101)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 102) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 103) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 104) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 105) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 106)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 107) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 108)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 109)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 110)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 111
 if (sym_number == 112) hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 if (sym_number == 113)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 if (sym_number == 114)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 115
 if (sym_number == 116)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 117)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 118)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 119) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 120)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 121) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 122)
 {
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 123
 if (sym_number == 124)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 125)
 {
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & +hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & +hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 126)
 {
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 127)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
 }
 if (sym_number == 128)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 129)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
 }
 if (sym_number == 130)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 131)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 132)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 133)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 134)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 135)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 136)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 137)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 138)
 {
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 139) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 if (sym_number == 140)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 141)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 if (sym_number == 142)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
 }
 #########################################################################################################################
 #########
 #########            TRIGONAL: 143 to 167 
 #########
 #########################################################################################################################
 # Nothing to do for sym_number 143
 if (sym_number == 144) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
 if (sym_number == 145) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
 if (sym_number == 146) 
 {
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
  # Nothing to do for setting 2
 }
 # Nothing to do for sym_number 147
 if (sym_number == 148) 
 {
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
  # Nothing to do for setting 2
 }
 # Nothing to do for sym_number 149
 # Nothing to do for sym_number 150
 if (sym_number == 151) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
 if (sym_number == 152) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
 if (sym_number == 153) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
 if (sym_number == 154) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
 if (sym_number == 155) 
 {
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
  # Nothing to do for setting 2
 }
 # Nothing to do for sym_number 156
 # Nothing to do for sym_number 157
 if (sym_number == 158)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 159)
 {
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 160) 
 {
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
  # Nothing to do for setting 2
 }
 if (sym_number == 161) 
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
  # Nothing else for setting 2
 }
 # Nothing to do for sym_number 162
 if (sym_number == 163)
 {
  hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 # Nothing to do for sym_number 164
 if (sym_number == 165)
 {
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 if (sym_number == 166) 
 {
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
 }
 if (sym_number == 167) 
 {
  if (setting == 1)
  {
   hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
  }
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
 }
 #########################################################################################################################
 #########
 #########            HEXAGONAL: 168 to 194 
 #########
 #########################################################################################################################
 #########################################################################################################################
 #########
 #########            CUBIC: 195 to 230 
 #########
 #########################################################################################################################

 # Delete systematic absences
 idx <- which(complete.cases(hkl))

 return(idx)
}
