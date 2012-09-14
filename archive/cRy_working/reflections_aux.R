sysabs  <- function(hkl,sym_number,setting)
{
 # Input: hkl is a m X 3 matrix whose lines are Miller indices
 #        sym_number is the symmetry number in the International Tables for Crystallography
 #        setting is the setting as in the "syminfo.lib" file
 # Output is a vector containing indices of Miller indices obeying systematic absences

 # Long list for all space groups and settings

 # Nothing to do for sym_number 1,2,3
 if (sym_number == 4)
 {
  if (setting == 1) hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  if (setting == 2) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 3) hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
 }
 if (sym_number == 5)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 10) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 11) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
 }
 # Nothing to do for sym_number 6
 if (sym_number == 7)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 2,3] <- NA
  if (setting == 2) hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  if (setting == 4) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  if (setting == 5) hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  if (setting == 7) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  if (setting == 8) hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 9) hkl[hkl[,1] == 0 & hkl[,3]%%2,3] <- NA
 }
 if (sym_number == 8)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA 
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 10) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
 }
 if (sym_number == 9)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 7)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 8)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 9)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 10)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 11)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 12)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 13)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 14)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 15)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 16)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 17)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 18)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 # Nothing to do for sym_number 10, 11
 if (sym_number == 12)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
 }
 if (sym_number == 13)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 2) hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  if (setting == 4) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  if (setting == 5) hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  if (setting == 7) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  if (setting == 8) hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 9) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
 }
 if (sym_number == 14)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 7)
  {
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 8)
  {
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 9)
  {
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 15)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 7)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 8)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 9)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 10)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 11)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 12)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 13)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 14)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 15)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 16)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 17)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 18)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 # Nothing to do for sym_number 16
 if (sym_number == 17)
 {
  if (setting == 1) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 2) hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  if (setting == 3) hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
 }
 # Nothing to do for sym_number 18
 if (sym_number == 19)
 {
  hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
 }
 if (sym_number == 20)
 { 
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 21)
 { 
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- 3
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- 3
  if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- 3
 }
 if (sym_number == 22)
 {
  hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
 }
 if (sym_number == 23) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
 if (sym_number == 24)
 {
  hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
 }
 # Nothing to do for sym_number 25
 if (sym_number == 26)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 2) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 3) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
 }
 if (sym_number == 27)
 {
  if (setting == 1) 
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2) 
  {
   hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3) 
  {
   hkl[hkl[,3] == 0 & hkl[,1] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 28)
 {
  if (setting == 1) hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,3] <- NA
  if (setting == 2 | setting == 3) hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,3] <- NA
  if (setting == 4) hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 5) hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  if (setting == 6) hkl[hkl[,3] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,3] <- NA
 }
 if (sym_number == 29)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,3] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 30)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 31)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 32)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 33)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 34)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 35)
 {
  if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
 }
 if (sym_number == 36)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 37)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 38)
 {
  if (setting == 1) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 2) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  if (setting == 4) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 5) hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  if (setting == 6) hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
 }
 if (sym_number == 39)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 40)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 41)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 42)
 {
  hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
 }
 if (sym_number == 43)
 {
  if (setting == 1)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%4 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%4 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%4 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%4 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[(hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[(hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%4 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%4 != 0,3] <- NA
  }
 }
 if (sym_number == 44) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
 if (sym_number == 45)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 46)
 {
  if (setting == 1)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 # Nothing to do for sym_number 47
 if (sym_number == 48)
 {
  hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
 }
 if (sym_number == 49)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 50)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 51)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 52)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 53)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 54)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 55)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 56)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 57)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 58)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 59)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 60)
 {
  if (setting == 1)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {  
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {  
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {  
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 61)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
  }
 }
 if (sym_number == 62)
 {
  if (setting == 1)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 2)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 3)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,3] <- NA
  }
  if (setting == 4)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 5)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,3] <- NA
  }
  if (setting == 6)
  {
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,3] <- NA
   hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,3] <- NA
   hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,3] <- NA
  }
 }

 # Delete systematic absences
 idx <- which(complete.cases(hkl))

 return(idx)
}
