# R code to implement the crystallographic ideas related to merged.
# and unmerged reflections.
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

## MergedReflections class
#
# MergedReflections
setClass(
         Class="MergedReflections",
         representation=representation(cell="UnitCell",symmetry="Symmetry",records="data.frame",dtypes="character"),
         validity=function(object)
         {
          # Remember, in every check you'll have to take into account the possibility that the slot is empty (its length is zero)

          # records slot needs to have more than 3 columns. It always starts with "H", "K", "L"
          if (length(object@records) != 0)
          {
           if (length(object@records) <= 3) stop("Slot records in object of class MergedReflections needs to have more than 3 columns")
           hklnames <- colnames(object@records)[1:3]
           if (hklnames[1] != "H" & hklnames[1] != "h")
           {
            msg <- "First column of slot records in object of class MergedReflections can have label name H or h"
            stop(msg)
           }
           if (hklnames[2] != "K" & hklnames[2] != "k")
           {
            msg <- "Second column of slot records in object of class MergedReflections can have label name K or k"
            stop(msg)
           }
           if (hklnames[3] != "L" & hklnames[3] != "l")
           {
            msg <- "Third column of slot records in object of class MergedReflections can have label name L or l"
            stop(msg)
           }
          }

          # First 3 characters of slot dtypes have to be "H", "H", "H"
          if (length(object@dtypes) != 0)
          {
           if (length(object@dtypes) <= 3) 
           {
            msg <- "dtypes slot for object of class MergedReflections needs to have more than 3 columns"
            stop(msg)
           }
           if (object@dtypes[1] != "H" | object@dtypes[2] != "H" | object@dtypes[3] != "H")
           {
            msg <- "The first dtypes for object of class MergedReflections need always to be H H H"
            stop(msg)
           }
          }

          # Length of dtypes needs to be equal to length of records
          if (length(object@records) != 0 & length(object@dtypes) != 0)
          {
           if (length(object@records) != length(object@dtypes))
           {
            msg <- paste("Number of columns of slot records in object of class MergedReflections has to be equal to number",
                         "of elements of slot dtypes")
            stop(msg)
           }
          }

          # Green light: object appears to be correct
          return(TRUE)
         }
        )
