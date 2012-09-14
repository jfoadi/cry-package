# R code for package cRy


#
# *** [getFields] *** Extract all object's slots
setGeneric(
           name="getFields",
           def=function(object){standardGeneric("getFields")}
          )
#
# *** [getParameters] *** Get the essential parameters from the specific object
setGeneric(
           name="getParameters",
           def=function(object){standardGeneric("getParameters")}
          )
#
# *** [getData] *** Get data from the specific object
setGeneric(
           name="getData",
           def=function(object){standardGeneric("getData")}
          )
#
# *** [getCentringOps] *** Return centring operators. Applies to Lattice and Symmetry classes
setGeneric(
           name="getCentringOps",
           def=function(object){standardGeneric("getCentringOps")}
          )
#
# *** [checkSymmetryWithCell] *** Check symmetry is compatible with unit cell parameters
setGeneric(
           name="checkSymmetryWithCell",
           def=function(obSymmetry,obUnitCell){standardGeneric("checkSymmetryWithCell")}
          )
