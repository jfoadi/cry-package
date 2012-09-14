###################################################################
###################################################################
### This module is part of the cRy package for crystallography. ###
### Authors: J. Foadi & D. G. Waterman                          ###
### MPL - Imperial College London / Diamond Light Source Ltd    ###
### CCP4 - Research Complex at Harwell                          ###
###                                                             ###
### All generics for methods included in methods code.          ###
###################################################################
###################################################################



## Valid for all classes:

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

## User-friendly constructors

## For Angle class
#
# *** [angle] *** 2 possible variables as input
setGeneric(
           name="angle",
           def=function(ang,rad_flag,...){standardGeneric("angle")}
          )

## For UnitCell class
#
# *** [unitcell] *** 6 possible variables as input (6 cell parameters or the first can be a file name)
setGeneric(
           name="unitcell",
           def=function(a,b,c,alpha,beta,gamma,...){standardGeneric("unitcell")}
          )

## For ReciprocalUnitCell class
#
# *** [reciprocalunitcell] *** 6 possible variables as input (6 reciprocal-cell parameters or the first can be a file name)
setGeneric(
           name="reciprocalunitcell",
           def=function(ar,br,cr,alphar,betar,gammar,...){standardGeneric("reciprocalunitcell")}
          )

## For BravaisType class
#
# *** [bravaistype] *** 1 possible variable as input (either a bl character or a file name)
setGeneric(
           name="bravaistype",
           def=function(x,...){standardGeneric("bravaistype")}
          )

## For Lattice class
#
# *** [lattice] *** 6 possible variables as input
setGeneric(
           name="lattice",
           def=function(x1,x2,x3,x4,x5,x6,...){standardGeneric("lattice")}
          )