###################################################################
###################################################################
### This module is part of the cRy package for crystallography. ###
### Authors: J. Foadi & D. G. Waterman                          ###
###                                                             ###
### Secondary generics for methods included in methods code.    ###
###################################################################
###################################################################



######### For Angle class
#
# *** [degToRad] *** Convert degrees to radians. Angle object as input.
setGeneric(
           name="degToRad",
           def=function(object){standardGeneric("degToRad")} 
          )
#
# *** [radToDeg] *** Convert radians to degrees. Angle object as input.
setGeneric(
           name="radToDeg",
           def=function(object){standardGeneric("radToDeg")}
          )


######### For UnitCell class
#
# *** [computeReciprocalUnitCell] *** Compute reciprocal cell starting from direct cell
setGeneric(
           name="computeReciprocalUnitCell",
           def=function(object){standardGeneric("computeReciprocalUnitCell")}
          )


######### For ReciprocalUnitCell class
#
# *** [computeUnitCell] *** Compute direct cell starting from reciprocal cell
setGeneric(
           name="computeUnitCell",
           def=function(object){standardGeneric("computeUnitCell")}
          )


######### For both UnitCell and  ReciprocalUnitCell class
#
# *** [computeCellVolume] *** Compute cell volume (valid for both direct and reciprocal cells)
setGeneric(
           name="computeCellVolume",
           def=function(object){standardGeneric("computeCellVolume")}
          )


######### For UnitCell and  Symmetry class
#
# *** [checkSymmetryWithCell] *** Check symmetry is compatible with unit cell parameters
setGeneric(
           name="checkSymmetryWithCell",
           def=function(obUnitCell,obSymmetry){standardGeneric("checkSymmetryWithCell")}
          )


######### For BravaisType, Lattice and  Symmetry class
#
# *** [getCentringOps] *** Returns centering operators
setGeneric(
           name="getCentringOps",
           def=function(object){standardGeneric("getCentringOps")}
          )
