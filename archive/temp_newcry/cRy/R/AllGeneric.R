# R code for package cRy


#
# *** [isCorrect] *** Check any object is valid.
setGeneric(
           name="isCorrect",
           def=function(object,name,message=FALSE){standardGeneric("isCorrect")}
          )
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

# Specific generic functions for class Angle
#
# Degrees to radians conversion
setGeneric(
           name="degToRad",
           def=function(object){standardGeneric("degToRad")} 
          )
#
# Radians to degrees conversion
setGeneric(
           name="radToDeg",
           def=function(object){standardGeneric("radToDeg")} 
          )

# Specific generic functions for class UnitCell and ReciprocalUnitCell
#
# Compute reciprocal cell starting from direct cell
setGeneric(
           name="computeReciprocalUnitCell",
           def=function(object){standardGeneric("computeReciprocalUnitCell")}
          )
#
# Compute direct cell starting from reciprocal cell
setGeneric(
           name="computeUnitCell",
           def=function(object){standardGeneric("computeUnitCell")}
          )
#
# Compute cell volume (valid for both direct and reciprocal cells)
setGeneric(
           name="computeCellVolume",
           def=function(object){standardGeneric("computeCellVolume")}
          )


# Specific generic functions for class Lattice and Symmetry
#
# Returns centering operators
setGeneric(
           name="getCentringOps",
           def=function(object){standardGeneric("getCentringOps")} 
          )
