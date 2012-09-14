# R code for classes used in cRy
# This code uses S4 classes formalism.
#
# J. Foadi (Imperial College) and D. Waterman (Diamond Light Source)
#

#
# Angle
setClass(
         Class="Angle",
         representation=representation(ang="numeric",rad_flag="logical")
        )
#
# UnitCell
setClass(
         Class="UnitCell",
         representation=representation(a="numeric",b="numeric",c="numeric",alpha="Angle",beta="Angle",gamma="Angle")
        )
#
# ReciprocalUnitCell
setClass(
         Class="ReciprocalUnitCell",
         representation=representation(ar="numeric",br="numeric",cr="numeric",alphar="Angle",betar="Angle",gammar="Angle")
        )
#
# BravaisType
setClass(
         Class="BravaisType",
         representation=representation(bl="character")
        )
#
# Lattice
setClass(
         Class="Lattice",
         representation=representation(cell="UnitCell",bl="BravaisType")
        )
