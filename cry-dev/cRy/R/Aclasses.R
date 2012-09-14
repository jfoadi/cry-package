###################################################################
###################################################################
### This module is part of the cRy package for crystallography. ###
### Authors: J. Foadi & D. G. Waterman                          ###
### MPL - Imperial College London / Diamond Light Source Ltd    ###
### CCP4 - Research Complex at Harwell                          ###
###                                                             ###
### Classes not needing anything else to function. Thus, they   ###
### are collated at the top.                                    ###
###################################################################
###################################################################


## Angle
#
setClass(
         Class="Angle",
         representation=representation(ang="numeric",rad_flag="logical"),
         validity=function(object)
                  {
                   # Extract slots' content, even if they're empty
                   value <- object@ang
                   flag <- object@rad_flag

                   # Check angle is a number between 0 and 180 (or 0 and pi)
                   if (length(value) != 0 & length(flag) == 0)
                   {
                    if (object@ang < 0) return("Numerical value for this angle is not allowed because negative")
                   }
                   if (length(value) != 0 & length(flag) != 0)
                   {
                    if (!object@rad_flag & (value <= 0 | value >= 180)) return("Numerical value for this angle is not in the 0 - 180 range")
                    if (object@rad_flag & (value <= 0 | value >= pi)) return("Numerical value for this angle is not in the 0 - pi range")
                   }

                   # Green light: object appears to be correct
                   return(TRUE)
                  }
        )

## UnitCell (angles always stored in degrees)
#
setClass(
         Class="UnitCell",
         representation=representation(a="numeric",b="numeric",c="numeric",alpha="Angle",beta="Angle",gamma="Angle"),
         validity=function(object)
                  {
                   # Angles are always stored in degrees
                   if (length(object@alpha@ang) != 0)
                   {
                    if (object@alpha@rad_flag) stop("Cell angles always stored in degrees for object of class UnitCell.")
                   }
                   if (length(object@beta@ang) != 0)
                   {
                    if (object@beta@rad_flag) stop("Cell angles always stored in degrees for object of class UnitCell.")
                   }
                   if (length(object@gamma@ang) != 0)
                   {
                    if (object@gamma@rad_flag) stop("Cell angles always stored in degrees for object of class UnitCell.")
                   }

                   # Cell sides can be any number different from zero or infinity
                   if (length(object@a) != 0)
                   {
                    if (object@a <= 0 | !is.finite(object@a)) stop("Cell side a can be neither zero nor infinity.")
                   }
                   if (length(object@b) != 0)
                   {
                    if (object@b <= 0 | !is.finite(object@b)) stop("Cell side b can be neither zero nor infinity.")
                   }
                   if (length(object@c) != 0)
                   {
                    if (object@c <= 0 | !is.finite(object@c)) stop("Cell side c can be neither zero nor infinity.")
                   }

                   # Angles are objects of class Angles, so they all have values between 0 and 180 degrees. But combinations
                   # of alpha beta and gamma need to obey certain rules (see J. Foadi & G. Evans (2011))
                   if (length(object@alpha@ang) != 0 & length(object@beta@ang) != 0 & length(object@gamma@ang) != 0)
                   {
                    aa <- object@alpha@ang
                    bb <- object@beta@ang
                    cc <- object@gamma@ang
                    ss <- aa+bb+cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                    ss <- aa+bb-cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                    ss <- aa-bb+cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                    ss <- -aa+bb+cc
                    if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist")
                   }

                   # Green light: object appears to be correct
                   return(TRUE)
                  }
        )

## ReciprocalUnitCell (angles always stored in degrees)
#
setClass(
         Class="ReciprocalUnitCell",
         representation=representation(ar="numeric",br="numeric",cr="numeric",alphar="Angle",betar="Angle",gammar="Angle"),
         validity=function(object)
                  {
                   # Angles are always stored in degrees
                   if (length(object@alphar@ang) != 0)
                   {
                    if (object@alphar@rad_flag) stop("Cell angles always stored in degrees for object of class ReciprocalUnitCell.")
                   }
                   if (length(object@betar@ang) != 0)
                   {
                    if (object@betar@rad_flag) stop("Cell angles always stored in degrees for object of class ReciprocalUnitCell.")
                   }
                   if (length(object@gammar@ang) != 0)
                   {
                    if (object@gammar@rad_flag) stop("Cell angles always stored in degrees for object of class ReciprocalUnitCell.")
                   }

                   # Cell sides can be any number different from zero or infinity
                   if (length(object@ar) != 0)
                   {
                    if (object@ar <= 0 | !is.finite(object@ar)) stop("Cell side ar can be neither zero nor infinity.")
                   }
                   if (length(object@br) != 0)
                   {
                    if (object@br <= 0 | !is.finite(object@br)) stop("Cell side br can be neither zero nor infinity.")
                   }
                   if (length(object@cr) != 0)
                   {
                    if (object@cr <= 0 | !is.finite(object@cr)) stop("Cell side cr can be neither zero nor infinity.")
                   }

                   # Angles are objects of class Angles, so they all have values between 0 and 180 degrees. But combinations
                   # of alphar betar and gammar need to obey certain rules
                   if (length(object@alphar@ang) != 0 & length(object@betar@ang) != 0 & length(object@gammar@ang) != 0)
                   {
                    aar <- object@alphar@ang
                    bbr <- object@betar@ang
                    ccr <- object@gammar@ang
                    ssr <- aar+bbr+ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                    ssr <- aar+bbr-ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                    ssr <- aar-bbr+ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocal unit cell with these angles cannot exist")
                    ssr <- -aar+bbr+ccr
                    if (ssr <= 0 | ssr >= 360) stop("A reciprocalunit cell with these angles cannot exist")
                   }

                   # Green light: object appears to be correct
                   return(TRUE)
                  }
        )

## BravaisType
#
setClass(
         Class="BravaisType",
         representation=representation(bl="character"),
         validity=function(object)
         {
          # Extract lattice and centering (also check on symbol validity)
          if (length(object@bl) != 0)
          {
           latt <- substr(object@bl,1,1)
           centring <- substr(object@bl,2,2)
           if (latt != "a" & latt !="m" & latt != "o" & latt != "t" & latt != "h" & latt != "c") stop("Invalid lattice type.")
           if (centring == "S") stop("Please, specify if it is an A, B or C centring")
           if (centring != "P" & centring != "A" & centring != "B" & centring != "C" & centring != "F" & centring != "I" & centring != "R") stop("Invalid centring type")

           # Now check validity for specific families
           # Triclinic family
           if (latt == "a" & centring != "P") stop("Invalid Bravais lattice type")

           # Monoclinic family
           #if ((latt == "m" & centring == "I") | (latt == "m" & centring == "F") | (latt == "m" & centring == "R")) stop("Invalid Bravais lattice type")
           if ((latt == "m" & centring == "F") | (latt == "m" & centring == "R")) stop("Invalid Bravais lattice type")

           # Orthorombic family
           if ((latt == "o" & centring != "P") & (latt == "o" & centring != "A") & (latt == "o" & centring != "B") & (latt == "o" & centring != "C") &
               (latt == "o" & centring != "F") & (latt == "o" & centring != "I")) stop("Invalid Bravais lattice type")

           # Tetragonal family
           if ((latt == "t" & centring != "P") & (latt == "t" & centring != "I")) stop("Invalid Bravais lattice type")

           # Hexagonal and Trigonal families
           if (latt == "h"  & (centring != "P" & centring != "R")) stop("Invalid Bravais lattice type")

           # Cubic family
           if (latt == "c" & (centring == "A" | centring == "B" | centring == "C")) stop("Invalid Bravais lattice type")
          }

          # Green light: object appears to be correct
          return(TRUE)
         }
        )

## Lattice
#
setClass(
         Class="Lattice",
         representation=representation(cell="UnitCell",bl="BravaisType"),
         validity=function(object)
         {
          if (length(object@cell@a) != 0 & length(object@cell@b) != 0 & length(object@cell@c) != 0 &
              length(object@cell@alpha@ang) != 0 & length(object@cell@beta@ang) != 0 & length(object@cell@gamma@ang) != 0 &
              length(object@bl@bl) != 0)
          {
           # Extract cell parameters
           a <- object@cell@a
           b <- object@cell@b
           c <- object@cell@c
           aa <- object@cell@alpha@ang
           bb <- object@cell@beta@ang
           cc <- object@cell@gamma@ang

           # Extract lattice and centring
           latt <- substr(object@bl@bl,1,1)
           centring <- substr(object@bl@bl,2,2)

           # Compatibility between unit cell and lattice type
           if (latt == "m")                                    # Monoclinic family
           {
            erang <- 0.000001      # Finite accuracy of binary representation of numbers means we have to
            diff_a <- abs(aa-90)   # test number == 90 in this way
            diff_b <- abs(bb-90)
            diff_c <- abs(cc-90)
            if (diff_a > erang & (diff_b > erang | diff_c > erang)) stop ("Unit cell angles non compatible with lattice type")
            if (diff_b > erang & (diff_a > erang | diff_c > erang)) stop ("Unit cell angles non compatible with lattice type")
            if (diff_c > erang & (diff_a > erang | diff_b > erang)) stop ("Unit cell angles non compatible with lattice type")
            if (diff_a > erang & centring == "A") stop ("Unit cell parameters non compatible with lattice type")
            if (diff_b > erang & centring == "B") stop ("Unit cell parameters non compatible with lattice type")
            if (diff_c > erang & centring == "C") stop ("Unit cell parameters non compatible with lattice type")
           }
           if (latt == "o")                                    # Orthorombic family
           {
            if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
            {
             stop ("Unit cell parameters non compatible with lattice type")
            }
           }
           if (latt == "t")                                    # Tetragonal family
           {
            if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
            {
             stop ("Unit cell parameters non compatible with lattice type")
            }
            if (abs(a-b) > 0.000001 & abs(a-c) > 0.000001 & abs(b-c) > 0.000001)
            {
             stop ("Unit cell parameters non compatible with lattice type")
            }
           }
           if (latt == "c")                                    # Cubic family
           {
            if (abs(a-b) > 0.000001 | abs(a-c) > 0.000001 | abs(b-c) > 0.000001)
            {
             stop ("Unit cell parameters non compatible with lattice type")
            }
            if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
            {
             stop ("Unit cell parameters non compatible with lattice type")
            }
           }
           if (latt == "h" & centring == "P")                    # Hexagonal family
           {
            if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
            if (abs(cc-120) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
            if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
           }
           if (latt == "h" & centring == "R")                    # Rombohedral family
           {
            if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
            if (abs(a-b) < 0.000001)
            {
             if (abs(a-c) < 0.000001)
             {
              if (abs(aa-bb) > 0.000001 | abs(aa-cc) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
              if (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)
              {
               messaggio <- paste("Settings for this lattice are for a cubic, rather than a rombohedral system.",
                                  "Change BravaisType if you wish to carry on with this object")
               stop(messaggio)
              }
             }
             if (abs(a-c) > 0.000001)
             {
              if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 & abs(cc-120) > 0.000001)
                  stop ("Unit cell parameters non compatible with lattice type")
             }
            }
           }
          }

          # Green light: object appears to be correct
          return(TRUE)
         }
        )

## Symmetry class
#
# Symmetry
setClass(
         Class="Symmetry",
         representation=representation(sym_xHM="character"),
         validity=function(object)
         {
          # If object has been created empty, or with some empty slots, fine.
          if (length(object@sym_xHM) == 0)
          {
           # Green light: object appears to be correct
           return(TRUE)
          }

          # Check extended Herman-Maguin symbol is correct
          if (length(object@sym_xHM) != 0)
          {
           # To be used throughout this block
           tmp <- .translate_SG(object@sym_xHM)

           lista <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
           if (!lista$ans) stop(lista$msg)
          }

          # Green light: object appears to be correct
          return(TRUE)
         }
        )
