###################################################################
###################################################################
### This module is part of the cRy package for crystallography. ###
### Authors: J. Foadi & D. G. Waterman                          ###
### MPL - Imperial College London / Diamond Light Source Ltd    ###
### CCP4 - Research Complex at Harwell                          ###
###                                                             ###
### Methods collated at the top of code.                        ###
###################################################################
###################################################################



######### PRINT ######### PRINT ######### PRINT ######### PRINT ######### PRINT ######### PRINT ######### PRINT ######### PRINT #########

## For Angle class
#
setMethod(
          f="print",
          signature="Angle",
          definition=function(x,...)
                     {
                      if (length(x@rad_flag) != 0 & length(x@rad_flag) != 0)
                      {
                       if (!x@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.2f degrees. ***\n",x@ang)
                       if (x@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.4f radians. ***\n",x@ang)
                       cat(stringa)
                      }
                      if (length(x@rad_flag) == 0 & length(x@ang) != 0)
                      {
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is %6.2f, but it is unclear",x@ang)
                       stringa2 <- sprintf("whether this value is in degrees or radians. ***\n")
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(x@rad_flag) != 0 & length(x@ang) == 0)
                      {
                       if (x@rad_flag) fchar <- "radians"
                       if (!x@rad_flag) fchar <- "degrees"
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is not known, but")
                       stringa2 <- sprintf("it is measured in %s. ***\n",fchar)
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(x@rad_flag) == 0 & length(x@ang) == 0)
                      {
                       stringa <- sprintf("*** This is an object of class Angle. Its slots are both empty. ***\n")
                       cat(stringa)
                      }
                     }
         )

## For UnitCell class
#
setMethod(
          f="print",
          signature="UnitCell",
          definition=function(x,...)
                     {
                      cat("*** This is a UnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(x@a) == 0) cat(" Slot for parameter a is empty.\n")
                      if (length(x@a) != 0) cat(sprintf(" Parameter a = %8.3f\n",x@a))
                      if (length(x@b) == 0) cat(" Slot for parameter b is empty.\n")
                      if (length(x@b) != 0) cat(sprintf(" Parameter b = %8.3f\n",x@b))
                      if (length(x@c) == 0) cat(" Slot for parameter c is empty.\n")
                      if (length(x@c) != 0) cat(sprintf(" Parameter c = %8.3f\n",x@c))
                      if (length(x@alpha@ang) == 0) cat(" Slot for parameter alpha is empty.\n")
                      if (length(x@alpha@ang) != 0) cat(sprintf(" Parameter alpha = %6.2f degrees\n",x@alpha@ang))
                      if (length(x@beta@ang) == 0) cat(" Slot for parameter beta is empty.\n")
                      if (length(x@beta@ang) != 0) cat(sprintf(" Parameter beta = %6.2f degrees\n",x@beta@ang))
                      if (length(x@gamma@ang) == 0) cat(" Slot for parameter gamma is empty.\n")
                      if (length(x@gamma@ang) != 0) cat(sprintf(" Parameter gamma = %6.2f degrees\n",x@gamma@ang))
                     }
         )

## For ReciprocalUnitCell class
#
setMethod(
          f="print",
          signature="ReciprocalUnitCell",
          definition=function(x,...)
                     {
                      cat("*** This is a ReciprocalUnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(x@ar) == 0) cat(" Slot for parameter ar is empty.\n")
                      if (length(x@ar) != 0) cat(sprintf(" Parameter ar = %8.3f\n",x@ar))
                      if (length(x@br) == 0) cat(" Slot for parameter br is empty.\n")
                      if (length(x@br) != 0) cat(sprintf(" Parameter br = %8.3f\n",x@br))
                      if (length(x@cr) == 0) cat(" Slot for parameter cr is empty.\n")
                      if (length(x@cr) != 0) cat(sprintf(" Parameter cr = %8.3f\n",x@cr))
                      if (length(x@alphar@ang) == 0) cat(" Slot for parameter alphar is empty.\n")
                      if (length(x@alphar@ang) != 0) cat(sprintf(" Parameter alphar = %6.2f degrees\n",x@alphar@ang))
                      if (length(x@betar@ang) == 0) cat(" Slot for parameter betar is empty.\n")
                      if (length(x@betar@ang) != 0) cat(sprintf(" Parameter betar = %6.2f degrees\n",x@betar@ang))
                      if (length(x@gammar@ang) == 0) cat(" Slot for parameter gammar is empty.\n")
                      if (length(x@gammar@ang) != 0) cat(sprintf(" Parameter gammar = %6.2f degrees\n",x@gammar@ang))
                     }
         )

## For BravaisType class
#
setMethod(
          f="print",
          signature="BravaisType",
          definition=function(x,...)
                     {
                      # Initial message
                      cat("*** This is a BravaisType object. ***\n")

                      # If object is empty print message
                      if (length(x@bl) == 0) cat(" The slot for parameter bl is empty.\n")

                      # Extract bl character
                      if (length(x@bl) != 0)
                      {
                       bl <- x@bl
                       cat(paste(" The slot for parameter bl contains character ",bl,".\n",sep=""))
                      }
                     }
         )

## For Lattice class
#
setMethod(
          f="print",
          signature="Lattice",
          definition=function(x,...)
                     {
                      # Initial message
                      cat("*** This is a Lattice object. ***\n")

                      # Print empty or cell values
                      cat(" The unit cell has the following parameters:\n")
                      if (length(x@cell@a) == 0) cat(" a is empty.\n")
                      if (length(x@cell@a) != 0) cat(sprintf(" a = %7.2f\n",x@cell@a))
                      if (length(x@cell@b) == 0) cat(" b is empty.\n")
                      if (length(x@cell@b) != 0) cat(sprintf(" b = %7.2f\n",x@cell@b))
                      if (length(x@cell@c) == 0) cat(" c is empty.\n")
                      if (length(x@cell@c) != 0) cat(sprintf(" c = %7.2f\n",x@cell@c))
                      if (length(x@cell@alpha@ang) == 0) cat(" alpha is empty.\n")
                      if (length(x@cell@alpha@ang) != 0) cat(sprintf(" alpha = %6.2f degrees.\n",x@cell@alpha@ang))
                      if (length(x@cell@beta@ang) == 0) cat(" beta is empty.\n")
                      if (length(x@cell@beta@ang) != 0) cat(sprintf(" beta = %6.2f degrees.\n",x@cell@beta@ang))
                      if (length(x@cell@gamma@ang) == 0) cat(" gamma is empty.\n")
                      if (length(x@cell@gamma@ang) != 0) cat(sprintf(" gamma = %6.2f degrees.\n",x@cell@gamma@ang))
                      if (length(x@bl@bl) == 0) cat(" Information on the lattice is not vailable.\n")

                      # Extract lattice information (if lattice has been defined)
                      if (length(x@bl@bl) != 0)
                      {
                       latt_info <- .extractLatticeStuff(x)

                       # Print
                       stringa <- sprintf("\n The Bravais lattice is %s",latt_info[1])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal family is %s",latt_info[2])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal system is %s",latt_info[3])
                       cat(stringa)
                       stringa <- sprintf("\n The lattice system is %s\n",latt_info[4])
                       cat(stringa)
                      }
                     }
         )

## For Symmetry class
#
setMethod(
          f="print",
          signature="Symmetry",
          definition=function(x,...)
                     {
                      # Initial message
                      cat("*** This is a Symmetry object ***\n")

                      # Empty object
                      if (length(x@sym_xHM) == 0)
                      {
                       stringa <- sprintf(" The slot for its extended Herman-Maguin symbol is empty.\n")
                       cat(stringa)
                      }

                      if (length(x@sym_xHM) != 0)
                      {
                       # Find space group number
                       lista <- .translate_SG(value=x@sym_xHM,SG_in="xHM",SG_out="number")
                       sgn <- lista$msg

                       # Find crystal system
                       cr_sys <- .crystal_system(sgn)

                       # Find symmetry and centring operations
                       lista <- .syminfo_to_op_xyz_list(x@sym_xHM)
                       xyz <- lista[[1]]
                       cenop <- lista[[2]]

                       # Print
                       stringa <- sprintf("    The space group extended Hermann-Maguin symbol is %s\n",x@sym_xHM)
                       cat(stringa)
                       stringa <- sprintf("    The space group number from the International Tables is %d\n",sgn)
                       cat(stringa)
                       stringa <- sprintf("    This space group belongs to the %s crystal system\n",cr_sys)
                       cat(stringa)
                       stringa <- sprintf("    This space group includes the following symmetry operations:\n")
                       cat(stringa)
                       for (i in 1:length(xyz))
                       {
                        stringa <- sprintf("    %2d) %s\n",i,xyz[i])
                        cat(stringa)
                       }
                       if (length(cenop) == 1)
                       {
                        stringa <- sprintf("    There are no centring operators associated with this group\n")
                        cat(stringa)
                       }
                       if (length(cenop) != 1)
                       {
                        stringa <- sprintf("    This space group is associated with a centred cell. The centring operators are:\n")
                        cat(stringa)
                        for (i in 2:length(cenop))
                        {
                         stringa <- sprintf("    %2d) %s\n",(i-1),cenop[i])
                         cat(stringa)
                        }
                       }
                      }
                     }
         )


######### SHOW ######### SHOW ######### SHOW ######### SHOW ######### SHOW ######### SHOW ######### SHOW ######### SHOW #########

## For Angle class
#
setMethod(
          f="show",
          signature="Angle",
          definition=function(object)
                     {
                      if (length(object@rad_flag) != 0 & length(object@rad_flag) != 0)
                      {
                       if (!object@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.2f degrees. ***\n",object@ang)
                       if (object@rad_flag) stringa <- sprintf("*** This is an object of class Angle. Its value is %6.4f radians. ***\n",object@ang)
                       cat(stringa)
                      }
                      if (length(object@rad_flag) == 0 & length(object@ang) != 0)
                      {
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is %6.2f, but it is unclear",object@ang)
                       stringa2 <- sprintf("whether this value is in degrees or radians. ***\n")
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(object@rad_flag) != 0 & length(object@ang) == 0)
                      {
                       if (object@rad_flag) fchar <- "radians"
                       if (!object@rad_flag) fchar <- "degrees"
                       stringa1 <- sprintf("*** This is an object of class Angle. Its numerical value is not known, but")
                       stringa2 <- sprintf("it is measured in %s. ***\n",fchar)
                       cat(paste(stringa1,stringa2))
                      }
                      if (length(object@rad_flag) == 0 & length(object@ang) == 0)
                      {
                       stringa <- sprintf("*** This is an object of class Angle. Its slots are both empty. ***\n")
                       cat(stringa)
                      }
                     }
         )

## For UnitCell class
#
setMethod(
          f="show",
          signature="UnitCell",
          definition=function(object)
                     {
                      cat("*** This is a UnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(object@a) == 0) cat(" Slot for parameter a is empty.\n")
                      if (length(object@a) != 0) cat(sprintf(" Parameter a = %8.3f\n",object@a))
                      if (length(object@b) == 0) cat(" Slot for parameter b is empty.\n")
                      if (length(object@b) != 0) cat(sprintf(" Parameter b = %8.3f\n",object@b))
                      if (length(object@c) == 0) cat(" Slot for parameter c is empty.\n")
                      if (length(object@c) != 0) cat(sprintf(" Parameter c = %8.3f\n",object@c))
                      if (length(object@alpha@ang) == 0) cat(" Slot for parameter alpha is empty.\n")
                      if (length(object@alpha@ang) != 0) cat(sprintf(" Parameter alpha = %6.2f degrees\n",object@alpha@ang))
                      if (length(object@beta@ang) == 0) cat(" Slot for parameter beta is empty.\n")
                      if (length(object@beta@ang) != 0) cat(sprintf(" Parameter beta = %6.2f degrees\n",object@beta@ang))
                      if (length(object@gamma@ang) == 0) cat(" Slot for parameter gamma is empty.\n")
                      if (length(object@gamma@ang) != 0) cat(sprintf(" Parameter gamma = %6.2f degrees\n",object@gamma@ang))
                     }
         )

## For ReciprocalUnitCell class
#
setMethod(
          f="show",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      cat("*** This is a ReciprocalUnitCell object. ***\n")

                      # List all values. Consider what to write if some are empty
                      if (length(object@ar) == 0) cat(" Slot for parameter ar is empty.\n")
                      if (length(object@ar) != 0) cat(sprintf(" Parameter ar = %8.3f\n",object@ar))
                      if (length(object@br) == 0) cat(" Slot for parameter br is empty.\n")
                      if (length(object@br) != 0) cat(sprintf(" Parameter br = %8.3f\n",object@br))
                      if (length(object@cr) == 0) cat(" Slot for parameter cr is empty.\n")
                      if (length(object@cr) != 0) cat(sprintf(" Parameter cr = %8.3f\n",object@cr))
                      if (length(object@alphar@ang) == 0) cat(" Slot for parameter alphar is empty.\n")
                      if (length(object@alphar@ang) != 0) cat(sprintf(" Parameter alphar = %6.2f degrees\n",object@alphar@ang))
                      if (length(object@betar@ang) == 0) cat(" Slot for parameter betar is empty.\n")
                      if (length(object@betar@ang) != 0) cat(sprintf(" Parameter betar = %6.2f degrees\n",object@betar@ang))
                      if (length(object@gammar@ang) == 0) cat(" Slot for parameter gammar is empty.\n")
                      if (length(object@gammar@ang) != 0) cat(sprintf(" Parameter gammar = %6.2f degrees\n",object@gammar@ang))
                     }
         )

## For BravaisType class
#
setMethod(
          f="show",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Initial message
                      cat("*** This is a BravaisType object. ***\n")

                      # If object is empty print message
                      if (length(object@bl) == 0) cat(" The slot for parameter bl is empty.\n")

                      # Extract bl character
                      if (length(object@bl) != 0)
                      {
                       bl <- object@bl
                       cat(paste(" The slot for parameter bl contains character ",bl,".\n",sep=""))
                      }
                     }
         )

## For Lattice class
#
setMethod(
          f="show",
          signature="Lattice",
          definition=function(object)
                     {
                      # Initial message
                      cat("*** This is a Lattice object. ***\n")

                      # Print empty or cell values
                      cat(" The unit cell has the following parameters:\n")
                      if (length(object@cell@a) == 0) cat(" a is empty.\n")
                      if (length(object@cell@a) != 0) cat(sprintf(" a = %7.2f\n",object@cell@a))
                      if (length(object@cell@b) == 0) cat(" b is empty.\n")
                      if (length(object@cell@b) != 0) cat(sprintf(" b = %7.2f\n",object@cell@b))
                      if (length(object@cell@c) == 0) cat(" c is empty.\n")
                      if (length(object@cell@c) != 0) cat(sprintf(" c = %7.2f\n",object@cell@c))
                      if (length(object@cell@alpha@ang) == 0) cat(" alpha is empty.\n")
                      if (length(object@cell@alpha@ang) != 0) cat(sprintf(" alpha = %6.2f degrees.\n",object@cell@alpha@ang))
                      if (length(object@cell@beta@ang) == 0) cat(" beta is empty.\n")
                      if (length(object@cell@beta@ang) != 0) cat(sprintf(" beta = %6.2f degrees.\n",object@cell@beta@ang))
                      if (length(object@cell@gamma@ang) == 0) cat(" gamma is empty.\n")
                      if (length(object@cell@gamma@ang) != 0) cat(sprintf(" gamma = %6.2f degrees.\n",object@cell@gamma@ang))
                      if (length(object@bl@bl) == 0) cat(" Information on the lattice is not vailable.\n")

                      # Extract lattice information (if lattice has been defined)
                      if (length(object@bl@bl) != 0)
                      {
                       latt_info <- .extractLatticeStuff(object)

                       # Print
                       stringa <- sprintf("\n The Bravais lattice is %s",latt_info[1])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal family is %s",latt_info[2])
                       cat(stringa)
                       stringa <- sprintf("\n The crystal system is %s",latt_info[3])
                       cat(stringa)
                       stringa <- sprintf("\n The lattice system is %s\n",latt_info[4])
                       cat(stringa)
                      }
                     }
         )

## For Symmetry class
#
setMethod(
          f="show",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Initial message
                      cat("*** This is a Symmetry object ***\n")

                      # Empty object
                      if (length(object@sym_xHM) == 0)
                      {
                       stringa <- sprintf(" The slot for its extended Herman-Maguin symbol is empty.\n")
                       cat(stringa)
                      }

                      if (length(object@sym_xHM) != 0)
                      {
                       # Find space group number
                       lista <- .translate_SG(value=object@sym_xHM,SG_in="xHM",SG_out="number")
                       sgn <- lista$msg

                       # Find crystal system
                       cr_sys <- .crystal_system(sgn)

                       # Find symmetry and centring operations
                       lista <- .syminfo_to_op_xyz_list(object@sym_xHM)
                       xyz <- lista[[1]]
                       cenop <- lista[[2]]

                       # Print
                       stringa <- sprintf("    The space group extended Hermann-Maguin symbol is %s\n",object@sym_xHM)
                       cat(stringa)
                       stringa <- sprintf("    The space group number from the International Tables is %d\n",sgn)
                       cat(stringa)
                       stringa <- sprintf("    This space group belongs to the %s crystal system\n",cr_sys)
                       cat(stringa)
                       stringa <- sprintf("    This space group includes the following symmetry operations:\n")
                       cat(stringa)
                       for (i in 1:length(xyz))
                       {
                        stringa <- sprintf("    %2d) %s\n",i,xyz[i])
                        cat(stringa)
                       }
                       if (length(cenop) == 1)
                       {
                        stringa <- sprintf("    There are no centring operators associated with this group\n")
                        cat(stringa)
                       }
                       if (length(cenop) != 1)
                       {
                        stringa <- sprintf("    This space group is associated with a centred cell. The centring operators are:\n")
                        cat(stringa)
                        for (i in 2:length(cenop))
                        {
                         stringa <- sprintf("    %2d) %s\n",(i-1),cenop[i])
                         cat(stringa)
                        }
                       }
                      }
                     }
         )


######### GETFIELDS ######### GETFIELDS ######### GETFIELDS ######### GETFIELDS ######### GETFIELDS ######### GETFIELDS #########

## For Angle class
#
setMethod(
          f="getFields",
          signature="Angle",
          definition=function(object)
                     {
                      lista <- list(ang=object@ang,rad_flag=object@rad_flag)

                      return(lista)
                     }
         )

## For UnitCell class
#
setMethod(
          f="getFields",
          signature="UnitCell",
          definition=function(object)
                     {
                      # Extract slots
                      a <- object@a
                      b <- object@b
                      c <- object@c
                      alpha <- object@alpha
                      beta <- object@beta
                      gamma <- object@gamma
                      lista <- list(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

                      return(lista)
                     }
         )

## For ReciprocalUnitCell class
#
setMethod(
          f="getFields",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # Extract slots
                      ar <- object@ar
                      br <- object@br
                      cr <- object@cr
                      alphar <- object@alphar
                      betar <- object@betar
                      gammar <- object@gammar
                      lista <- list(ar=ar,br=br,cr=cr,alphar=alphar,betar=betar,gammar=gammar)

                      return(lista)
                     }
         )

## For BravaisType class
#
setMethod(
          f="getFields",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Extract slots
                      bl <- object@bl
                      lista <- list(bl=bl)

                      return(lista)
                     }
         )

## For Lattice class
#
setMethod(
          f="getFields",
          signature="Lattice",
          definition=function(object)
                     {
                      # Extract slots
                      cell <- object@cell
                      bl <- object@bl
                      lista <- list(cell=cell,bl=bl)

                      return(lista)
                     }
         )

## For Symmetry class
#
setMethod(
          f="getFields",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Extract slot
                      sym_xHM <- object@sym_xHM
                      lista <- list(sym_xHM=sym_xHM)

                      return(lista)
                     }
         )


######### GETPARAMETERS ######### GETPARAMETERS ######### GETPARAMETERS ######### GETPARAMETERS ######### GETPARAMETERS #########

## For Angle class
#
setMethod(
          f="getParameters",
          signature="Angle",
          definition=function(object)
                     {
                      # For this object the parameters are simply the radians flag and the angle value
                      lista <- list(ang=object@ang,rad_flag=object@rad_flag)

                      return(lista)
                     }
         )

## For UnitCell class
#
setMethod(
          f="getParameters",
          signature="UnitCell",
          definition=function(object)
                     {
                      # Extract parameters
                      a <- object@a
                      b <- object@b
                      c <- object@c
                      alpha <- object@alpha@ang
                      beta <- object@beta@ang
                      gamma <- object@gamma@ang
                      lista <- list(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

                      return(lista)
                     }
         )

## For ReciprocalUnitCell class
#
setMethod(
          f="getParameters",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # Extract parameters
                      ar <- object@ar
                      br <- object@br
                      cr <- object@cr
                      alphar <- object@alphar@ang
                      betar <- object@betar@ang
                      gammar <- object@gammar@ang
                      lista <- list(ar=ar,br=br,cr=cr,alphar=alphar,betar=betar,gammar=gammar)

                      return(lista)
                     }
         )

## For BravaisType class
#
setMethod(
          f="getParameters",
          signature="BravaisType",
          definition=function(object)
                     {
                      # Build list to be returned
                      lista <- list(bl=object@bl)

                      return(lista)
                     }
         )

## For Lattice class
#
setMethod(
          f="getParameters",
          signature="Lattice",
          definition=function(object)
                     {
                      # Build list to be returned
                      lista <- list(a=object@cell@a,b=object@cell@b,c=object@cell@c,
                                    alpha=object@cell@alpha@ang,beta=object@cell@beta@ang,gamma=object@cell@gamma@ang,
                                    bl=object@bl@bl)

                      return(lista)
                     }
         )

## For Symmetry class
#
setMethod(
          f="getParameters",
          signature="Symmetry",
          definition=function(object)
                     {
                      # Build list to be returned
                      lista <- list(sym_xHM=object@sym_xHM)

                      return(lista)
                     }
         )


######### GETDATA ######### GETDATA ######### GETDATA ######### GETDATA ######### GETDATA ######### GETDATA #########

## For Angle class
#
setMethod(
          f="getData",
          signature="Angle",
          definition=function(object)
                     {
                      # This class does not contain data and, thus, it returns NULL

                      return(NULL)
                     }
         )

## For UnitCell class
#
setMethod(
          f="getData",
          signature="UnitCell",
          definition=function(object)
                     {
                      # No data are contained in an object of class UnitCell, so this function returns NULL

                      return(NULL)
                     }
         )

## For ReciprocalUnitCell class
#
setMethod(
          f="getData",
          signature="ReciprocalUnitCell",
          definition=function(object)
                     {
                      # No data are contained in an object of class ReciprocalUnitCell, so this function returns NULL

                      return(NULL)
                     }
         )

## For BravaisType class
#
setMethod(
          f="getData",
          signature="BravaisType",
          definition=function(object)
                     {
                      # No data are contained in an object of class BravaisType, so this function returns NULL

                      return(NULL)
                     }
         ) 

## For Lattice class
#
setMethod(
          f="getData",
          signature="Lattice",
          definition=function(object)
                     {
                      # No data are contained in an object of class Lattice, so this function returns NULL

                      return(NULL)
                     }
         )

## For Symmetry class
#
setMethod(
          f="getData",
          signature="Symmetry",
          definition=function(object)
                     {
                      # No data are contained in an object of class Symmetry

                      return(NULL)
                     }
         )


######### EXTRACT[] ######### EXTRACT[] ######### EXTRACT[] ######### EXTRACT[] ######### EXTRACT[] ######### EXTRACT[] #########

## For Angle class
#
setMethod(
          f="[",
          signature="Angle",
          definition=function(x,i,j,drop)
                     {
                      # One or two fields?
                      n <- length(i)

                      # If n < 1 or n > 2 stop
                      if (n < 1 | n > 2) stop("This object includes only two values in its slots")

                      # Valid slots are: ang (numeric), rad_flag (logical)
                      flist <- list()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ang") flist <- c(flist,list(ang=x@ang))
                       if (i[ii] == 1) flist <- c(flist,list(ang=x@ang))
                       if (i[ii] == "rad_flag") flist <- c(flist,list(rad_flag=x@rad_flag))
                       if (i[ii] == 2) flist <- c(flist,list(rad_flag=x@rad_flag))
                       if (i[ii] != "ang" & i[ii] != "rad_flag" &
                           i[ii] != 1 & i[ii] != 2) stop("Slot is not included in object range")
                      }

                      return(flist)
                     }
         )

## For UnitCell class
#
setMethod(
          f="[",
          signature="UnitCell",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") flist <- c(flist,list(a=x@a))
                       if (i[ii] == "b") flist <- c(flist,list(b=x@b))
                       if (i[ii] == "c") flist <- c(flist,list(c=x@c))
                       if (i[ii] == "alpha") flist <- c(flist,list(alpha=x@alpha@ang))
                       if (i[ii] == "beta") flist <- c(flist,list(beta=x@beta@ang))
                       if (i[ii] == "gamma") flist <- c(flist,list(gamma=x@gamma@ang))
                       if (i[ii] == 1) flist <- c(flist,list(a=x@a))
                       if (i[ii] == 2) flist <- c(flist,list(b=x@b))
                       if (i[ii] == 3) flist <- c(flist,list(c=x@c))
                       if (i[ii] == 4) flist <- c(flist,list(alpha=x@alpha@ang))
                       if (i[ii] == 5) flist <- c(flist,list(beta=x@beta@ang))
                       if (i[ii] == 6) flist <- c(flist,list(gamma=x@gamma@ang))
                       if (i[ii] != "a" & i[ii] != "b" & i[ii] != "c" &
                           i[ii] != "alpha" & i[ii] != "beta" & i[ii] != "gamma" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6)
                           stop("Slot not included in object range")
                      }

                      return(flist)
                     }
         )

## For ReciprocalUnitCell class
#
setMethod(
          f="[",
          signature="ReciprocalUnitCell",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ar") flist <- c(flist,list(ar=x@ar))
                       if (i[ii] == "br") flist <- c(flist,list(br=x@br))
                       if (i[ii] == "cr") flist <- c(flist,list(cr=x@cr))
                       if (i[ii] == "alphar") flist <- c(flist,list(alphar=x@alphar@ang))
                       if (i[ii] == "betar") flist <- c(flist,list(betar=x@betar@ang))
                       if (i[ii] == "gammar") flist <- c(flist,list(gammar=x@gammar@ang))
                       if (i[ii] == 1) flist <- c(flist,list(ar=x@ar))
                       if (i[ii] == 2) flist <- c(flist,list(br=x@br))
                       if (i[ii] == 3) flist <- c(flist,list(cr=x@cr))
                       if (i[ii] == 4) flist <- c(flist,list(alphar=x@alphar@ang))
                       if (i[ii] == 5) flist <- c(flist,list(betar=x@betar@ang))
                       if (i[ii] == 6) flist <- c(flist,list(gammar=x@gammar@ang))
                       if (i[ii] != "ar" & i[ii] != "br" & i[ii] != "cr" &
                           i[ii] != "alphar" & i[ii] != "betar" & i[ii] != "gammar" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6)
                           stop("Slot not included in object range")
                      }

                      return(flist)
                     }
         )

## For BravaisType class
#
setMethod(
          f="[",
          signature="BravaisType",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Only 1 parameter for object of class "BravaisType")
                      if (n < 1 | n > 1) stop("This object does not include so many parameters")

                      # Create list with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "bl") flist <- c(flist,list(bl=x@bl))
                       if (i[ii] == 1) flist <- c(flist,list(bl=x@bl))
                       if (i[ii] != "bl" & i[ii] != 1) stop("Parameter not included in list for this object")
                      }

                      return(flist)
                     }
         )

## For Lattice class
#
setMethod(
          f="[",
          signature="Lattice",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Only 7 parameters for object of class "Lattice")
                      if (n < 1 | n > 7) stop("This object does not include so many parameters")

                      # Create vector with required content in it
                      flist <- c()
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") flist <- c(flist,list(a=x@cell@a))
                       if (i[ii] == 1) flist <- c(flist,list(a=x@cell@a))
                       if (i[ii] == "b") flist <- c(flist,list(b=x@cell@b))
                       if (i[ii] == 2) flist <- c(flist,list(b=x@cell@b))
                       if (i[ii] == "c") flist <- c(flist,list(c=x@cell@c))
                       if (i[ii] == 3) flist <- c(flist,list(c=x@cell@c))
                       if (i[ii] == "alpha") flist <- c(flist,list(alpha=x@cell@alpha@ang))
                       if (i[ii] == 4) flist <- c(flist,list(alpha=x@cell@alpha@ang))
                       if (i[ii] == "beta") flist <- c(flist,list(beta=x@cell@beta@ang))
                       if (i[ii] == 5) flist <- c(flist,list(beta=x@cell@beta@ang))
                       if (i[ii] == "gamma") flist <- c(flist,list(gamma=x@cell@gamma@ang))
                       if (i[ii] == 6) flist <- c(flist,list(gamma=x@cell@gamma@ang))
                       if (i[ii] == "bl") flist <- c(flist,list(bl=x@bl@bl))
                       if (i[ii] == 7) flist <- c(flist,list(bl=x@bl@bl))
                       if (i[ii] != "a" & i[ii] != 1 &
                           i[ii] != "b" & i[ii] != 2 &
                           i[ii] != "c" & i[ii] != 3 &
                           i[ii] != "alpha" & i[ii] != 4 &
                           i[ii] != "beta" & i[ii] != 5 &
                           i[ii] != "gamma" & i[ii] != 6 &
                           i[ii] != "bl" & i[ii] != 7
                          ) stop("Parameter not included in list for this object")
                      }

                      return(flist)
                     }
         )

## For Symmetry class
#
setMethod(
          f="[",
          signature="Symmetry",
          definition=function(x,i,j,drop)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Only 1 parameter for object of class Symmetry
                      if (n != 1) stop("Object of class Symmetry has only 1 parameter")

                      # Create list with required content in it
                      if (i[1] == "sym_xHM") flist <- list(sym_xHM=x@sym_xHM)
                      if (i[1] == 1) flist <- list(sym_xHM=x@sym_xHM)
                      if (i[1] != "sym_xHM" & i[1] != 1) stop("Parameter not included in list for this object")

                      return(flist)
                     }
         )


######### REPLACE[] ######### REPLACE[] ######### REPLACE[] ######### REPLACE[] ######### REPLACE[] ######### REPLACE[] #########

## For Angle class
#
setReplaceMethod(
          f="[",
          signature="Angle",
          definition=function(x,i,j,value)
                     {
                      # One or two slots?
                      n <- length(i)

                      # Only one or two slots accepted to replace values
                      if (n != 1 & n != 2) stop("This object of class Angle has only two slots")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Only one value is required to change
                      if (n == 1)
                      {
                       if (i == "ang") x@ang <- value
                       if (i == 1) x@ang <- value
                       if (i == "rad_flag")
                       {
                        if (!is.logical(value)) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value
                       }
                       if (i == 2)
                       {
                        if (!is.logical(value)) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value
                       }
                       if (i != 1 & i != 2 & i != "ang" & i != "rad_flag") stop("Slot to be replaced is not included in object range")
                      }

                      # Two values required to change
                      if (n == 2)
                      {
                       # Check that vlues are included in a list
                       if (!is.list(value)) stop("When replacing more than one value for object of class Angle, use a list")

                       if (i[1] == "ang") x@ang <- value[[1]]
                       if (i[1] == 1) x@ang <- value[[1]]
                       if (i[2] == "rad_flag")
                       {
                        if (!is.logical(value[[2]])) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value[[2]]
                       }
                       if (i[2] == 2)
                       {
                        if (!is.logical(value[[2]])) stop("Only logical values can replace slot rad_flag of Angle class")
                        x@rad_flag <- value[[2]]
                       }
                       if (i[1] != 1 & i[1] != 2 & i[1] != "ang" & i[1] != "rad_flag") stop("Slot to be replaced is not included in object range")
                       if (i[2] != 1 & i[2] != 2 & i[2] != "ang" & i[2] != "rad_flag") stop("Slot to be replaced is not included in object range")
                      }

                      # Check changed object is correct
                      validObject(x)

                      return(x)
                     }
         )

## For UnitCell class
#
setReplaceMethod(
          f="[",
          signature="UnitCell",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # If n > 1 and value is not a list, stop
                      if (n > 1 & !is.list(value)) stop("For multiple replacement, replacing values need to be included in a list")

                      # Valid slots are: a (numeric), b (numeric), c (numeric), alpha (Angle), beta (Angle) gamma (Angle)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") x@a <- value[[ii]]
                       if (i[ii] == 1) x@a <- value[[ii]]
                       if (i[ii] == "b") x@b <- value[[ii]]
                       if (i[ii] == 2) x@b <- value[[ii]]
                       if (i[ii] == "c") x@c <- value[[ii]]
                       if (i[ii] == 3) x@c <- value[[ii]]
                       if (i[ii] == "alpha") x@alpha <- angle(value[[ii]])
                       if (i[ii] == 4) x@alpha <- angle(value[[ii]])
                       if (i[ii] == "beta") x@beta <- angle(value[[ii]])
                       if (i[ii] == 5) x@beta <- angle(value[[ii]])
                       if (i[ii] == "gamma") x@gamma <- angle(value[[ii]])
                       if (i[ii] == 6) x@gamma <- angle(value[[ii]])
                       if (i[ii] != "a" & i[ii] != "b" & i[ii] != "c" &
                           i[ii] != "alpha" & i[ii] != "beta" & i[ii] != "gamma" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      if (validObject(x)) return(x)
                     }
         )

## For ReciprocalUnitCell class
#
setReplaceMethod(
          f="[",
          signature="ReciprocalUnitCell",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Stop if n < 1 or n > 6
                      if (n < 1 | n > 6) stop("This object does not contain that many slots")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # If n > 1 and value is not a list, stop
                      if (n > 1 & !is.list(value)) stop("For multiple replacement, replacing values need to be included in a list")

                      # Valid slots are: ar (numeric), br (numeric), cr (numeric), alphar (Angle), betar (Angle) gammar (Angle)
                      # Remember: Angle ang values are supposed to be always in degrees
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "ar") x@ar <- value[[ii]]
                       if (i[ii] == 1) x@ar <- value[[ii]]
                       if (i[ii] == "br") x@br <- value[[ii]]
                       if (i[ii] == 2) x@br <- value[[ii]]
                       if (i[ii] == "cr") x@cr <- value[[ii]]
                       if (i[ii] == 3) x@cr <- value[[ii]]
                       if (i[ii] == "alphar") x@alphar <- angle(value[[ii]])
                       if (i[ii] == 4) x@alphar <- angle(value[[ii]])
                       if (i[ii] == "betar") x@betar <- angle(value[[ii]])
                       if (i[ii] == 5) x@betar <- angle(value[[ii]])
                       if (i[ii] == "gammar") x@gammar <- angle(value[[ii]])
                       if (i[ii] == 6) x@gammar <- angle(value[[ii]])
                       if (i[ii] != "ar" & i[ii] != "br" & i[ii] != "cr" &
                           i[ii] != "alphar" & i[ii] != "betar" & i[ii] != "gammar" &
                           i[ii] != 1 & i[ii] != 2 & i[ii] != 3 & i[ii] != 4 & i[ii] != 5 & i[ii] != 6) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      # Check changed object is correct
                      if (validObject(x)) return(x)
                     }
         )

## For BravaisType class
#
setReplaceMethod(
          f="[",
          signature="BravaisType",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Only 1 parameter for object of class "BravaisType")
                      if (n < 1 | n > 1) stop("This object does not include so many parameters")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Valid slots are: bl (character)
                      wflag <- FALSE
                      if (is.list(value)) value <- value[[1]]
                      if (i == "bl") x@bl <- value
                      if (i == 1) x@bl <- value
                      if (i != "bl" & i != 1) wflag <- TRUE
                      if (wflag) warning("Requested slots does not exist")

                      if (validObject(x)) return(x)
                     }
         )

## For Lattice class
#
setReplaceMethod(
          f="[",
          signature="Lattice",
          definition=function(x,i,j,value)
                     {
                      # Find out length of i
                      n <- length(i)

                      # Only 7 parameters for object of class "Lattice")
                      if (n < 1 | n > 7) stop("This object does not include so many parameters")

                      # If length of i is different from length of value stop
                      if (n != length(value)) stop("Number of items to be replaced does not match number of items replacing them.")

                      # Check if input is a list
                      if (n > 1 & !is.list(value)) stop("For multiple replacement, values need to be included in a list")

                      # Parameters types are: a (numeric), b (numeric), c (numeric), alpha (numeric), beta (numeric), gamma (numeric),
                      # bl (character)
                      wflag <- FALSE
                      for (ii in 1:n)
                      {
                       if (i[ii] == "a") x@cell@a <- value[[ii]]
                       if (i[ii] == 1) x@cell@a <- value[[ii]]
                       if (i[ii] == "b") x@cell@b <- value[[ii]]
                       if (i[ii] == 2) x@cell@b <- value[[ii]]
                       if (i[ii] == "c") x@cell@c <- value[[ii]]
                       if (i[ii] == 3) x@cell@c <- value[[ii]]
                       if (i[ii] == "alpha") x@cell@alpha@ang <- value[[ii]]
                       if (i[ii] == 4) x@cell@alpha@ang <- value[[ii]]
                       if (i[ii] == "beta") x@cell@beta@ang <- value[[ii]]
                       if (i[ii] == 5) x@cell@beta@ang <- value[[ii]]
                       if (i[ii] == "gamma") x@cell@gamma@ang <- value[[ii]]
                       if (i[ii] == 6) x@cell@gamma@ang <- value[[ii]]
                       if (i[ii] == "bl") x@bl@bl <- value[[ii]]
                       if (i[ii] == 7) x@bl@bl <- value[[ii]]
                       if (i[ii] != "a" & i[ii] != 1 &
                           i[ii] != "b" & i[ii] != 2 &
                           i[ii] != "c" & i[ii] != 3 &
                           i[ii] != "alpha" & i[ii] != 4 &
                           i[ii] != "beta" & i[ii] != 5 &
                           i[ii] != "gamma" & i[ii] != 6 &
                           i[ii] != "bl" & i[ii] != 7
                          ) wflag <- TRUE
                      }
                      if (wflag) warning("One or more requested slots do not exist.")

                      if (validObject(x)) return(x)
                     }
         )

## For Symmetry class
#
setReplaceMethod(
                 f="[",
                 signature="Symmetry",
                 definition=function(x,i,j,value)
                            {
                             # Find out length of i
                             n <- length(i)

                             # Only 1 parameter for object of class Symmetry
                             if (n != 1) stop("Object of class Symmetry has only 1 parameter")

                             # Parameter type is: sym_xHM (character)
                             if (i[1] == "sym_xHM")
                             {
                              if (is.list(value))
                              {
                               value[[1]] <- .findHM(value[[1]])
                               x@sym_xHM <- value[[1]]
                              }
                              if (!is.list(value))
                              {
                               value <- .findHM(value)
                               x@sym_xHM <- value
                              }
                             }
                             if (i[1] == 1)
                             {
                              if (is.list(value))
                              {
                               value[[1]] <- .findHM(value[[1]])
                               x@sym_xHM <- value[[1]]
                              }
                              if (!is.list(value))
                              {
                               value <- .findHM(value)
                               x@sym_xHM <- value
                              }
                             }
                             if (i[1] != "sym_xHM" & i[1] != 1) stop("One or more requested slots do not exist")

                             if (validObject(x)) return(x)
                            }
                )


######### CONSTRUCTORS ######### CONSTRUCTORS ######### CONSTRUCTORS ######### CONSTRUCTORS ######### CONSTRUCTORS #########

## For Angle class
#
# Both ang and rad_flag (no default assignment)
setMethod(
          f="angle",
          signature=c("numeric","logical"),
          function(ang,rad_flag,...)
          {
           # Assign values to old code's variables
           value <- ang

           # Turn value in the 0 - 180 range
           if (rad_flag)
           {
            value <- value%%(2*pi)
            if (value > pi) value <- 2*pi-value
           }
           if (!rad_flag)
           {
            value <- value%%360
            if (value > 180) value <- 360-value
           }

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=value,rad_flag=rad_flag)

           return(object)
          }
         )
# Only ang (default assignment - FALSE - for rad_flag)
setMethod(
          f="angle",
          signature=c("numeric","missing"),
          function(ang,rad_flag,...)
          {
           # Assign values to old code's variables
           value <- ang

           # Values between 0 and 180 degrees
           value <- value%%360
           if (value > 180) value <- 360-value

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=value,rad_flag=FALSE)

           return(object)
          }
         )
# Only rad_flag (default assignment - 90 - for ang)
setMethod(
          f="angle",
          signature=c("logical","missing"),
          function(ang,rad_flag,...)
          {
           # ang value is fixed at 90.0 degrees
           value <- 90.0
           rad_flag <- ang

           # If rad_flag is TRUE, turn value in radians
           if (rad_flag) value <- 0.5*pi

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=value,rad_flag=rad_flag)

           return(object)
          }
         )
# Only rad_flag,with rad_flag= declared (default assignment - 90 - for ang)
#setMethod(
#          f="angle",
#          signature=c("ANY","logical"),
#          function(ang,rad_flag,...)
#          {
#           # ang value is fixed at 90.0 degrees
#           value <- 90.0
#
#           # If rad_flag is TRUE, turn value in radians
#           if (rad_flag) value <- 0.5*pi
#
#           # Now create an instance of Angle class
#           object <- new(Class="Angle",ang=value,rad_flag=rad_flag)
#
#           return(object)
#          }
#         )
# No values (default assignment - 90 - for ang and FALSE for rad_flag)
setMethod(
          f="angle",
          signature=c("missing","missing"),
          function(ang,rad_flag,...)
          {
           # ang value is fixed at 90.0 degrees
           # rad_flag is fixed at FALSE

           # Now create an instance of Angle class
           object <- new(Class="Angle",ang=90,rad_flag=FALSE)

           return(object)
          }
         )

## For UnitCell class
#
# All 6 parameters given (no default assignment)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","numeric","numeric","numeric"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # Turn angles into Angle objects
           aa <- angle(alpha)
           bb <- angle(beta)
           cc <- angle(gamma)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Three cell sides and one angle (beta)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","numeric","missing","missing"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # Turn angles into Angle objects (entry alpha is copied into exit beta)
           bb <- angle(alpha)
           aa <- angle(90)
           cc <- angle(90)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Three cell sides and no angles (orthorombic, tetragonal or cubic)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","numeric","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Two cell sides and no angles (tetragonal or cubic)
setMethod(
          f="unitcell",
          signature=c("numeric","numeric","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=a,c=b,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic - a)
setMethod(
          f="unitcell",
          signature=c("numeric","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           b <- a
           c <- a

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# No cell sides and no angles (cubic with default side 1)
setMethod(
          f="unitcell",
          signature=c("missing","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           a <- 1
           b <- 1
           c <- 1

           # Create UnitCell object
           object <- new(Class="UnitCell",a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc)

           return(object)
          }
         )
# Values taken from file
# File covered are: MTZ, PDB
setMethod(
          f="unitcell",
          signature=c("character","missing","missing","missing","missing","missing"),
          function(a,b,c,alpha,beta,gamma,...)
          {
           # In this case a is really containing a file name
           file <- a
           ans <- file.exists(file)
           if (!ans)
           {
            stop("Character input is not an existing file name")
           }

           # Check file type (only MTZ, PDB allowed at present)
           tipo <- NULL
           if (.try_read_MTZ(file))                                      # MTZ case
           {
            # Load mtz file into a named list
            lmtz <- .readMTZ(file,messages=FALSE)
 
            # Extract cell parameters
            cpar <- lmtz$header$CELL

            # Create new object
            object <- new("UnitCell",a=cpar[1],b=cpar[2],c=cpar[3],
                          alpha=new("Angle",ang=cpar[4],rad_flag=FALSE),
                          beta=new("Angle",ang=cpar[5],rad_flag=FALSE),
                          gamma=new("Angle",ang=cpar[6],rad_flag=FALSE))
            return(object) 
           }
           if (.try_read_PDB(file))                                      # PDB case
           {
            # Load pdb file into a named list
            lpdb <- .readPDB(file)

            # Cell parameters
            cpar <- lpdb$cryst1$cell_par

            # Create new object
            object <- new("UnitCell",a=cpar[1],b=cpar[2],c=cpar[3],
                          alpha=new("Angle",ang=cpar[4],rad_flag=FALSE),
                          beta=new("Angle",ang=cpar[5],rad_flag=FALSE),
                          gamma=new("Angle",ang=cpar[6],rad_flag=FALSE))
            return(object) 
           }

           # Stop if file is not in one of allowed formats
           if (is.null(tipo)) stop("This data file is corrupted, or in a format not covered by this package")
          }
         )

## For ReciprocalUnitCell class
#
# All 6 parameters given (no default assignment)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","numeric","numeric","numeric"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # Turn angles into Angle objects
           aa <- angle(alphar)
           bb <- angle(betar)
           cc <- angle(gammar)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Three cell sides and one angle (betar)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","numeric","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # Turn angles into Angle objects (entry alphar is copied to exit betar)
           bb <- angle(alphar)
           aa <- angle(90)
           cc <- angle(90)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Three cell sides and no angles (orthorombic, tetragonal or cubic)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","numeric","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Two cell sides and no angles (tetragonal or cubic)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","numeric","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=ar,cr=br,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# One cell side and no angles (cubic)
setMethod(
          f="reciprocalunitcell",
          signature=c("numeric","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           br <- ar
           cr <- ar

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# No cell sides and no angles (cubic with default side 1)
setMethod(
          f="reciprocalunitcell",
          signature=c("missing","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # Turn angles into Angle objects
           aa <- angle(90)
           bb <- angle(90)
           cc <- angle(90)

           # Missing sides
           ar <- 1
           br <- 1
           cr <- 1

           # Create ReciprocalUnitCell object
           object <- new(Class="ReciprocalUnitCell",ar=ar,br=br,cr=cr,alphar=aa,betar=bb,gammar=cc)

           return(object)
          }
         )
# Values taken from file
# File covered are: MTZ, PDB
setMethod(
          f="reciprocalunitcell",
          signature=c("character","missing","missing","missing","missing","missing"),
          function(ar,br,cr,alphar,betar,gammar,...)
          {
           # In this case ar is really containing a file name
           file <- ar

           # Check file exists
           ans <- file.exists(file)
           if (!ans)
           {
            stop("Character input is not an existing file name")
           }

           # Check file type (only MTZ, PDB allowed at present)
           tipo <- NULL
           if (.try_read_MTZ(file))                                      # MTZ case
           {
            # Load mtz file into a named list
            lmtz <- .readMTZ(file,messages=FALSE)
 
            # Extract cell parameters
            cpar <- lmtz$header$CELL

            # Now compute the 6-parameters list
            lista <- .dcl(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])

            # Calculate all 6 unit cell parameters
            ar <- lista[[7]]
            br <- lista[[8]]
            cr <- lista[[9]]
            alphar <- atan2(lista[[10]],lista[[13]])*180/pi
            betar  <- atan2(lista[[11]],lista[[14]])*180/pi
            gammar <- atan2(lista[[12]],lista[[15]])*180/pi

            # Create new object
            object <- new("ReciprocalUnitCell",ar=ar,br=br,cr=cr,
                          alphar=new("Angle",ang=alphar,rad_flag=FALSE),
                          betar=new("Angle",ang=betar,rad_flag=FALSE),
                          gammar=new("Angle",ang=gammar,rad_flag=FALSE))
            return(object) 
           }
           if (.try_read_PDB(file))                                      # PDB case
           {
            # Load pdb file into a named list
            lpdb <- .readPDB(file)

            # Cell parameters
            cpar <- lpdb$cryst1$cell_par

            # Now compute the 6-parameters list
            lista <- .dcl(cpar[1],cpar[2],cpar[3],cpar[4],cpar[5],cpar[6])

            # Calculate all 6 unit cell parameters
            ar <- lista[[7]]
            br <- lista[[8]]
            cr <- lista[[9]]
            alphar <- atan2(lista[[10]],lista[[13]])*180/pi
            betar  <- atan2(lista[[11]],lista[[14]])*180/pi
            gammar <- atan2(lista[[12]],lista[[15]])*180/pi

            # Create new object
            object <- new("ReciprocalUnitCell",ar=ar,br=br,cr=cr,
                          alphar=new("Angle",ang=alphar,rad_flag=FALSE),
                          betar=new("Angle",ang=betar,rad_flag=FALSE),
                          gammar=new("Angle",ang=gammar,rad_flag=FALSE))
            return(object) 
           }

           # Stop if file is not in one of allowed formats
           if (is.null(tipo)) stop("This data file is corrupted, or in a format not covered by this package")
          }
         )

## For BravaisType class
#
# Create a BravaisType object starting from a string (could be Bravais Type or file name)
setMethod(
          f="bravaistype",
          signature="character",
          function(x,...)
          {
           # If character indicates file name call proper specification of generic function
           if (x != "aP" & x != "mS" & x != "mA" & x != "mB" & x != "mC" & x != "mI" & x != "mP" &
               x != "oP" & x != "oS" & x != "oA" & x != "oB" & x != "oC" & x != "oI" & x != "oF" &
               x != "tP" & x != "tI" &
               x != "hP" & x != "hR" &
               x != "cP" & x != "cI" & x != "cF")
           {
            ans <- file.exists(x)
            if (ans)
            {
             file <- x
             # Check file type (only MTZ, PDB allowed at present)
             tipo <- NULL
             if (.try_read_MTZ(file))                                      # MTZ case
             {
              # Load mtz file into a named list
              lmtz <- .readMTZ(file,message=FALSE)

              # Extract SG number
              sgn <- lmtz$header$SYMINF[[4]]
              lista <- .translate_SG(sgn)
              if (!lista$ans)
              {
               cat("*** Space Group symbol read in MTZ file header has formatting problems ***\n")
               stop(lista$msg)
              }
              SG_name <- lista$msg

              # Build bl string in two steps
              b <- .crystal_system(sgn)
              if (b == "TRICLINIC") b <- "a"
              if (b == "MONOCLINIC") b <- "m"
              if (b == "ORTHOROMBIC") b <- "o"
              if (b == "TETRAGONAL") b <- "t"
              if (b == "TRIGONAL") b <- "h"
              if (b == "HEXAGONAL") b <- "h"
              if (b == "CUBIC") b <- "c"
              l <- substr(SG_name,1,1)
              bl <- paste(b,l,sep="")

              object <- new("BravaisType",bl=bl)

              return(object)
             }
             if (.try_read_PDB(file))                                      # PDB case
             {
              # Load pdb file into a named list
              lpdb <- .readPDB(file)

              # Extract SG name in xHM format and translate it in space group number
              SG_name <- lpdb$cryst1$SG
              # There are cases where xHM in syminfo is different from xHM in the PDB file
              SG_name <- .findHM(SG_name)
              lista <- .translate_SG(SG_name,SG_in="xHM",SG_out="number")
              if (!lista$ans) stop("Extended Herman-Maguin symbol from PDB file appears to be wrongly formatted")
              sgn <- lista$msg

              # Build bl string in two steps
              b <- .crystal_system(sgn)
              if (b == "TRICLINIC") b <- "a"
              if (b == "MONOCLINIC") b <- "m"
              if (b == "ORTHOROMBIC") b <- "o"
              if (b == "TETRAGONAL") b <- "t"
              if (b == "TRIGONAL") b <- "h"
              if (b == "HEXAGONAL") b <- "h"
              if (b == "CUBIC") b <- "c"
              l <- substr(SG_name,1,1)
              bl <- paste(b,l,sep="")

              object <- new("BravaisType",bl=bl)

              return(object)
             }

             # Stop if file is not in one of allowed formats
             if (is.null(tipo)) stop("This data file is corrupted, or in a format not covered by this package")
            }
            if (!ans)
            {
             stop("Character input is neither a valid Bravais Type nor an existing file name")
            }
           }

           # Otherwise is a Bravais Type
           object <- new(Class="BravaisType",bl=x)

           return(object)
          }
         )
# Create a BravaisType object starting from nothing (default is cP)
setMethod(
          f="bravaistype",
          signature="missing",
          function(x,...)
          {
           object <- new(Class="BravaisType",bl="cP")

           return(object)
          }
         )

## For Lattice class
#
# Create a Lattice object starting from UnitCell and BravaisType
setMethod(
          f="lattice",
          signature=c("UnitCell","BravaisType","missing","missing","missing","missing"),
          function(x1,x2,x3,x4,x5,x6,...)
          {
           # No need to do anything else
           object <- new("Lattice",cell=x1,bl=x2)

           return(object)
          }
         )
# Create a Lattice object starting from BravaisType and UnitCell
setMethod(
          f="lattice",
          signature=c("BravaisType","UnitCell","missing","missing","missing","missing"),
          function(x1,x2,x3,x4,x5,x6,...)
          {
           # No need to do anything else
           object <- new("Lattice",cell=x2,bl=x1)

           return(object)
          }
         )
# Starting from UnitCell only and assigning arbitrary, but compatible, BravaisType
setMethod(
          f="lattice",
          signature=c("UnitCell","missing","missing","missing","missing","missing"),
          function(x1,x2,x3,x4,x5,x6,...)
          {
           # Extract values from cell
           a <- x1@a
           b <- x1@b
           c <- x1@c
           aa <- x1@alpha@ang
           bb <- x1@beta@ang
           cc <- x1@gamma@ang

           # Default value for bbl is triclinic
           bbl <- "aP"

           # Cubic
           if (abs(a-b) < 0.000001 & abs(a-c) < 0.000001 & abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001) bbl <- "cP"

           # Tetragonal
           if ((abs(a-b) < 0.000001 & abs(a-c) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "tP"
           if ((abs(a-c) < 0.000001 & abs(a-b) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "tP"
           if ((abs(b-c) < 0.000001 & abs(b-a) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "tP"

           # Hexagonal
           if ((abs(a-b) < 0.000001 & abs(a-c) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-120) < 0.000001)) bbl <- "hP"
           if ((abs(a-c) < 0.000001 & abs(a-b) >= 0.000001) & (abs(aa-90) < 0.000001 & abs(bb-120) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "hP"
           if ((abs(b-c) < 0.000001 & abs(b-a) >= 0.000001) & (abs(aa-120) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "hP"

           # Rombohedral
           if ((abs(a-b) < 0.000001 & abs(a-c) < 0.000001) & (abs(aa-bb) < 0.000001 & abs(aa-cc) < 0.000001)) bbl <- "hR"

           # Orthorombic
           if ((abs(a-b) > 0.000001 & abs(a-c) > 0.000001 & abs(b-c) > 0.000001) &
               (abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001)) bbl <- "oP"

           # Monoclinic
           if (abs(aa-120) > 0.000001 & abs(aa-90) > 0.000001 & abs(bb-90) < 0.000001 & abs(cc-90) < 0.000001) bbl <- "mP"
           if (abs(bb-120) > 0.000001 & abs(bb-90) > 0.000001 & abs(aa-90) < 0.000001 & abs(cc-90) < 0.000001) bbl <- "mP"
           if (abs(cc-120) > 0.000001 & abs(cc-90) > 0.000001 & abs(aa-90) < 0.000001 & abs(bb-90) < 0.000001) bbl <- "mP"

           # Create BravaisType object
           bl <- bravaistype(bbl)

           # Create Lattice object
           object <- new(Class="Lattice",cell=x1,bl=bl)

           return(object)
          }
         )
# Create a Lattice object starting from BravaisType only, and assigning default, compatible unit cell
setMethod(
          f="lattice",
          signature=c("BravaisType","missing","missing","missing","missing","missing"),
          function(x1,x2,x3,x4,x5,x6,...)
          {
           # Extract character from BravaisType object (x1 in this case, because first in the list)
           bbl <- x1@bl

           # Cubic (default 1 1 1 90 90 90)
           if (bbl == "cP" | bbl == "cI" | bbl == "cF") cella <- new(Class="UnitCell",a=1,b=1,c=1,
               alpha=new(Class="Angle",ang=90,rad_flag=FALSE),
               beta=new(Class="Angle",ang=90,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=90,rad_flag=FALSE))

           # Hexagonal
           if (bbl == "hP") cella <- new(Class="UnitCell",a=1,b=1,c=2,
               alpha=new(Class="Angle",ang=90,rad_flag=FALSE),
               beta=new(Class="Angle",ang=90,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=120,rad_flag=FALSE))

           # Rombohedral
           if (bbl == "hR") cella <- new(Class="UnitCell",a=1,b=1,c=1,
               alpha=new(Class="Angle",ang=80,rad_flag=FALSE),
               beta=new(Class="Angle",ang=80,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=80,rad_flag=FALSE))

           # Tetragonal
           if (bbl == "tP") cella <- new(Class="UnitCell",a=1,b=1,c=2,
               alpha=new(Class="Angle",ang=90,rad_flag=FALSE),
               beta=new(Class="Angle",ang=90,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=90,rad_flag=FALSE))

           # Orthorombic
           if (bbl == "oP") cella <- new(Class="UnitCell",a=1,b=2,c=3,
               alpha=new(Class="Angle",ang=90,rad_flag=FALSE),
               beta=new(Class="Angle",ang=90,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=90,rad_flag=FALSE))

           # Monoclinic
           if (bbl == "mP" | bbl == "mS" | bbl == "mI" | bbl == "mA" | bbl == "mC") cella <- new(Class="UnitCell",a=1,b=2,c=3,
               alpha=new(Class="Angle",ang=90,rad_flag=FALSE),
               beta=new(Class="Angle",ang=110,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=90,rad_flag=FALSE))
           if (bbl == "mB") cella <- new(Class="UnitCell",a=1,b=2,c=3,
               alpha=new(Class="Angle",ang=110,rad_flag=FALSE),
               beta=new(Class="Angle",ang=90,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=90,rad_flag=FALSE))

           # Triclinic
           if (bbl == "aP") cella <- new(Class="UnitCell",a=1,b=2,c=3,
               alpha=new(Class="Angle",ang=80,rad_flag=FALSE),
               beta=new(Class="Angle",ang=70,rad_flag=FALSE),
               gamma=new(Class="Angle",ang=75,rad_flag=FALSE))

           # Create lattice object
           object <- new(Class="Lattice",cell=cella,bl=x1)

           return(object)
          }
         )
