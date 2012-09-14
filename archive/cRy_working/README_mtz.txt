
README file for module "mtz.R" (J. Foadi - 14/12/2009)

WARNING!!! I have problem, at present, with user loading modules or files, like syminfo.lib and symmetry.R. Please amend as needed

1) mtz.R

The functions to be used are:
  a) readMTZ
  b) readMTZHeader
  c) writeMTZ

  a) readMTZ

  USAGE: readMTZ(name_of_MTZ_file,messages=TRUE)

  This function returns a list of 3 elements. It prints also COLUMN content, unless messages=FALSE.

  The list is organized as follows:
  [[1]] contains all numeric data structured as a data frame;
  [[2]] is a list of 14 or 15 elements called:
                                        NCOL (numeric vector of size 3, the number of columns, the number of reflections and the number of batches)
                                        CELL (a numeric vector of size 6, the six cell parameters)
                                        SORT (an vector of integers of size 5, telling in respect to which columns the file has been sorted)
                                        SYMINF (a list of 6 elements, all having to do with symmetry. The first element is the number of symmetry
                                                operations; the second is the number of primitive operations; the third is a character denoting the 
                                                lattice type (e.g. "P", "I", etc); the fourth is the space group number; the fifth is the space group 
                                                name; the sixth is the point group name )
                                        RESO (a vector of two numbers, the highest and lowest resolution as (1/d)^2)
                                        NDIF number of datasets represented in the file
                                        SYMM vector of string characters. Its size depends on the number of symmetry operators. One of these strings
                                             looks similar to the following: "SYMM Y,  -X,  Z+3/4     ". The length of each string is 80 characters
                                        PROJECT is a data.frame of 2 variables, a project id ("id") and a project name ("pname"). The number of observations
                                                is NDIF
                                        CRYSTAL is a data frame of 2 variables, a crystal id ("id") and a crystal name ("cname"). The number of observations
                                                is NDIF
                                        DATASET is a data frame of 2 variables, a dataset id ("id") and a dataset name ("dname"). Th enumber of observations
                                                is NDIF
                                        DCELL is a data frame of 7 variables, a dataset id ("id"), 6 cell parameters: "a", "b", "c", "alpha", "beta", "gamma".
                                                The number of observations is NDIF
                                        DWAVEL is a data frame of 2 variables, a dataset id ("id") and a wavelength ("lambda"). Th enumber of observations is
                                                NDIF
                                        COLUMN a data frame with following variables: "labels" (e.g. H, K, L, BATCH, etc...),
                                               "types (e.g. H, H, H, B, etc...), "min", "max", "id"). The number of observations depends on the specific file.
                                        BATCH is a vector containing batch numbers. This is present only if NCOL[3] is not zero
                                        HISTORY is vector of up to 30 string characters.
  [[3]] contains the Batch Header info. It is a list containing as many elements as there are batches in the MTZ file. Obviously there will be no elements if 
        the MTZ file is not an unmerged file. Each element of this list is a list of its own, consisting of 48 elements, basically the stuff described by the 
        orientation block data for multi-record files. They are:
                                        TITLE a character string
                                        NWORDS number of words in orientation block
                                        NINTGR number of integers in these blocks (first part of block includes these counts)
                                        NREALS number of reals in these blocks
                                        IORTYP type of orientation block (for possible future use, now = 0)
                                        LBCELL 6 integers used for cell parameters refinement (refinement flags for cell dimensions)
                                        MISFLG an integer. Status of missetting angles (PHIXYZ). If = 0 PHIXYZ is not used, all orientation is in UMAT; if = 1
                                               there is only one set of missetting angles, i PHIXYZ(,1) (PHIXYZ(,2)=PHIXYZ(,2)); if = 2 there are 2 sets of
                                               missetting angles (PHIXYZ(,1) != PHIXYZ(,2))
                                        JUMPAX is an integer, 1, 2 or 3; it indicates the reciprocal axis closest to principal goniostat axis E1
                                        NCRYST is the crystal number. A crystal can contain several batches
                                        LCRFLG is an integer indicating the type of crystal mosaicity information. = 0 for isotropic, = 1 for anisotropic
                                        LDTYPE is an integer indcating the type of data. = 1 oscillation data (2D), = 2 area detector data (3D), = 3 Laue data
                                        JSCAX is the goniostat scan axis number (= 1,2,3 or = 0 for multiple axis scan)
                                        NBSCAL number of batch scales and B-factors plus their standard deviations (4 at present, BSCALE, BBFAC and their
                                               standard deviations). It is set = 0 if batch scale is unset
                                        NGONAX is the number of goniostat axis (normally 1 or 3)
                                        LBMFLG is a flag for the type of beam information. = 0 for ALAMBD, DELAMB only; = 1 for ALAMBD, DELAMB, DELCOR, DIVHD,
                                               DIVVD (other options could include white beam)
                                        NDET is the number of detectors (normally 1, but a maximum of 2)
                                        LBSETID is the dataset number
                                        INTPAD is a vector with 8 numbers. At the moment they are all 0s. They could change in the future
                                        CELL a vector of size 6, the 6 cell parameters
                                        UMAT 3X3 orientation matrix
                                        PHIXYZ 3X2 matrix describing missetting angles at beginning and end of rotation
                                        CRYDAT is a list of 4, to describe the crystal mosaicity
                                        CRYDAT$ETAD reflection width (full width in degrees), to do with mosaicity
                                        CRYDAT$ETADH horizontal reflection width (to do with mosaicity)
                                        CRYDAT$ETADV vertical reflection width (to do with mosaicity)
                                        CRYDAT$GENERIC 3X3 matrix for more generic mosaicity info
                                        DATUM datum values of goniostat axes, from which Phi is measured (degrees)
                                        PHISTT start value of phi
                                        PHIEND stop value of phi
                                        SCANAX 3 numbers giving the rotation axis in laboratory frame (not yet implemented. Only relevant if JSCAXS=0)
                                        TIME1 start time in minutes
                                        TIME2 stop time in minutes
                                        BSCALE batch scale
                                        BBFAC batch temperature factor (corresponding scale is exp(-2B(sin(theta)/lambda)^2))
                                        SDBSCL sd(Bscale)
                                        SDBFAC sd(BBfac)
                                        PHIRANGE range of phi values: typically this will be PHIEND-PHISTT, but storing this explicitly allows to 
                                                 discriminate between + or - degrees rots
                                        BATPAD vector of size 11. At the moment they are all 0s, but could change in the future
                                        E1,E2,E3 3 vectors (of size 3 each) defining the NGONAX goniostat axes (in Cambridge laboratory frame)
                                        GONPAD vector of size 12. At the moment they are all 0s, but could change in the future
                                        SOURCE vector of size 3. Idealized (ie excluding tilts) source vector (antiparallel to beam), in Cambridge laboratory frame
                                        S0 vector of size 3. Source vector (antiparallel to beam) in Cambridge laboratory frame, including tilts
                                        BEMDAT is a list of 6, to describe the beam
                                        BEMDAT$ALAMBD wavelength in angstroms
                                        BEMDAT$DELAMB dispersion DeltaLambda/Lambda
                                        BEMDAT$DELCOR correlated component of wavelength dispersion
                                        BEMDAT$DIVHD horizontal beam divergence in degrees
                                        BEMDAT$DIVVD vertical beam divergence in degrees
                                        BEMDAT$rest vector of size 20. At the moment they are all 0s, but could change in the future
                                        DX1 crystal to detector distance (in mm) for first detector
                                        THETA1 first detector tilt angle (degrees)
                                        DETLM1 2X2 matrix minimum and maximum values of first detector coordinates (pixels), along YDet and ZDet
                                        DX2 crystal to detector distance (in mm) for second detector (if it exists)
                                        THETA2 second detector (if it exists) tilt angle (degrees)
                                        DETLM2 2X2 matrix minimum and maximum values of second detector (if it exists) coordinates (pixels), along YDet and ZDet
                                        DETPAD vector of size 33. At the moment they are all 0s, but could change in the future
                                        GONLAB vector of character strings, size 3, if they have been assigned. Otherwise it is simply NULL

  b) readMTZHeader

  USAGE: readMTZHeader(name_of_MTZ_file,messages=TRUE)

  This returns a list of 14 or 15 elements, exactly the same as [[2]] of readMTZ above. It prints COLUMN content unless messages=FALSE

  c) writeMTZ

  USAGE: writeMTZ(data,name_of_MTZ_file,title="Untitled")
  data is a list of 3, containing MTZ data structure (see readMTZ)

  This function writes the content of MTZ data list to a binary MTZ file. If a title is known it needs to be replaced to "Untitled" through last argument, title.
