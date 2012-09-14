These are still not-finalized. But all functions are ready to be used.
Here they are:


1) File symmetry.R

   ******** translate_SG(value,SG_in="number",SG_out="xHM",setting=1) ********
     # Input: space group in one of known formats; output: space group in another
     # known format. Possible formats are:
     # 1) Space group number
     # 2) Hall symbol (e.g. ' P 2yb (z,x,y)')
     # 3) Extended Hermann-Maguin symbol (e.g. 'P 1 1 21')
     # If more than one setting is implied in an unbiguous way in the input value,
     # then the first setting will be selected by default for the output value, unless
     # argument "setting" is set to another value
     # Input: SG_in is a string representing the space group format. Possible values are:
     # 1) "number"
     # 2) "ccp4"
     # 3) "Hall"
     # 4) "xHM"
     # 5) "old"
     # Input: SG_out is a string representing the space group format needed as output. Possible
     #        values are the same as those for SG_in
     # Input: setting is a number like 1,2,...If for a same symbol there are more settings, use this.
    
   ******** extract_symmetry_info(SG) ********
     # Input is spacegroup symbol in xHM format.
     # Output is a list with all symmetry information
     # for the specific SG group, as contained in syminfo.lib. 
    
   ******** syminfo_to_op_xyz_list(SG) ********
     # Input is path for CCP4 file syminfo.lib and spacegroup symbol
     # in xHM format. Output is a list of 2 vector of strings,
     # The first string describes point group and translation; the
     # second string describes cell centering
    
   ******** op_xyz_list_to_matrix_list(op_xyz_list) ********
     # Input: "syminfo_to_op_xyz_list" output, i.e. a list of 2 character vectors.
     # The first one contains symmetry operators in x,y,z format; the second one centering operators
     # in x,y,z format.
     # Returns a list consisting of 3 lists. The first list contains 3X3 point group matrices;
     # the second list contains the same number of 3X1 translation vectors. First matrix is always the
     # identity matrix, the first translation vector is always the null vector. the third
     # list consists of centering vectors; the first centering vector is always the null
     # vector. To summarize, the output looks like the following:
     # [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
     # where:
     # I            = identity    3X3 matrix
     # 0            = null        3X1 vector
     # M2,M3,...,Mn = point group 3X3 matrices
     # V2,V3,...,Cn = translation 3X1 vectors
     # C2,C3,...,Cm = centering   3X1 vectors
    
   ******** syminfo_to_matrix_list(SG) ********
     # Input: space group name in xHM format.
     # Returns a list consisting of 3 lists. The first list contains 3X3 point group matrices;
     # the second list contains the same number of 3X1 translation vectors. First matrix is always the
     # identity matrix, the first translation vector is always the null vector. the third
     # list consists of centering vectors; the first centering vector is always the null
     # vector. To summarize, the output looks like the following:
     # [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
     # where:
     # I            = identity    3X3 matrix
     # 0            = null        3X1 vector
     # M2,M3,...,Mn = point group 3X3 matrices
     # V2,V3,...,Cn = translation 3X1 vectors
     # C2,C3,...,Cm = centering   3X1 vectors
     #
     # This function is simply a wrapper for op_xyz_list_to_matrix_list.
    
   ******** equiv_xyzf(xyzf,SG,in_cell=FALSE) ********
     # Given fractional coordinates for a 3D point, returns all equivalent points
     # according to symmetry operators of space group with xHM symbol SG.
     # If "in_cell" is TRUE, all points outside cell are translated back into cell.
     # Input: xyzf can be a vector, matrix or data.frame; if a matrix or a data.frame,
     #        it needs to have 3 columns
     # Input: SG xHM format for space group
     # Input: in_cell, TRUE or FALSE
     # Output: data frame with 6 columns: xf yf zf sym c_sym orig_xyzf. "sym" indicates which
     #         operator produced the specific "xf", "yf", "zf". "c_sym" specifies the centering.
     #         "orig_xyzf" indicates to which initial point the produced point is related.
    
   ******** equiv_hkl(hkl,SG) ********
     # Given a point in reciprocal space with Miller indices h, k and l
     # and a space group with xHM symbol SG, produces a data frame of equivalent
     # points with Miller indices grouped in columns named "h","k","l","pg_sym", "orig_hkl".
     # The last columns indicate which symmetry operator was responsible for the specific
     # equivalent reflection and to which set of original Miller indices is the specific
     # index related.
     # Input: hkl can be a vector, matrix or data.frame; if a matrix or a data.frame,
     #        it needs to have 3 columns
     # Input: SG is the space group symbol in xHM format
     # Output: a data frame with 5 columns: h k l pg_sym orig_hkl
    
   ******** equiv_xyzf(xyzf,MT,in_cell=FALSE) ********
     # Given fractional coordinates for a 3D point, returns all equivalent points
     # according to symmetry operators MT. If "in_cell" is TRUE, all points outside
     # cell are translated back into cell.
    
   ******** asu_ccp4_limits(SG) ********
     # Given space group xHM symbol, returns asymmetric unit limits along
     # x, y and z (according to ccp4 convention).
    
   ******** op_xyz_to_matrix(op_xyz) ********
     # Reads in a symmetry or centering operator in character form and output it in
     # matrix or vector form. If input is a symmetry operator, output is a list of
     # a 3X3 matrix and a 3X1 vector; if input is a centering operator, output is
     # still a matrix and a vector, but the matrix is always the identity matrix.


2) File lattice.R

   ******** dcl(a,b,c,aa,bb,gg) ********
     # Input: cell parameters. Output: values useful to all crystallographic
     # calculations. These are:
     # 1) sa = sin(alpha), sb = sin(beta), sc = sin(gamma)
     # 2) ca = cos(alpha), cb = cos(beta). cc = cos(gamma)
     # 3) sides of reciprocal cell: ar = a*, br = b*, cr = c*
     # 4) sines of angles of reciprocal cell: sar = sin(alpha*), sbr = sin(beta*), scr = sin(gamma*)
     # 5) cosines of angles of reciprocal cell: car = cos(alpha*), cbr = cos(beta*), ccr = cos(gamma*)
     # 6) Volume of unit cell: V

   ******** triclinic_to_orthogonal_01(a,b,c,aa,bb,cc) ********
     # Given cell parameters, return cell orthogonalization matrix (First choice in Giacovazzo's book)

   ******** triclinic_to_orthogonal_02(a,b,c,aa,bb,cc) ********
     # Given cell parameters, return cell orthogonalization matrix (MOSFLM choice, second choice in Giacovazzo's book)

   ******** crystal_orientation(chi,psi,fi) ********
     # Return matrix to rotate crystal in a given orientation

   ******** orientation(chi,psi,fi) ********
     # Product of orthogonalization and orientation matrices

   ******** rotation_matrix(axis,w) ********
     # Rotation matrix around x, y or z axis.
     # w is in degrees. axis is a string; either "x", "y", or "z".

   ******** d_hkl(h,k,l,a,b,c,aa,bb,cc) ********
     # Given Miller indices and cell parameters, this function returns
     # resolution corresponding to the specific Miller indices.

   ******** bravais(gn) ********
     # Given the space group number (gn), returns the bravais system
     # bs = 1 TRICLINIC
     # bs = 2 MONOCLINIC
     # bs = 3 ORTHOROMBIC
     # bs = 4 TETRAGONAL
     # bs = 5 CUBIC
     # bs = 6 HEXAGONAL
     # bs = 7 TRIGONAL

   ******** crystal_face(p1x,p1y,p1z,p2x,p2y,p2z) ********
     # Compute coefficients for plane parallel to crystal axes p1 and p2, whose
     # coordinates are (p1x,p1y,p1z) and (p2x,p2y,p2z) respectively. Returns a
     # unit vector, (nx,ny,nz)

   ******** cell_vertex(a,b,c,aa,bb,cc,ochoice=1) ********
     # Returns 8 3D vectors having coordinates corresponding to the
     # 8 vertex of a unit cell. The parameter ochoice controls which
     # convention is being used to collocate the cell in an orthonormal
     # cartesian frame.
     # ochoice = 1: X axis along a; Y axis normal to a, in the (a,b) plane;
     #              Z axis normal to X and Y (and therefore parallel to
     #              c*).
     # ochoice = 2: this is also called "Cambridge setting". The X axis is
     #              along a*; the Y axis lies in the (a*,b*) plane; the Z
     #              axis is, consequently, along c.
     # a,b,c = cell sides in angstroms; aa,bb,cc = cell angles in degrees

   ******** draw_unit_cell(a,b,c,aa,bb,cc,ochoice=1) ********
     # Draw unit cell using RGL
     # a,b,c = cell sides in angstroms; aa,bb,cc = cell angles in degrees

   ******** frac_to_orth(xf,yf,zf,a,b,c,aa,bb,cc,ochoice=1) ********
     # a,b,c = cell sides in angstroms; aa,bb,cc = cell angles in degrees
     # The parameter ochoice controls which convention is being used to
     # collocate the cell in an orthonormal cartesian frame.
     # ochoice = 1: X axis along a; Y axis normal to a, in the (a,b) plane;
     #              Z axis normal to X and Y (and therefore parallel to
     #              c*).
     # ochoice = 2: this is also called "Cambridge setting". The X axis is
     #              along a*; the Y axis lies in the (a*,b*) plane; the Z
     #              axis is, consequently, along c.
     # xf,yf,zf are fractional crystal coordinates. This function returns a
     # 3D vector with cartesian orthonormal coordinates.

   ******** orth_to_frac(x,y,z,a,b,c,aa,bb,cc,ochoice=1) ********
     # Does the inverse job of "frac_to_orth"

   ******** draw_map_cell(map,ochoice=1) ********
     # This function is a wrapper to call function "draw_unit_cell"
     # without having to input all cell parameters when a "map"
     # structure is available.
     # Input: a map structure.
     # Input: ochoice is an integer parameter indicating choice of
     # orthogonalization system (see function "frac_to_orth")

   ******** map_frac_to_orth(xf,yf,zf,map,ochoice=1) ********
     # Given a 3D vector, "xyzf", in fractional coordinates, and a map object, "map",
     # this function returns a 3D vector in cartesian coordinates. The cartesian system
     # has settings decided by "ochoice" (see "frac_to_orth" function).
