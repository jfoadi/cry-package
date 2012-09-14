# For use with crystallographic packages

# Return matrix to rotate crystal in a given orientation
crystal_orientation <- function(chi,psi,fi)
{
 chi <- chi*pi/180
 psi <- psi*pi/180
 fi <- fi*pi/180
 o11 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(cos(fi))^2
 o12 <- -sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o13 <- -sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o21 <- sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o22 <- cos(chi)+(1-cos(chi))*(cos(psi))^2
 o23 <- sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o31 <- sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o32 <- -sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o33 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(sin(fi))^2
 oline <- c(o11,o12,o13,o21,o22,o23,o31,o32,o33)
 Om <- matrix(oline,nrow=3,ncol=3,byrow=TRUE)
 return(Om)
}

# Product of orthogonalization and orientation matrices
orientation <- function(m_1,Om)
{
 return(m_1%*%t(Om))
}

rotation_matrix <- function(axis,w)
{
 # Rotation matrix around x, y or z axis.
 # w is in degrees. axis is a string; either "x", "y", or "z".
 c <- cos(w*pi/180)
 s <- sin(w*pi/180)

 if (axis == "x")
 {
  Rw <- matrix(c(1,0,0,0,c,s,0,-s,c),nrow=3,ncol=3)
 }
 if (axis == "y")
 {
  Rw <- matrix(c(c,0,-s,0,1,0,s,0,c),nrow=3,ncol=3)
 }
 if (axis == "z")
 {
  Rw <- matrix(c(c,s,0,-s,c,0,0,0,1),nrow=3,ncol=3)
 }
 return(Rw)
}

crystal_face <- function(p1x,p1y,p1z,p2x,p2y,p2z)
{
# Compute coefficients for plane parallel to crystal axes p1 and p2, whose
# coordinates are (p1x,p1y,p1z) and (p2x,p2y,p2z) respectively. Returns a
# unit vector, (nx,ny,nz)

 nx <- p1y*p2z-p1z*p2y
 ny <- p1z*p2x-p1x*p2z
 nz <- p1x*p2y-p1y*p2x
 n <- sqrt(nx^2+ny^2+nz^2)
 nx <- nx/n
 ny <- ny/n
 nz <- nz/n

 return(c(nx,ny,nz))
}

cell_vertex <- function(a,b,c,aa,bb,cc,ochoice=1)
{
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

 # Orthogonalization matrix
 if (ochoice == 1) M <- triclinic_to_orthogonal_01(a,b,c,aa,bb,cc)
 if (ochoice == 2) M <- triclinic_to_orthogonal_02(a,b,c,aa,bb,cc)

 # Crystal axis
 av <- c(M[1,1],M[1,2],M[1,3])
 bv <- c(M[2,1],M[2,2],M[2,3])
 cv <- c(M[3,1],M[3,2],M[3,3])

 # Cell vertices
 v1 <- c(0,0,0)
 v2 <- av
 v3 <- av+bv
 v4 <- bv
 v5 <- cv
 v6 <- av+cv
 v7 <- av+bv+cv
 v8 <- bv+cv

 return(list(v1,v2,v3,v4,v5,v6,v7,v8))
}
