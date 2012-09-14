# For use in conjunction with S4-classes based modules in cRy.
# The word "dot" has been adopted because all functions here start with a ".", so
# to be invisible.

#
# Input word "wordin" is stripped of blanks before, after or before and after the word itself.
# side can be "left","right", "both".
.strip_blanks <- function(wordin,side="both")
{
 if (side != "left" & side != "right" & side != "both") stop("Input variable side can only be 'left', 'right' or 'both'")

 if (side == "left" | side == "both")
 {
  nb <- 1
  while(substr(wordin,nb,nb) == " ") nb <- nb+1
  wordin <- substr(wordin,nb,nchar(wordin))
 }
  
 if (side == "right" | side == "both")
 {
  nb <- nchar(wordin)
  while(substr(wordin,nb,nb) == " ") nb <- nb-1
  wordin <- substr(wordin,1,nb)
 }

 return(wordin)
}
#
# Given 3 Euler angles returns 3X3 rotation matrix
.euler_to_matrix <- function(alpha,beta,gamma)
{
 # Input angles are read in degrees
 aa <- alpha*pi/180
 bb <- beta*pi/180
 cc <- gamma*pi/180
 R1 <- matrix(c(cos(cc),sin(cc),0,-sin(cc),cos(cc),0,0,0,1),nrow=3,ncol=3)
 R2 <- matrix(c(1,0,0,0,cos(bb),sin(bb),0,-sin(bb),cos(bb)),nrow=3,ncol=3)
 R3 <- matrix(c(cos(aa),sin(aa),0,-sin(aa),cos(aa),0,0,0,1),nrow=3,ncol=3)
 Rfinal <- R3%*%R2%*%R1

 return(Rfinal)
}
#
# Given 3 polar angles returns 3X3 rotation matrix
.polar_to_matrix <- function(chi,psi,fi)
{
 # Input angles are in degrees
 chi <- chi*pi/180
 psi <- psi*pi/180
 fi <- fi*pi/180

 # Matrix components
 o11 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(cos(fi))^2
 o12 <- -sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o13 <- -sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o21 <- sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o22 <- cos(chi)+(1-cos(chi))*(cos(psi))^2
 o23 <- sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o31 <- sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o32 <- -sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o33 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(sin(fi))^2

 # Matrix
 oline <- c(o11,o12,o13,o21,o22,o23,o31,o32,o33)
 Om <- matrix(oline,nrow=3,ncol=3,byrow=TRUE)

 return(Om)
}
#
# Cross product u X v in 3D. Returns a 3 components vector
"%X%" <- function(u,v)
{
 x <- u[2]*v[3]-u[3]*v[2]
 y <- u[3]*v[1]-u[1]*v[3]
 z <- u[1]*v[2]-u[2]*v[1]
 return(c(x,y,z))
}
