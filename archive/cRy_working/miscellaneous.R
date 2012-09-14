# Functions needed here and there throughout main code 

# Extract all strings between empty spaces in a string and return a list containing them
noBlanks <- function(stringa)
{
 parole <- strsplit(stringa," ",fixed=TRUE)[[1]]
 posizioni <- which(nchar(parole) != 0)

 return(parole[posizioni])
}

#
crs_to_xyz <- function(crs,prm,grid,crsstart=c(0,0,0))
{
 # Turns a grid point of a map into its fractional coordinates
 # crs = c(c,r,s)
 # prm = c(mapc,mapr,maps)
 # grid = c(nx,ny,nz)
 # crsstart = c(NCSTART,NRSTART,NSSTART)

 crs <- crs-c(1,1,1)
 ixyz <- c(crs[which(prm == 1)],crs[which(prm == 2)],crs[which(prm == 3)])
 istart <- c(crsstart[which(prm == 1)],crsstart[which(prm == 2)],crsstart[which(prm == 3)])
 xyz <- c((ixyz[1]+crsstart[1])/grid[1],(ixyz[2]+crsstart[2])/grid[2],(ixyz[3]+crsstart[3])/grid[3])

 return(xyz)
}

#
xyz_to_crs <- function(x,y,z,prm,grid,crsstart=c(0,0,0))
{
 # Turns fractional coordinates into grid point of a map (3 integers)
 # prm = c(mapc,mapr,maps)
 # grid = c(nx,ny,nz)
 # crsstart = c(NCSTART,NRSTART,NSSTART)

 ixyz <- c(floor(x*grid[1]+0.5),floor(y*grid[2]+0.5),floor(z*grid[3]+0.5))
 crs <- c(ixyz[prm[1]]+1-crsstart[1],ixyz[prm[2]]+1-crsstart[2],ixyz[prm[3]]+1-crsstart[3])

 return(crs)
}

#
connect_atoms <- function(xyz,threshold=2)
{
 # Input: xyz is a matrix or dataframe with 3 columns. These are cartesian coordinates of atoms
 #        composing a fragment.
 # Output: data frame with 2 columns: index of xyz rows indicating starting and ending atoms
 #         between which collocate bonding stick
 # This function considers connected atoms closer than a threshold distance

 # A few checks
 if (!is.matrix(xyz) & !is.data.frame(xyz)) stop("Wrong input object. It has to be either a matrix or a data frame with 3 columns")
 if (dim(xyz)[2] != 3) stop("Wrong size for input object. It needs to have 3 columns")
  
 # Compute distances between all atoms
 dd <- dist(xyz)
 mm <- as.matrix(dd)

 # Create output data frame
 conn <- data.frame(B=NA,E=NA)

 # Loop over all atoms in the fragment
 for (i in 1:dim(mm)[2])
 {
  gg <- unname(which(mm[i,] < threshold & mm[i,] != 0))
 
  # Small loop over connected atoms
  for (j in gg) conn <- rbind(conn,data.frame(B=i,E=j))
 }

 # Get rid of NA
 conn <- na.omit(conn)

 return(conn)
}

#
draw_sticks_and_balls <- function(xyz,conn,radius=0.25)
{
 # Input: xyz is a matrix or dataframe with 3 columns. These are cartesian coordinates of atoms
 #        composing a fragment.
 # Input: conn is a data frame with 2 columns, created with function "connect_atoms"
 # Output: rgl view of fragment as balls and sticks model

 # Load RGL package if not already loaded
 reans <- require(rgl,quietly=T,warn.conflicts=F)
 if (!reans) print("Package rgl not installed. Install it and start again.")

 # First draw balls (atoms)
 rgl.spheres(xyz[,1],xyz[,2],xyz[,3],radius=radius,col="white")

 # Main loop
 for (i in 1:length(conn[,1]))
 {
  vx <- c(xyz[conn$B[i],1],xyz[conn$E[i],1])
  vy <- c(xyz[conn$B[i],2],xyz[conn$E[i],2])
  vz <- c(xyz[conn$B[i],3],xyz[conn$E[i],3])
  rgl.lines(vx,vy,vz,lwd=2,col="black")
 }
}
#
# Given 3 Euler angles returns 3X3 rotation matrix
euler_to_matrix <- function(alpha,beta,gamma)
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
polar_to_matrix <- function(chi,psi,fi)
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
