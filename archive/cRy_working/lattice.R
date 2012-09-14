# For use with crystallographic packages

dcl <- function(a,b,c,aa,bb,gg)
{
# Input: cell parameters. Output: values useful to all crystallographic
# calculations. These are: 
# 1) sa = sin(alpha), sb = sin(beta), sc = sin(gamma)
# 2) ca = cos(alpha), cb = cos(beta). cc = cos(gamma)
# 3) sides of reciprocal cell: ar = a*, br = b*, cr = c*
# 4) sines of angles of reciprocal cell: sar = sin(alpha*), sbr = sin(beta*), scr = sin(gamma*)
# 5) cosines of angles of reciprocal cell: car = cos(alpha*), cbr = cos(beta*), ccr = cos(gamma*)
# 6) Volume of unit cell: V
 aa <- aa*pi/180
 bb <- bb*pi/180
 gg <- gg*pi/180
 sa <- sin(aa)
 sb <- sin(bb)
 sc <- sin(gg)
 ca <- cos(aa)
 cb <- cos(bb)
 cc <- cos(gg)
 V <- a*b*c*sqrt(1-ca^2-cb^2-cc^2+2*ca*cb*cc)
 ar <- b*c*sa/V
 br <- a*c*sb/V
 cr <- a*b*sc/V
 sar <- V/(a*b*c*sb*sc)
 sbr <- V/(a*b*c*sa*sc)
 scr <- V/(a*b*c*sa*sb)
 car <- (cb*cc-ca)/(sb*sc)
 cbr <- (ca*cc-cb)/(sa*sc)
 ccr <- (ca*cb-cc)/(sa*sb)
 l <- c(sa,sb,sc,ca,cb,cc,ar,br,cr,sar,sbr,scr,car,cbr,ccr,V)
 names(l) <- c("SIN_ALPHA","SIN_BETA","SIN_GAMMA","COS_ALPHA","COS_BETA","COS_GAMMA","A*","B*","C*","SIN_ALPHA*","SIN_BETA*","SIN_GAMMA*",
               "COS_ALPHA*","COS_BETA*","COS_GAMMA*","V")

 return(l)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (First choice in Giacovazzo's book)
triclinic_to_orthogonal_01 <- function(a,b,c,aa,bb,cc)
{
 lp <- dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(a,0,0,b*lp[6],b*lp[3],0,c*lp[5],-c*lp[2]*lp[13],1/lp[9]),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (MOSFLM choice, second choice in Giacovazzo's book)
triclinic_to_orthogonal_02 <- function(a,b,c,aa,bb,cc)
{
 lp <- dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(1/lp[7],-lp[15]/(lp[7]*lp[12]),a*lp[5],0,1/(lp[8]*lp[12]),b*lp[4],0,0,c),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

d_hkl <- function(h,k,l,a,b,c,aa,bb,cc)
{
 # Given Miller indices and cell parameters, this function returns
 # resolution corresponding to the specific Miller indices.

 aa <- aa*pi/180
 bb <- bb*pi/180
 cc <- cc*pi/180
 top <- 1-(cos(aa))^2-(cos(bb))^2-(cos(cc))^2+2*cos(aa)*cos(bb)*cos(cc)
 b1 <- h^2*(sin(aa))^2/a^2
 b2 <- k^2*(sin(bb))^2/b^2
 b3 <- l^2*(sin(cc))^2/c^2
 b4 <- 2*h*k*(cos(aa)*cos(bb)-cos(cc))/(a*b)
 b5 <- 2*h*l*(cos(aa)*cos(cc)-cos(bb))/(a*c)
 b6 <- 2*k*l*(cos(bb)*cos(cc)-cos(aa))/(b*c)
 d2 <- top/(b1+b2+b3+b4+b5+b6)
 return(sqrt(d2))
}

bravais <- function(gn)
{
 # Given the space group number (gn), returns the bravais system
 # bs = 1 TRICLINIC
 # bs = 2 MONOCLINIC
 # bs = 3 ORTHOROMBIC
 # bs = 4 TETRAGONAL
 # bs = 5 CUBIC
 # bs = 6 HEXAGONAL
 # bs = 7 TRIGONAL
 if (gn == 1 | gn == 2) bs <- 1
 if (gn == 3 | gn == 4 | gn == 5) bs <- 2
 if (gn == 16 | gn == 17 | gn == 18 | gn == 19 | gn == 20 | gn == 21 | gn == 22 | gn == 23 | gn == 24) bs <- 3
 if (gn == 75 | gn == 76 | gn == 77 | gn == 78 | gn == 79 | gn == 80 |
     gn == 89 | gn == 90 | gn == 91 | gn == 92 | gn == 93 | gn == 94 | gn == 95 | gn == 96 | gn == 97 | gn == 98) bs <- 4
 if (gn == 195 | gn == 196 | gn == 197 | gn == 198 | gn == 199 |
     gn == 207 | gn == 208 | gn == 209 | gn == 210 | gn == 211 | gn == 212 | gn == 213 | gn == 214) bs <- 5
 if (gn == 168 | gn == 169 | gn == 170 | gn == 171 | gn == 172 | gn == 173 |
     gn == 177 | gn == 178 | gn == 179 | gn == 180 | gn == 181 | gn == 182) bs <- 6
 if (gn == 143 | gn == 144 | gn == 145 | gn == 146 |
     gn == 149 | gn == 150 | gn == 151 | gn == 152 | gn == 153 | gn == 154 | gn == 155) bs <- 7 
 bs_name <- c("TRICLINIC","MONOCLINIC","ORTHOROMBIC","TETRAGONAL","CUBIC","HEXAGONAL","TRIGONAL")

 return(c(bs,bs_name[bs]))
}

# Draw unit cell using RGL
draw_unit_cell <- function(a,b,c,aa,bb,cc,ochoice=1)
{
 # a,b,c = cell sides in angstroms; aa,bb,cc = cell angles in degrees

 # Load RGL package if not already loaded
 reans <- require(rgl,quietly=T,warn.conflicts=F)
 if (!reans) print("Package rgl not installed. Install it and start again.")

 # Orthogonalization matrix
 if (ochoice == 1) M <- triclinic_to_orthogonal_01(a,b,c,aa,bb,cc)
 if (ochoice == 2) M <- triclinic_to_orthogonal_02(a,b,c,aa,bb,cc)

 # Crystal axis
 av <- c(M[1,1],M[1,2],M[1,3])
 bv <- c(M[2,1],M[2,2],M[2,3])
 cv <- c(M[3,1],M[3,2],M[3,3])

 tf <- 0.2  # A small segment near origin, along a, b, c, will be painted in black
 # Crystal side a
 rgl.lines(c(0,tf*av[1]),c(0,tf*av[2]),c(0,tf*av[3]),col="black")
 rgl.lines(c(tf*av[1],av[1]),c(tf*av[2],av[2]),c(tf*av[3],av[3]),col="red")
 rgl.lines(c(cv[1],av[1]+cv[1]),c(cv[2],av[2]+cv[2]),c(cv[3],av[3]+cv[3]),col="red")
 rgl.lines(c(bv[1],av[1]+bv[1]),c(bv[2],av[2]+bv[2]),c(bv[3],av[3]+bv[3]),col="red")
 rgl.lines(c(bv[1]+cv[1],av[1]+bv[1]+cv[1]),c(bv[2]+cv[2],av[2]+bv[2]+cv[2]),c(bv[3]+cv[3],av[3]+bv[3]+cv[3]),col="red")

 # Crystal side b
 rgl.lines(c(0,tf*bv[1]),c(0,tf*bv[2]),c(0,tf*bv[3]),col="black")
 rgl.lines(c(tf*bv[1],bv[1]),c(tf*bv[2],bv[2]),c(tf*bv[3],bv[3]),col="green")
 rgl.lines(c(av[1],av[1]+bv[1]),c(av[2],av[2]+bv[2]),c(av[3],av[3]+bv[3]),col="green")
 rgl.lines(c(cv[1],bv[1]+cv[1]),c(cv[2],bv[2]+cv[2]),c(cv[3],bv[3]+cv[3]),col="green")
 rgl.lines(c(av[1]+cv[1],av[1]+bv[1]+cv[1]),c(av[2]+cv[2],av[2]+bv[2]+cv[2]),c(av[3]+cv[3],av[3]+bv[3]+cv[3]),col="green")

 # Crystal side c
 rgl.lines(c(0,tf*cv[1]),c(0,tf*cv[2]),c(0,tf*cv[3]),col="black")
 rgl.lines(c(tf*cv[1],cv[1]),c(tf*cv[2],cv[2]),c(tf*cv[3],cv[3]),col="blue")
 rgl.lines(c(av[1],av[1]+cv[1]),c(av[2],av[2]+cv[2]),c(av[3],av[3]+cv[3]),col="blue")
 rgl.lines(c(bv[1],bv[1]+cv[1]),c(bv[2],bv[2]+cv[2]),c(bv[3],bv[3]+cv[3]),col="blue")
 rgl.lines(c(av[1]+bv[1],av[1]+bv[1]+cv[1]),c(av[2]+bv[2],av[2]+bv[2]+cv[2]),c(av[3]+bv[3],av[3]+bv[3]+cv[3]),col="blue")

 return(NULL)
}

frac_to_orth <- function(xyzf,a,b,c,aa,bb,cc,ochoice=1)
{
 # a,b,c = cell sides in angstroms; aa,bb,cc = cell angles in degrees
 # The parameter ochoice controls which convention is being used to
 # collocate the cell in an orthonormal cartesian frame.
 # ochoice = 1: X axis along a; Y axis normal to a, in the (a,b) plane;
 #              Z axis normal to X and Y (and therefore parallel to
 #              c*).
 # ochoice = 2: this is also called "Cambridge setting". The X axis is
 #              along a*; the Y axis lies in the (a*,b*) plane; the Z
 #              axis is, consequently, along c. 
 # xyzf is a vector, matrix or data frame of fractional crystal coordinates.
 # If matrix or data frame it needs to have 3 columns. 
 # This function returns a data frame with 3 columns, the cartesian orthonormal coordinates.

 # Check xyzf object is either vector, matrix or data.frame
 if (!is.vector(xyzf) & !is.matrix(xyzf) & !is.data.frame(xyzf)) stop("Input object is not of a valid type (vector, matrix or data.frame)")

 # If xyzf is a vector, turn it into a matrix
 if (is.vector(xyzf)) dim(xyzf) <- c(1,length(xyzf))

 # If xyzf has less or more than 3 columns stop
 if (dim(xyzf)[2] != 3) stop("Input object has more or less than 3 columns")

 # If xyzf is a data.frame turn it into a matrix
 if (is.data.frame(xyzf)) xyzf <- as.matrix(xyzf)

 # Orthogonalization matrix
 if (ochoice == 1) M_1 <- triclinic_to_orthogonal_01(a,b,c,aa,bb,cc)
 if (ochoice == 2) M_1 <- triclinic_to_orthogonal_02(a,b,c,aa,bb,cc)

 # Transformed coordinates
 #x <- av[1]*xf+bv[1]*yf+cv[1]*zf
 #y <- av[2]*xf+bv[2]*yf+cv[2]*zf
 #z <- av[3]*xf+bv[3]*yf+cv[3]*zf
 xyz <- xyzf%*%M_1 

 # Turn matrix back into a data frame
 xyz <- as.data.frame(xyz)
 colnames(xyz) <- c("x","y","z")

 return(xyz)
}

#
orth_to_frac <- function(xyz,a,b,c,aa,bb,cc,ochoice=1)
{
 # Does the inverse job of "frac_to_orth"

 # Check xyz object is either vector, matrix or data.frame
 if (!is.vector(xyz) & !is.matrix(xyz) & !is.data.frame(xyz)) stop("Input object is not of a valid type (vector, matrix or data.frame)")

 # If xyz is a vector, turn it into a matrix
 if (is.vector(xyz)) dim(xyz) <- c(1,length(xyz))

 # If xyz has less or more than 3 columns stop
 if (dim(xyz)[2] != 3) stop("Input object has more or less than 3 columns")

 # If xyz is a data.frame turn it into a matrix
 if (is.data.frame(xyz)) xyz <- as.matrix(xyz)

  # Orthogonalization matrix
 if (ochoice == 1) M_1 <- triclinic_to_orthogonal_01(a,b,c,aa,bb,cc)
 if (ochoice == 2) M_1 <- triclinic_to_orthogonal_02(a,b,c,aa,bb,cc)

 # Inverse matrix
 M <- solve(M_1)

 # Transform coordinates
 xyzf <- xyz%*%M

 # Turn matrix back into a data frame
 xyzf <- as.data.frame(xyzf)
 colnames(xyzf) <- c("xf","yf","zf")

 return(xyzf)
}

#
draw_map_cell <- function(map,ochoice=1)
{
 # This function is a wrapper to call function "draw_unit_cell"
 # without having to input all cell parameters when a "map"
 # structure is available.
 # Input: a map structure.
 # Input: ochoice is an integer parameter indicating choice of
 # orthogonalization system (see function "frac_to_orth")

 # Extract cell parameters
 a <- map$header$X
 b <- map$header$Y
 c <- map$header$Z
 aa <- map$header$Alpha
 bb <- map$header$Beta
 cc <- map$header$Gamma

 # Call function
 draw_unit_cell(a,b,c,aa,bb,cc,ochoice=ochoice)
}

#
map_frac_to_orth <- function(xyzf,map,ochoice=1)
{
 # Wrapper function for "frac_to_orth". Cell parameters are here
 # extracted by the map object, "map",

 # Extract cell parameters
 a <- map$header$X
 b <- map$header$Y
 c <- map$header$Z
 aa <- map$header$Alpha
 bb <- map$header$Beta
 cc <- map$header$Gamma

 # Call "frac_to_orth"
 xyz <- frac_to_orth(xyzf,a,b,c,aa,bb,cc,ochoice=ochoice)

 return(xyz)
}

#
map_orth_to_frac <- function(xyz,map,ochoice=1)
{
 # Wrapper function for "orth_to_frac". Cell parameters are here
 # extracted by the map object, "map",

 # Extract cell parameters
 a <- map$header$X
 b <- map$header$Y
 c <- map$header$Z
 aa <- map$header$Alpha
 bb <- map$header$Beta
 cc <- map$header$Gamma

 # Call "orth_to_frac"
 xyzf <- orth_to_frac(xyz,a,b,c,aa,bb,cc,ochoice=1)

 return(xyzf)
}

#
squared_resolution_coeffs <- function(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
{
 # Given cell parameters return coefficients for the expression of squared resolution:
 #
 #  s^2 = a*h^2+b*k^2+c*l^2+2*d*h*k+2*e*h*l+2*f*k*l

 ctmp <- unname(dcl(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma))
 den <- 1-ctmp[4]^2-ctmp[5]^2-ctmp[6]^2+2*ctmp[4]*ctmp[5]*ctmp[6]
 a <- ctmp[1]^2/(cell_a^2*den)
 b <- ctmp[2]^2/(cell_b^2*den)
 c <- ctmp[3]^2/(cell_c^2*den)
 d <- (ctmp[4]*ctmp[5]-ctmp[6])/(cell_a*cell_b*den)
 e <- (ctmp[4]*ctmp[6]-ctmp[5])/(cell_a*cell_c*den)
 f <- (ctmp[5]*ctmp[6]-ctmp[4])/(cell_b*cell_c*den)

 return(c(a,b,c,d,e,f))
}

#
generate_reflections_up_to_resolution <- function(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,reso)
{
 # Generate all reflections (Miller indices) with resolutions lower or equal
 # to reso.

 # Inverse, squared resolution
 ss <- (1/reso)^2

 # Coefficients of squared resolution function
 ctmp <- squared_resolution_coeffs(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
 a <- ctmp[1]
 b <- ctmp[2]
 c <- ctmp[3]
 d <- ctmp[4]
 e <- ctmp[5]
 f <- ctmp[6]

 # Find max h, k, l to generate all (h,k,l) within minimal box

 # Coefficient determinant (common to both h, and k and l)
 D <- a*(b*c-f^2)-d*(c*d-e*f)+e*(d*f-b*e)

 # Max h
 A <- (b*c-f^2)/D 
 B <- (e*f-c*d)/D
 C <- (d*f-b*e)/D
 t <- sqrt(ss/(a*A^2+b*B^2+c*C^2+2*d*A*B+2*e*A*C+2*f*B*C))
 max_h <- ceiling(A*t)

 # Max k
 A <- (e*f-c*d)/D 
 B <- (a*c-e^2)/D
 C <- (d*e-a*f)/D
 t <- sqrt(ss/(a*A^2+b*B^2+c*C^2+2*d*A*B+2*e*A*C+2*f*B*C))
 max_k <- ceiling(B*t)

 # Max l
 A <- (d*f-b*e)/D 
 B <- (d*e-a*f)/D
 C <- (a*b-d^2)/D
 t <- sqrt(ss/(a*A^2+b*B^2+c*C^2+2*d*A*B+2*e*A*C+2*f*B*C))
 max_l <- ceiling(C*t)
 #print(c(max_h,max_k,max_l))

 # Dataframe with all h, k, l
 hkl <- expand.grid(h=-max_h:max_h,k=-max_k:max_k,l=-max_l:max_l)

 # Strictly reflections with resolution lower than reso
 new_hkl <- hkl[1/d_hkl(hkl$h,hkl$k,hkl$l,cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma) <= 1/reso,]
 
 # Add inverse of resolution (s) as last column of data frame
 s <- 1/d_hkl(new_hkl$h,new_hkl$k,new_hkl$l,cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
 new_hkl <- cbind(new_hkl,data.frame(s=s))

 return(new_hkl)
}
