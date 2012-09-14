# Functions to perform 2D and 3D FFT. Direct transform maps reals to complex numbers.
# Inverse transform maps complex to reals; thus each F(h,k) needs to equal F(-h,-k)*.
# Mostly used for crystallographic maps, but can be used to perform standard orthogonal
# FFTs if angles are 90 degrees.
# Symmetry is used when available.

# Miscellaneous functions
source(paste(CRY_HOME,"miscellaneous.R",sep="/"))

# Produces data frame with full set of Miller indices.
# Input is 2 or 3 integers indicating max frequency
# for each Miller index (2D and 3D case)
generateMiller <- function(rH)
{
 ndim <- length(rH)
 if (ndim == 2)
 {
  mill <- expand.grid(h=-rH[1]:rH[1],k=-rH[2]:rH[2])
 }
 if (ndim == 3)
 {
  mill <- expand.grid(h=-rH[1]:rH[1],k=-rH[2]:rH[2],l=-rH[3]:rH[3])
 }

 return(mill)
}

# Generic 2D FFT. Input is a 2D real matrix
dfour2D <- function(ro)
{
 # Extract matrix dimensions
 nx <- dim(ro)[1]
 ny <- dim(ro)[2]

 # 2D fft as 2 1D ffts
 tmp <- apply(ro,1,fft)/length(ro[1,])
 F_cplx <- apply(tmp,1,fft)/length(tmp[1,])
 rm(tmp)

 return(F_cplx)
}

# Given a 2D or 3D complex transform, and a matrix (dimension nX2 or nX3) of Miller indices
 #(Fourier coefficients), returns corresponding modules and arguments
calign2D <- function(F_cplx)
{
 # Check dimensions of both F_cplx and Miller are compatible
 #nM <- dim(Miller)[2]
 nF <- dim(F_cplx)
 
 return(nM)
}
