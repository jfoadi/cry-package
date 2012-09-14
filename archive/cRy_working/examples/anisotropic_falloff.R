# Example to represent anisotropic fall-off in scaled intensities.

# Packages needed
reans <- require(rgl,quietly=T,warn.conflicts=F)
if (!reans) stop("Package rgl not installed")
reans <- require(misc3d,quietly=T,warn.conflicts=F)
if (!reans) stop("Package misc3d not installed")

# Load cRy library
.CRY_HOME <- Sys.getenv("CRY_HOME")  # To locate where directory including all cRy stuff ("cRy_dir") is located
source(paste(.CRY_HOME,"all_load_S4.R",sep="/"))

# Local functions

#
acolori <- function(x,nlevels)
{
 # This function produces 4 vector of same length as input vector x.
 # Each vector associates one of nlevels colors to each member of x.
 # there are 4 possible colour schemes; all are included in the list returned.

 hst <- hist(x,breaks=nlevels,plot=FALSE)
 cheat <- heat.colors(nlevels)
 ctopo <- topo.colors(nlevels)
 cterr <- terrain.colors(nlevels)
 cmcol <- cm.colors(nlevels)
 v <- rep(NA,times=length(x))
 u <- rep(NA,times=length(x))
 w <- rep(NA,times=length(x))
 y <- rep(NA,times=length(x))
 for (i in 1:nlevels)
 {
  idx <- which(x >= hst$breaks[i] & x < hst$breaks[i+1])
  if (i == nlevels) idx <- which(x >= hst$breaks[i] & x <= hst$breaks[i+1])
  v[idx] <- cheat[i]
  u[idx] <- ctopo[i]
  w[idx] <- cterr[i]
  y[idx] <- cmcol[i]
 }

 return(list(v,u,w,y))
}

#
partitionCube <- function(L,n)
{
 # To partition a cube is easier than to partition a sphere,
 # as there are no volume differences. This cube is centred at
 # the origin; L is half of its side.

 sbreaks <- seq(-L,L,length=(n+1))

 return(sbreaks)
}

#
sampleReflectionsInCube <- function(refs,imean,n=20)
{
 # If X,Y,Z are not included in refs data frame, stop
 nomi <- names(refs)
 if (sum(nomi == "X" | nomi == "Y" | nomi == "Z") != 3) stop("Input data frame does not contain cartesian coordinates for reciprocal nodes")

 # Half-side of cube to fit totally inside reflections region
 L <- min(range(refs$X)[2],range(refs$Y)[2],range(refs$Z)[2])

 # Grid limits for sampling
 sbreaks <- partitionCube(L,n)
 mids <- 0.5*(sbreaks[1:n]+sbreaks[2:(n+1)])

 # Sampling
 statv <- c()
 for (i in 1:n)
 {
  for (j in 1:n)
  {
   for (k in 1:n)
   {
    idx <- which((refs$X > sbreaks[i] & refs$X < sbreaks[i+1]) &
                 (refs$Y > sbreaks[j] & refs$Y < sbreaks[j+1]) &
                 (refs$Z > sbreaks[k] & refs$Z < sbreaks[k+1]))
    if (length(idx) > 0)
    {
     statv <- rbind(statv,c(mids[i],mids[j],mids[k],mean(refs[idx,imean],na.rm=TRUE)))
    }
   }
  }
 }

 return(statv)
}

#
aveIvsRes <- function(nbin,ss,II=NULL,m=min(ss),M=max(ss))
{
 # Resize vectors if needed
 idx <- which(ss >= m & ss <= M)
 ss <- ss[idx]
 if (length(II) > 0) II <- II[idx]

 # Break points
 breaks <- seq(m,M,length=(nbin+1))

 # Histogram
 hst <- hist(ss,breaks=breaks,plot=FALSE)

 # Partition intensity in resolution bins
 if (length(II) > 0)
 {
  aveI <- c()
  sdI  <- c()
  for (i in 1:nbin)
  {
   idx <- which(ss >= breaks[i] & ss < breaks[i+1])
   if (i == nbin) idx <- which(ss >= breaks[i] & ss <= breaks[i+1])
   if (length(idx) > 0)
   {
    aveI <- c(aveI,mean(II[idx],na.rm=TRUE))
    sdI  <- c( sdI,  sd(II[idx],na.rm=TRUE))
   }
   if (length(idx) == 0)
   {
    aveI <- c(aveI,NA)
    sdI <- c(sdI,NA)
   }
  }
 }
 if (length(II) == 0)
 {
  aveI <- NULL
  sdI <- NULL
 }
 
 return(list(hst$counts,hst$mids,aveI,sdI))
}

#
#partitionSphere <- function(R,nr,nt,np)
#{
# # Partition a spheric volume of radius R in regions with same volume
#
# # Partition for r
# ii <- 1:nr
# rbreaks <- c(0,(ii/nr)^(1/3)*R)
#
# # Partition for t
# ii <- 1:nt
# tbreaks <- c(0,acos(1-(2/nt)*ii))
#
# # Partition for p
# pbreaks <- seq(0,2*pi,length=(np+1))
#
# return(list(rbreaks,tbreaks,pbreaks))
#}

######### Common operations for both examples #########

# Load data into a MergedReflection object
print("Loading mtz file ...")
#file <- paste(.CRY_HOME,"examples/insulin_scaled.mtz",sep="/")
#file <- paste(.CRY_HOME,"examples/rob_scaled.mtz",sep="/")
#file <- paste(.CRY_HOME,"examples/jcsg_scaled.mtz",sep="/")
#imean <- 4
file <- paste(.CRY_HOME,"examples/hassan_low.mtz",sep="/")
imean <- 5
refOb <- createMergedReflectionFromMTZ(file)

# Expand data to whole reciprocal space
print("Expanding reflections to cover whole reciprocal space ...")
refOb <- expandMergedReflection(refOb)

# Calculate cartesian coordinates for all reciprocal nodes (uses Cambridge convention)
print("Calculating cartesian coordinates for all reciprocal nodes ...")
refOb <- computeCartesianCoordinates(refOb)

# Extract reflection data frame from MergedReflection object
refs <- extractData(refOb)

######### Example 1 #########
# Represents anisotropic falloff as a series of encapsulated shells
# in a parametric 3d plot in spherical coordinates

# Calculate scaled intensities average in coarser cubic grid (intensities to be averaged are at column 4)
print("Calculating coarser dataset for faster interpolation ...")
coarseI <- sampleReflectionsInCube(refs,imean=imean,n=10)

# Selecting only positive values for interpolation, as a logarithm is involved
print("Selecting only positive values of averages for interpolation ...")
subs <- coarseI[coarseI[,4] > 0,]

# Preparing data for linear interpolation
print("Linear interpolation ...")
x <- subs[,1]
y <- subs[,2]
z <- subs[,3]
ff <- log(subs[,4])

# Linear interpolation
model <- lm(ff~I(x^2)+I(y^2)+I(z^2)+I(2*x*y)+I(2*x*z)+I(2*y*z))

# Check the goodness of fit (r-squared)
hownice <- summary(model)
print(paste("Goodness of fit. R-squared =",hownice$r.squared))

# Fitted surface
cfs <- model$coefficients
attributes(cfs) <- NULL
lsurface <- function(x,y,z)
{
 rr <- exp(cfs[1]+cfs[2]*x^2+cfs[3]*y^2+cfs[4]*z^2+2*cfs[5]*x*y+2*cfs[6]*x*z+2*cfs[7]*y*z)
 return(rr)
}

# Add here some quantitative measure of anisotropy

# Plot reciprocal nodes and fitted surface
highP <- lsurface(0,0,0)
lowP <- lsurface(max(abs(refs$X)),max(abs(refs$Y)),max(abs(refs$Z)))
width <- highP-lowP
levels <- c(lowP+0.75*width,lowP+0.5*width,lowP+0.25*width)  # Levels for surface contours
iran <- sample(1:length(refs[,1]),size=10000,replace=TRUE)  # Select random sample of reciprocal points for representation
ransel <- refs[iran,]
open3d()
points3d(ransel$X,ransel$Y,ransel$Z)
contour3d(lsurface,level=levels,x=seq(-1,1,length=100),y=seq(-1,1,length=100),z=seq(-1,1,length=100),
          color=c("black","orange","yellow"),alpha=c(0.8,0.5,0.3),add=TRUE)
lines3d(c(0,(1.2)*max(subs[,1])),c(0,0),c(0,0),color=2,lwd=3,add=TRUE)
lines3d(c(0,0),c(0,(1.2)*max(subs[,2])),c(0,0),color=3,lwd=3,add=TRUE)
lines3d(c(0,0),c(0,0),c(0,(1.2)*max(subs[,3])),color=4,lwd=3,add=TRUE)

######### Example 2 #########
# Plot standard 1D average-intensity vs resolution along x, y, z and
# compare them with the previous fit

# Select points around each axis in turn
idx <- which((refs$Y^2+refs$Z^2 < 0.001) & refs$X > 0)
xaxis <- refs[idx,]
idx <- which((refs$X^2+refs$Z^2 < 0.001) & refs$Y > 0)
yaxis <- refs[idx,]
idx <- which((refs$X^2+refs$Y^2 < 0.001) & refs$Z > 0)
zaxis <- refs[idx,]

# Compute averages in resolution bins and plot points and curves
x11()
tmp <- aveIvsRes(20,xaxis$X,xaxis$IMEAN)
tmpx <- seq(0,max(refs$X),length=100)
tmpy <- lsurface(tmpx,0,0)
plot(tmp[[2]],tmp[[3]],pch=16,xlab="inverse resolution (x)",ylab="average intensity")
points(tmpx,tmpy,type="l",lwd=2,col=2)
x11()
tmp <- aveIvsRes(20,yaxis$Y,yaxis$IMEAN)
tmpx <- seq(0,max(refs$X),length=100)
tmpy <- lsurface(0,tmpx,0)
plot(tmp[[2]],tmp[[3]],pch=16,xlab="inverse resolution (y)",ylab="average intensity")
points(tmpx,tmpy,type="l",lwd=2,col=3)
x11()
tmp <- aveIvsRes(20,zaxis$Z,zaxis$IMEAN)
tmpx <- seq(0,max(refs$X),length=100)
tmpy <- lsurface(0,0,tmpx)
plot(tmp[[2]],tmp[[3]],pch=16,xlab="inverse resolution (z)",ylab="average intensity")
points(tmpx,tmpy,type="l",lwd=2,col=4)

# To finish: all 3 curves overlapped in the same plot
x11()
tmpx <- seq(0,max(refs$X),length=100)
tmpy <- lsurface(tmpx,0,0)
xleg <- max(tmpx)
yleg <- max(tmpy)
legend <- c("x axis","y axis","z axis")
cols <- c(2,3,4)
plot(tmpx,tmpy,type="l",lwd=2,col=2,xlab="inverse resolution (x)",ylab="average intensity")
tmpx <- seq(0,max(refs$Y),length=100)
tmpy <- lsurface(0,tmpx,0)
points(tmpx,tmpy,type="l",lwd=2,col=3,xlab="inverse resolution (y)",ylab="average intensity")
tmpx <- seq(0,max(refs$Z),length=100)
tmpy <- lsurface(0,0,tmpx)
points(tmpx,tmpy,type="l",lwd=2,col=4,xlab="inverse resolution (z)",ylab="average intensity")
legend((0.75*xleg),(0.9*yleg),legend=legend,text.col=cols)
