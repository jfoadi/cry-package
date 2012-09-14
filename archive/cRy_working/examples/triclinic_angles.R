# Example to load a csv file of triclinic-cell parameters from all data from PDB
# and do a contour3d inside the tetrahedron of my paper on unit-cell angles' constraints
# Run this by sourcing it from an interactive R session, so you can orient the way you like and
# even save a snapshot

# Load libraries
require(rgl)
require(misc3d)

# Read csv file
data <- read.csv("triclinic_cellpars_in_PDB.csv")

# As the file contains weird data with 0 and 1 or NA for alpha, beta and gamma,
# let's get rid of them
data <- na.omit(data)
data <- data[data$Angle.Alpha > 2,]
data <- data[data$Angle.Beta > 2,]
data <- data[data$Angle.Gamma > 2,]

# Create 3D kernel out of the many points with coordinates (alpha, beta, gamma)
kernel <- kde3d(data$Angle.Alpha,data$Angle.Beta,data$Angle.Gamma)

# Plot 3 3D contour surfaces
delta <- range(kernel$d)[2]-range(kernel$d)[1]
livelli <- c(min(kernel$d)+0.6*delta,min(kernel$d)+0.1*delta,min(kernel$d)+0.02*delta)
contour3d(kernel$d,level=livelli,x=kernel$x,y=kernel$y,z=kernel$z,color=c("blue","blue","blue"),alpha=c(1,0.8,0.6))

# Build skeleton of tetrahedron
lines3d(c(0,180),c(0,0),c(0,180),color="black",lwd=3)
lines3d(c(0,0),c(0,180),c(0,180),color="black",lwd=3)
lines3d(c(0,180),c(0,180),c(0,0),color="black",lwd=3)
lines3d(c(180,0),c(0,180),c(180,180),color="black",lwd=3)
lines3d(c(180,180),c(0,180),c(180,0),color="black",lwd=3)
lines3d(c(0,180),c(180,180),c(180,0),color="black",lwd=3)

# Add the cube
cubo <- cube3d(color="blue",alpha=0.3)
cubo <- translate3d(cubo,1,1,1)
cubo <- scale3d(cubo,90,90,90)
shade3d(cubo)

# Add axes and labels
axes3d()
title3d(xlab="alpha",ylab="beta",zlab="gamma")
