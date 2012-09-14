# Script to demonstrate symmetry applied to atoms in a unit cell
# A group of atoms with cartesian coordinates is represented by dataframe "xyz".
# Symmetry is applied to the group and result is stored in data frame "xyz_eq".

# Load lattice ans symmetry modules
source("/home/james/workR/src_cRy/lattice.R")
source("/home/james/workR/src_cRy/symmetry.R")

# Create group of atoms
# 1) helix
data <- read.table("/home/james/workR/src_cRy/emp-helix-9.pdb")
xyz <- data[,7:9]

# 2) strand
#data <- read.table("/home/james/workR/src_cRy/emp-strand-9.pdb")
#xyz <- data[,7:9]

# Decide space group
SG <- translate_SG(5)

# To decide cell...
