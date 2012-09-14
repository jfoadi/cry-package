# File to source all R modules using S4 classes load all auxiliary files 
# needed to run cRy

# Root directory (full path to src_cRy). Change this according to your settings
#.CRY_HOME <- "/Users/james/Dropbox/source_code/CRY/cRy_working"   # Use this to load modules as in e.g. "source(paste(CRY_HOME,"symmetry.R"))"
.CRY_HOME <- Sys.getenv("CRY_HOME")

# Auxiliary files. You shouldn't need editing this, as files will be located in relation to CRY_HOME
.SYMMETRY_file <- paste(.CRY_HOME,"syminfo.lib",sep="/")
.ATOMS_file <- paste(.CRY_HOME,"elements_list.dat",sep="/")

# Other needed objects
.ATOMS_data.frame <- read.table(.ATOMS_file,header=TRUE,as.is=1:2)
.ALPHABET <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","W","Z"," ")
.EXAMPLES_dir <- paste(.CRY_HOME,"examples",sep="/")

# Load all ready non-S4 modules
source(paste(.CRY_HOME,"generic_dot.R",sep="/"))
source(paste(.CRY_HOME,"lattice_dot.R",sep="/"))
source(paste(.CRY_HOME,"symmetry_dot.R",sep="/"))
source(paste(.CRY_HOME,"atom_dot.R",sep="/"))
source(paste(.CRY_HOME,"reflections_dot.R",sep="/"))
source(paste(.CRY_HOME,"reflections_aux.R",sep="/"))

# Load all ready S4 modules
source(paste(.CRY_HOME,"angle_S4.R",sep="/"))
source(paste(.CRY_HOME,"cell_S4.R",sep="/"))
source(paste(.CRY_HOME,"lattice_S4.R",sep="/"))
source(paste(.CRY_HOME,"symmetry_S4.R",sep="/"))
source(paste(.CRY_HOME,"atom_S4.R",sep="/"))
source(paste(.CRY_HOME,"reflections_S4.R",sep="/"))
