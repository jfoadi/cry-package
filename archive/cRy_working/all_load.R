# File to source all R modules and load all auxiliary files 
# needed to run cRy

# Root directory (full path to src_cRy). Change this according to your settings
CRY_HOME <- "/Users/james/Dropbox/CRY/cRy_working"   # Use this to load modules as in e.g. "source(paste(CRY_HOME,"symmetry.R"))"

# Auxiliary files. You shouldn't need editing this, as files will be located in relation to CRY_HOME
SYMMETRY_file <- paste(CRY_HOME,"syminfo.lib",sep="/")

# Load all ready modules
source(paste(CRY_HOME,"lattice.R",sep="/"))
source(paste(CRY_HOME,"symmetry.R",sep="/"))
source(paste(CRY_HOME,"mtz.R",sep="/"))
source(paste(CRY_HOME,"map.R",sep="/"))
source(paste(CRY_HOME,"loggraph.R",sep="/"))
source(paste(CRY_HOME,"pdb.R",sep="/"))
source(paste(CRY_HOME,"miscellaneous.R",sep="/"))
