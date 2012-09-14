readLog <- function(filename, messages = TRUE) {
  
  ################################################################################
  ## a function for reading a text file containing CCP4 loggraph tables         ##
  ##                                                                            ##
  ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
  ################################################################################  
  
  # see http://www.ccp4.ac.uk/html/loggraphformat.html
  
  #read whole file as text lines into a character vector
  logText <- readLines(filename)
  
  #identify the lines that start tables
  tableLines <- grep(pattern = "$TABLE", x = logText, fixed = TRUE)
  if(length(tableLines) == 0) stop(paste("No tables found in", filename))
  
  if(messages) cat(length(tableLines),"tables have been found in", filename, "\n")
  
  #set up output list
  tableList <- vector(mode = "list", length = length(tableLines))
  tableNames <- sub(".*\\$TABLE[[:space:]]*:[[:space:]]*", "", logText[tableLines])
  tableNames <- sub("[[:space:]]*:$","", tableNames)
  names(tableList) <- tableNames
  
  for (i in seq_along(tableList)){
  
    fileLine <- tableLines[i] + 1
    #$GRAPHS or $SCATTER lines follow. Could plot these, but for the moment ignore them.
    
    #find next line containing '$$', defining the start of column headings
    while(!grepl("\\$\\$", logText[fileLine])) fileLine <- fileLine + 1
    
    #now concatenate lines from this point until all column headings, 'any character' lines
    #and numbers in this table are identified
    string <- sub("[^\\$]*\\$\\$", "", logText[fileLine])
    while(!grepl("\\$\\$.*\\$\\$.*\\$\\$", string)){
      fileLine <- fileLine + 1
      string <- paste(string, logText[fileLine])
    }
    
    #now split the string and remove leading white space from each element
    string <- strsplit(string, "\\$\\$")[[1]]
    string <- sub("^[[:space:]]*", "", string)
    
    columnHeadings <- strsplit(string[[1]], "[[:space:]]+")[[1]]
    anyChar <- strsplit(string[[2]], "[[:space:]]+")[[1]]
    numbers <- as.numeric(strsplit(string[[3]], "[[:space:]]+")[[1]])
    
    table <- matrix(data = numbers, byrow = TRUE, ncol = length(columnHeadings))
    colnames(table) <- columnHeadings
    tableList[[i]] <- as.data.frame(table)
  }
if(messages) summary(tableList)
return(tableList)
}