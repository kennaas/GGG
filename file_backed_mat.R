library(Matrix)
library(BGData)
library(ff)

# Help functions to create custom file backed matrices

# Function to initialize file backed system
ffNodeInitializer = function(nodeIndex, nrow, ncol, vmode, folderOut, ...) {
  filename <- paste0("geno_", nodeIndex, ".bin")
  node <- ff(dim = c(nrow, ncol), vmode = vmode, filename = paste0(folderOut, "/", filename), ...)
  # Change ff path to a relative one
  physical(node)[["filename"]] <- filename
  return(node)
}

# Initializes a file back matrix of a specific size (nrow x ncol),
# and data type. Data types:
# boolean = 1 bit without NA.
# logical = 2 bit (with NA)
# byte (integers up to 2^7)
# double
# etc.
initFileBackedMatrix = function(nrows, ncols, folderOut, outputType) {
  dir.create(folderOut)

  chunkSize = min(nrows, floor(.Machine[["integer.max"]] / ncols / 1.2))
  nNodes = ceiling(nrows / chunkSize)
  
  matrix = LinkedMatrix(nrow = nrows, ncol = ncols, nNodes = nNodes,
                        linkedBy = "rows", nodeInitializer = ffNodeInitializer,
                        vmode = outputType, folderOut = folderOut, dimorder = 2:1)
  return(matrix)
}
