library(Matrix)
library(BGData)
library(ff)


ffNodeInitializer = function(nodeIndex, nrow, ncol, vmode, folderOut, ...) {
  filename <- paste0("geno_", nodeIndex, ".bin")
  node <- ff(dim = c(nrow, ncol), vmode = vmode, filename = paste0(folderOut, "/", filename), ...)
  # Change ff path to a relative one
  physical(node)[["filename"]] <- filename
  return(node)
}

# Default: boolean = 1 bit without NA.
initFileBackedMatrix = function(nrows, ncols, folderOut, outputType) {
  dir.create(folderOut)

  chunkSize = min(nrows, floor(.Machine[["integer.max"]] / ncols / 1.2))
  nNodes = ceiling(nrows / chunkSize)
  
  matrix = LinkedMatrix(nrow = nrows, ncol = ncols, nNodes = nNodes,
                        linkedBy = "rows", nodeInitializer = ffNodeInitializer,
                        vmode = outputType, folderOut = folderOut, dimorder = 2:1)
  return(matrix)
}
