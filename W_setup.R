suppressMessages({
if (!require(R.utils)) install.packages("R.utils")
library(R.utils)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(foreach)) install.packages("foreach")
library(foreach)
})
  
source("My_R_code/file_backed_mat.R")

  
loter_run = 1

get_pop = function(group) {
  groupFile = file(paste0("Data/loter/Run ", loter_run, "/Input/", group, "_phased.vcf.gz"),
                   open = "r")
  groupLine = scan(groupFile, nlines = 1, skip = 8,
                   what = character(), quiet = TRUE)
  close(groupFile)
  return(groupLine[-(1:9)])
}

# Population is grouped by inner, outer, other, admixed
innerInds = get_pop("inner")
outerInds = get_pop("outer")
otherInds = get_pop("other")
admixedInds = get_pop("admixed")
pop = c(innerInds, outerInds, otherInds, admixedInds)
# Number of inds. in each (sub)population
numInnerInds = length(innerInds)
numOuterInds = length(outerInds)
numOtherInds = length(otherInds)
numAdmixedInds = length(admixedInds)
popSize = length(pop)

SNPs = fread(file = paste0("Data/loter/Run ", loter_run, "/Input/admixed_phased.vcf.gz"),
             select = 3)$ID
numSNPs = length(SNPs)

# Create file-backed matrix to hold haplotypes.
# First popSize rows correspond to DNA strand 1 and next popSize correspond to DNA strand 2
W = BGData(pheno = data.frame(ind = rep(pop, 2), h = rep(c(1, 2), each = popSize)),
           geno = initFileBackedMatrix(nrows = 2 * popSize, ncols = numSNPs, outputType = "boolean",
                                       folderOut = paste0("Data/loter/Run ", loter_run, "/W")))


indHaplo = function(ind, firstInd, groupVCF) {
  indVCF = ind + 1 - firstInd
  W@geno[ind, ] = sapply(groupVCF[indVCF, ], function(x) as.numeric(grepl("^1", x)))
}

# Fill matrix with haplotype data from the VCF format
fillW = function(group, firstInd, lastInd) {
  print("Loading haplotypes...")
  groupVCF = fread(file = paste0("Data/loter/Run ", loter_run, "/Input/", group, "_phased.vcf.gz"),
                   drop = 1:9)
  print("Parsing haplotypes...")
  h1 = lapply(groupVCF, function(x) grepl("^1", x))
  h2 = lapply(groupVCF, function(x) grepl("1$", x))

  rm(groupVCF)
  
  print("Saving haplotypes...")
  for (ind in firstInd:lastInd) {
    iVCF = ind + 1 - firstInd
    W@geno[ind, ] = h1[[iVCF]]
    W@geno[ind + popSize, ] = h2[[iVCF]]
    # print(paste0("Individual ", ind))
  }
  rm(h1, h2)
  
  # mclapply(firstInd:lastInd, indHaplo, firstInd, groupVCF, mc.cores = numCores)
  # rm(groupVCF)
  print("Done.")
}

print("Parsing inner")
fillW("inner", 1, numInnerInds)

print("Parsing outer")
fillW("outer", numInnerInds + 1, numInnerInds + numOuterInds)

print("Parsing other")
fillW("other", numInnerInds + numOuterInds + 1, popSize - numAdmixedInds)

print("Parsing admixed")
fillW("admixed", numInnerInds + numOuterInds + numOtherInds + 1, popSize)

rownames(W@map) = colnames(W@geno) = SNPs
rownames(W@pheno) = rownames(W@pheno) = apply(W@pheno, 1, paste, collapse = "_")

# Save
save(W, file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))
