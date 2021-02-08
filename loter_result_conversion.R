library(data.table)
library(BGData)
source("My_R_code/file_backed_mat.R")

# Must have already run W_setup.R

loter_run = 1

# Find purebred and admixed populations  
get_pop = function(group) {
  groupFile = file(paste0("Data/loter/Run ", loter_run, "/Input/", group, "_phased.vcf.gz"),
                   open = "r")
  groupLine = scan(groupFile, nlines = 1, skip = 8, what = character(), quiet = TRUE)
  close(groupFile)
  return(groupLine[-(1:9)])
}

innerInds = get_pop("inner")
outerInds = get_pop("outer")
otherInds = get_pop("other")
admixedInds = get_pop("admixed")
numInnerInds = length(innerInds)
numOuterInds = length(outerInds)
numOtherInds = length(otherInds)
numAdmixedInds = length(admixedInds)
pop = c(innerInds, outerInds, otherInds, admixedInds)
popSize = length(pop)

SNPs = fread(file = paste0("Data/loter/Run ", loter_run, "/Input/admixed_phased.vcf.gz"),
             select = 3)$ID
numSNPs = length(SNPs)

# Create BGData objects for each local ancestry matrix.
# First popSize rows correspond to DNA strand 1 and next popSize correspond to DNA strand 2
AOuter = BGData(
  geno = initFileBackedMatrix(2 * popSize, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/AOuter"),
                              outputType = "boolean"),
  pheno = data.frame(ind = rep(pop, 2), h = rep(c(1, 2), each = popSize)))
AInner = BGData(
  geno = initFileBackedMatrix(2 * popSize, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/AInner"),
                              outputType = "boolean"),
  pheno = data.frame(ind = rep(pop, 2), h = rep(c(1, 2), each = popSize)))
AOther = BGData(
  geno = initFileBackedMatrix(2 * popSize, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/AOther"),
                              outputType = "boolean"),
  pheno = data.frame(ind = rep(pop, 2), h = rep(c(1, 2), each = popSize)))

# For inner rows, set inner LA to 1, outer and other to 0
print("Inner individuals...")
for (ind in 1:numInnerInds) {
  AOuter@geno[ind, ] = AOther@geno[ind, ] =
    AOuter@geno[ind + popSize, ] = AOther@geno[ind + popSize, ] = rep(FALSE, numSNPs)
  AInner@geno[ind, ] = AInner@geno[ind + popSize, ] = rep(TRUE, numSNPs)
}

print("Outer individuals...")
# For outer rows, set outer LA to 1, inner and other to 0
for (ind in (numInnerInds + 1):(numInnerInds + numOuterInds)) {
  AInner@geno[ind, ] = AOther@geno[ind, ] =
    AInner@geno[ind + popSize, ] = AOther@geno[ind + popSize, ] = rep(FALSE, numSNPs)
  AOuter@geno[ind, ] = AOuter@geno[ind + popSize, ] = rep(TRUE, numSNPs)
}

print("Other individuals...")
# For other rows, set other LA to 1, inner and outer to 0
for (ind in (numInnerInds + numOuterInds + 1):(popSize - numAdmixedInds)) {
  AInner@geno[ind, ] = AOuter@geno[ind, ] =
    AInner@geno[ind + popSize, ] = AOuter@geno[ind + popSize, ] = rep(FALSE, numSNPs)
  AOther@geno[ind, ] = AOther@geno[ind + popSize, ] = rep(TRUE, numSNPs)
}

# In admixed rows, use loter results to determine values
# Order in loter: Outer == 0, Inner == 1, Other == 2
# (since we happened to place them in this order in loter)
print("Admixed individuals...")
admixedFile = file(paste0("Data/loter/Run ", loter_run, "/loter_result.txt"), 
                   open = "r")
for (ind in (popSize - numAdmixedInds + 1):popSize) {
  line1  = scan(admixedFile, nlines = 1, what = integer(), quiet = TRUE)
  AOuter@geno[ind, ] = line1 == 0
  AInner@geno[ind, ] = line1 == 1
  AOther@geno[ind, ] = line1 == 2
  
  line2 = scan(admixedFile, nlines = 1, what = integer(), quiet = TRUE)
  AOuter@geno[ind + popSize, ] = line2 == 0
  AInner@geno[ind + popSize, ] = line2 == 1
  AOther@geno[ind + popSize, ] = line2 == 2
  # print(row)
}
close(admixedFile)

# Set dimnames
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))
dimnames(AInner@geno) = dimnames(AOuter@geno) = dimnames(AOther@geno) = dimnames(W@geno)
dimnames(AInner@pheno) = dimnames(AOuter@pheno) = dimnames(AOther@pheno) = dimnames(W@pheno)
dimnames(AInner@map) = dimnames(AOuter@map) = dimnames(AOther@map) = dimnames(W@map)

# Save converted data in three local ancestry matrices with 0/1 entries
save(AInner, file = paste0("Data/loter/Run ", loter_run, "/AInner/AInner.RData"))
save(AOuter, file = paste0("Data/loter/Run ", loter_run, "/AOuter/AOuter.RData"))
save(AOther, file = paste0("Data/loter/Run ", loter_run, "/AOther/AOther.RData"))
