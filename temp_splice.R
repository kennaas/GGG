library(Matrix)
library(data.table)
library(BGData)
source("My_R_code/file_backed_mat.R")

group = "Inner"

# load.BGData(file = "Data/loter/Run 1/BGData_W1_beagle1/BGData.RData")
# W1 = BGData
# load.BGData(file = "Data/loter/Run 1/BGData_W2_beagle1/BGData.RData")
# W2 = BGData
load.BGData(file = paste0(
  "Data/loter/Run 1/BGData_LA_allele1", group, "/BGData.RData"))
A1 = BGData
load.BGData(file = paste0(
  "Data/loter/Run 1/BGData_LA_allele2", group, "/BGData.RData"))
A2 = BGData
rm(BGData)

cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)

numInds = dim(A1@geno)[1]
numSNPs = dim(A1@geno)[2]
inds = dimnames(A1@geno)[[1]]
SNPs = dimnames(A1@geno)[[2]]

A = BGData(geno = initFileBackedMatrix(numInds, 2 * numSNPs, outputType = "logical",
                         folderOut = paste0("Data/loter/Run 1/A", group)),
           pheno = data.frame(ID = inds))

# W = BGData(geno = initFileBackedMatrix(numInds, 2 * numSNPs, outputType = "logical",
#                          folderOut = "Data/loter/Run 1/W"), 
#            pheno = data.frame(ID = inds))


SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

# for (ind in 1:numInds) {
#   W@geno[ind, SNPs1] = W1@geno[ind, ]
#   W@geno[ind, SNPs2] = W2@geno[ind, ]
#   print(paste0("Subject ", ind, " / ", numInds))
# }
# rm(W1)
# rm(W2)

# colnames(W@geno) = do.call("paste0", c(list(rep(SNPs, each = 2)), "_", list(1:2)))

for (ind in 1:numInds) {
  A@geno[ind, SNPs1] = A1@geno[ind, ]
  A@geno[ind, SNPs2] = A2@geno[ind, ]
  print(paste0("Subject ", ind, " / ", numInds))
}
rm(A1)
rm(A2)

load.BGData(file = "Data/loter/Run 1/W/W.RData")
dimnames(A@geno) = dimnames(W@geno)

# save(W, file = "Data/loter/Run 1/W/W.RData")
save(A, file = paste0("Data/loter/Run 1/A", group, "/A", group, ".RData"))
