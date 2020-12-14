library(Matrix)
library(data.table)
library(BGData)
source("My_R_code/file_backed_mat.R")


load.BGData(file = "Data/BGData_W1_beagle1/BGData.RData")
W1 = BGData
load.BGData(file = "Data/BGData_W2_beagle1/BGData.RData")
W2 = BGData
cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)
numInds = dim(W1@geno)[1]
numSNPs = dim(W1@geno)[2]
inds = dimnames(W1@geno)[[1]]
SNPs = dimnames(W1@geno)[[2]]

W = initFileBackedMatrix(numInds, 2 * numSNPs,
                         folderOut = "W")

SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

######################## Splice together haplotype matrices #######


for (ind in 1:numInds) {
  W[ind, SNPs1] = W1@geno[ind, ]
  W[ind, SNPs2] = W2@geno[ind, ]
  print(paste0("Subject ", ind, " / ", numInds))
}
rm(W1)
rm(W2)

########################### Finding F #############################

f_numerator = chunkedApply(W, 2, sum,
                           nCores = cores, verbose = TRUE)
f_numerator =
  unname(tapply(f_numerator,
                (seq_along(f_numerator) - 1) %/% 2, sum))

f = f_numerator / (2 * numInds)

###################### Find G #####################################

G = getG(W, center = rep(f, each = 2), scale = FALSE,
         scaleG = FALSE, nCores = cores, verbose = TRUE)

G_denominator = 2 * sum(f * (1 - f))

GRM = G / G_denominator
dimnames(GRM)[[1]] = inds
dimnames(GRM)[[2]] = inds

save(GRM, file = "Runs/GRMs/GRM_vanRaden_imputed&phased.RData")
