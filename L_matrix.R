##################### Libraries ##################################
if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(parallel)) install.packages("parallel")
library(parallel)
if (!require(remotes)) install.packages("remotes")
library(remotes)
if (!require(fs)) install.packages("fs")
library(fs)
source("My_R_code/file_backed_mat.R")

##################### Settings ####################################

group = "Other"
loter_run = 1
# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
cores = ifelse(Sys.info()[['sysname']] == "Windows", 
               1, detectCores())

##################### Load data ###################################

# Haplotype data
load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/W/W.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/A", group, "/A", group, ".RData"))

load(file = paste0("Runs/AlleleFreqs/alleleFreqs", loter_run, ".RData"))
p = get(paste0("p", group))

#################### Useful vectors ###############################

numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2] / 2
inds = W@pheno$ID
# Indices
SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

###################### Find V #####################################

# Temporary haplotype matrix
L1 = BGData(geno = initFileBackedMatrix(
  numInds, numSNPs,  
  folderOut = paste0("Data/loter/Run ", loter_run, "/L1", group),
  outputType = "double"), pheno = data.frame(ID = inds))

L2 = BGData(geno = initFileBackedMatrix(
  numInds, numSNPs,  
  folderOut = paste0("Data/loter/Run ", loter_run, "/L2", group),
  outputType = "double"), pheno = data.frame(ID = inds))

p_scaling = sqrt(p * (1 - p))

head(SNPs1)
SNPs1 = as.integer(SNPs1)
SNPs2 = as.integer(SNPs2)

# a = seq(from = 1, by = 35, to = numInds)
# for (i in seq_along(a)[-1]) {
#   L1@geno[(a[i-1]):(a[i]), ] = chunkedApply(get(paste0("A", group))@geno[(a[i-1]):(a[i]), ],
#                                             nCores = cores, MARGIN = 1,
#                                             j = SNPs1, verbose = TRUE,
#                                 function(row) row * p_scaling, chunkSize = 3)
#   
#   L2@geno[(a[i-1]):(a[i]), ] = chunkedApply(get(paste0("A", group))@geno[(a[i-1]):(a[i]), ],
#                                             nCores = cores, MARGIN = 1,
#                                             j = SNPs2, verbose = TRUE,
#                                 function(row) row * p_scaling, chunkSize = 3)
#   print(paste0("L row ", a[i-1], " to ", a[i]))
# }

for (row in 1:numInds) {
  L1@geno[row, ] = get(paste0("A", group))@geno[row, SNPs1] * p_scaling
  L2@geno[row, ] = get(paste0("A", group))@geno[row, SNPs2] * p_scaling
  print(paste0("L row ", row, " / ", numInds))
}

save(L1, 
     file = paste0("Data/loter/Run ", loter_run, "/L1", group, "/L1.RData"))

save(L2, 
     file = paste0("Data/loter/Run ", loter_run, "/L2", group, "/L2.RData"))
