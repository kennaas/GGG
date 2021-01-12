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

group = "Inner"
loter_run = 1
# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
cores = ifelse(Sys.info()[['sysname']] == "Windows", 
               1, detectCores())

##################### Load data ###################################

# Local ancestry data for the group in question
load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/W/W.RData"))
# Haplotype data
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
V1 = BGData(geno = initFileBackedMatrix(
  numInds, numSNPs,  
  folderOut = paste0("Data/loter/Run ", loter_run, "/V1", group),
  outputType = "double"), pheno = data.frame(ID = inds))

V2 = BGData(geno = initFileBackedMatrix(
  numInds, numSNPs,  
  folderOut = paste0("Data/loter/Run ", loter_run, "/V2", group),
  outputType = "double"), pheno = data.frame(ID = inds))

# Center and scale. Should also be parallelized
for (ind in 1:numInds) {
  V1@geno[ind, ] = 
    get(paste0("A", group))@geno[ind, SNPs1] * (W@geno[ind, SNPs1] - p)
  V2@geno[ind, ] = 
    get(paste0("A", group))@geno[ind, SNPs2] * (W@geno[ind, SNPs2] - p)
  print(paste0("V row ", ind, " / ", numInds))
}

save(V1, 
     file = paste0("Data/loter/Run ", loter_run, "/V1", group, "/V1.RData"))

save(V2, 
     file = paste0("Data/loter/Run ", loter_run, "/V2", group, "/V2.RData"))
