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

packageVersion("BGData")
old_lib = ifelse(Sys.info()[['sysname']] == "Windows",
                 path_home_r("R/win-library/old-versions/"),
                 "/usr/local/lib/R/old-versions/")
# One time:
# install_version("BGData", version = "1.0", lib = old_lib)
packageVersion("BGData")
source("My_R_code/file_backed_mat.R")

##################### Settings ####################################

group = "Other"
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

####################### Find gamma denominator ####################

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/L1", group, "/L1.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/L2", group, "/L2.RData"))

gamma_denominator_VR_11 = getG(L1@geno,
                               center = FALSE, scale = FALSE,
                               scaleG = FALSE, nCores = cores,
                               verbose = TRUE)

save(gamma_denominator_VR_11, file = "safe1.RData")

gamma_denominator_VR_22 = getG(L2@geno,
                               center = FALSE, scale = FALSE,
                               scaleG = FALSE, nCores = cores,
                               verbose = TRUE)

save(gamma_denominator_VR_22, file = "safe2.RData")

detach("package:BGData")
library("BGData", lib.loc = old_lib)
packageVersion("BGData")

# Split up computation in two steps to alleviate memory issues
print("b1:")
b1 = tcrossprod_parallel(x = L1@geno, 
                         y = L2@geno[1:(numInds %/% 2), ],
                         nCores = cores)
print("b2:")
b2 = tcrossprod_parallel(x = L1@geno, 
                         y = L2@geno[(numInds %/% 2 + 1):numInds, ],
                         nCores = cores)
print("cbind:")
gamma_denominator_VR_12 = cbind(b1, b2)
print("transpose:")
gamma_denominator_VR_21 = t(gamma_denominator_VR_12)

# Sum the four matrix products
gamma_den = 
  gamma_denominator_VR_11 + gamma_denominator_VR_12 +
  gamma_denominator_VR_21 + gamma_denominator_VR_22

gamma_den = ifelse(gamma_den == 0, 1e-12, gamma_den)

save(gamma_den, 
     file = paste0("Data/loter/Run ", loter_run, "/gammaDenominator", group, ".RData"))
