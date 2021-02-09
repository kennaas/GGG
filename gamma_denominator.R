##################### Libraries ##################################
if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(parallel)) install.packages("parallel")
library(parallel)
source("My_R_code/file_backed_mat.R")

##################### Settings ####################################

group = "Inner"
loter_run = 1
# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, detectCores())

##################### Load data ###################################

# Local ancestry data for the group in question
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/A", group, "/A", group, ".RData"))

load(file = paste0("Runs/AlleleFreqs/alleleFreqs", loter_run, ".RData"))
p = get(paste0("p", group))
pScaling = 1 / sqrt(p * (1 - p))
rm(p, pInner, pOuter, pOther)

#################### Useful vectors ###############################

numInds = dim(get(paste0("A", group))@geno)[1] / 2
numSNPs = dim(get(paste0("A", group))@geno)[2]
inds = as.character(get(paste0("A", group))@pheno$ind[1:numInds])

####################### Find gamma denominator ####################

gammaDenomVR11 = getG(get(paste0("A", group))@geno, i = 1:numInds,
                      center = FALSE, scale = pScaling, scaleG = FALSE,
                      nCores = cores, verbose = TRUE)

gammaDenomVR22 = getG(get(paste0("A", group))@geno, i = (numInds + 1):(2 * numInds),
                      center = FALSE, scale = pScaling, scaleG = FALSE,
                      nCores = cores, verbose = TRUE)

gammaDenomVR12 = getG(get(paste0("A", group))@geno, i = 1:numInds, i2 = (numInds + 1):(2 * numInds),
                      center = FALSE, scale = pScaling, scaleG = FALSE,
                      nCores = cores, verbose = TRUE)

gammaDenomVR21 = t(gammaDenomVR12)

# Sum the four matrix products
gammaDen = gammaDenomVR11 + gammaDenomVR12 + gammaDenomVR21 + gammaDenomVR22

print("zero?")
any(gammaDen == 0)
gammaDen = ifelse(gammaDen == 0, 1e-12, gammaDen)

save(gammaDen, file = paste0("Data/loter/Run ", loter_run, "/gammaDenominator", group, ".RData"))
