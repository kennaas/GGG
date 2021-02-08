##################### Libraries ##################################
if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(parallel)) install.packages("parallel")
library(parallel)
source("My_R_code/file_backed_mat.R")

####################### Setup #####################################


group = "Inner"
loter_run = 1

# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
numCores = ifelse(Sys.info()[['sysname']] == "Windows", 1, detectCores())

load.BGData(file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, "/A", group, "/A", group, ".RData"))

#################### Useful vectors ###############################

numInds = dim(W@geno)[1] / 2
numSNPs = dim(W@geno)[2]
inds = as.character(W@pheno$ind[1:numInds])

###################################################################


print("getG() for h=1")
theta_11 = getG(get(paste0("A", group))@geno, i = 1:numInds,
                center = FALSE, scale = FALSE, scaleG = FALSE,
                nCores = numCores, verbose = TRUE)

print("getG() for h=2")
theta_22 = getG(get(paste0("A", group))@geno, i = (numInds + 1):(2 * numInds),
                center = FALSE, scale = FALSE, scaleG = FALSE,
                nCores = numCores, verbose = TRUE)

print("getG() for cross-comparison")
theta_12 = getG(get(paste0("A", group))@geno, i = 1:numInds, i2 = (numInds + 1):(2 * numInds),
                center = FALSE, scale = FALSE, scaleG = FALSE,
                nCores = numCores, verbose = TRUE)

theta_21 = t(theta_12)

# Sum and average over all matrix products to get theta
theta = (theta_11 + theta_12 + theta_21 + theta_22) / (4 * numSNPs)
dimnames(theta) = list(inds, inds)

save(theta, file = paste0("Data/loter/Run ", loter_run, "/theta", group, ".RData"))
