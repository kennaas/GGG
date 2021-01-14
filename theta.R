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

group = "Other"
loter_run = 1
# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
cores = ifelse(Sys.info()[['sysname']] == "Windows", 
               1, detectCores())

load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/W/W.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/A", group, "/A", group, ".RData"))

#################### Useful vectors ###############################

numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2] / 2
inds = W@pheno$ID
# Indices
SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

##################################################################

detach("package:BGData")
library("BGData", lib.loc = old_lib)
packageVersion("BGData")

block1 = tcrossprod_parallel(x = get(paste0("A", group))@geno[, SNPs1],
                             y = get(paste0("A", group))@geno[
                               1:(numInds %/% 3), SNPs2],
                             nCores = cores)
block2 = tcrossprod_parallel(x = get(paste0("A", group))@geno[, SNPs1],
                             y = get(paste0("A", group))@geno[
                               (numInds %/% 3 + 1):(2 * (numInds %/% 3)), SNPs2],
                             nCores = cores)
block3 = tcrossprod_parallel(x = get(paste0("A", group))@geno[, SNPs1],
                             y = get(paste0("A", group))@geno[
                               (2 * (numInds %/% 3) + 1):(numInds), SNPs2],
                             nCores = cores)
theta_12 = cbind(block1, block2, block3)

theta_21 = t(theta_12)

save(theta_12, theta_21, file = "safe1.RData")

detach("package:BGData")
library("BGData")
packageVersion("BGData")

theta_11 = getG(get(paste0("A", group))@geno, 
                j = SNPs1,
                center = FALSE, scale = FALSE, scaleG = FALSE,
                nCores = cores, verbose = TRUE)
save(theta_11, file = "safe2.RData")

theta_22 = getG(get(paste0("A", group))@geno, 
                j = SNPs2,
                center = FALSE, scale = FALSE, scaleG = FALSE,
                nCores = cores, verbose = TRUE)

save(theta_22, file = "safe3.RData")

theta = (theta_11 + theta_12 + theta_21 + theta_22) / (4 * numSNPs)

save(theta, 
     file = paste0("Data/loter/Run ", loter_run, "/theta", group,".RData"))
