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

################### Function to flush matrix ######################

flush_ff = function(name) {
  rm(get(name))
  unlink(paste0("Data/loter/Run ", loter_run, "/", name, group, 
                "/geno_1.bin"), 
         recursive = TRUE)
  unlink(paste0("Data/loter/Run ", loter_run, "/", name, group), 
         recursive = TRUE)
}

#################### Find pi ######################################

print("Finding pi")

pi = chunkedApply(get(paste0("A", group))@geno, MARGIN = 1, mean,
                  nCores = cores, verbose = TRUE, chunkSize = 250L)


###################### Find gamma #################################

load(paste0("Data/loter/Run ", loter_run, "/gammaNumerator", group, ".RData"))
load(paste0("Data/loter/Run ", loter_run, "/gammaDenominator", group, ".RData"))

# The 1 / 2 is to give the proper standardization (similar to vanRaden)
gamma_VR = gamma_num / (gamma_den / 2)

# experiment: group-specific version of GCTA (rather than vanRaden)

# GCTA - version 1 (local ancestries in num and mult. by theta)

# gamma_GCTA1 = getG(V@geno, center = FALSE, scale = 1 / s, 
#                    scaleG = TRUE, nCores = cores, verbose = TRUE)

# GCTA - version 2 (only mult. by theta)

# gamma_GCTA2 = getG(W@geno, center = rep(p, each = 2), 
#                    scale = 1 / s, scaleG = TRUE, nCores = cores,
#                    verbose = TRUE)

###################### Find theta #################################

load(paste0("Data/loter/Run ", loter_run, "/theta", group,".RData"))

###################### Find Delta #################################

delta = theta - tcrossprod(pi)

##################### Find G ######################################

G_VR = gamma_VR * theta
# G_GCTA1 = gamma_GCTA1 * theta
# G_GCTA2 = gamma_GCTA2 * theta

dimnames(G_VR)[[1]] = dimnames(G_VR)[[2]] =
  # dimnames(G_GCTA1)[[1]] = dimnames(G_GCTA1)[[2]] = 
  # dimnames(G_GCTA2)[[1]] = dimnames(G_GCTA2)[[2]] =
  dimnames(delta)[[1]] = dimnames(delta)[[2]] =
  inds

# Check if GRMs are invertible
print("G_VR invertible:")
all(eigen(G_VR + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)
# print("G_GCTA1 invertible:")
# all(eigen(G_GCTA1 + diag(0.01, numInds, numInds),
#           only.values = TRUE)$values > 0)
# print("G_GCTA2 invertible:")
# all(eigen(G_GCTA2 + diag(0.01, numInds, numInds),
#           only.values = TRUE)$values > 0)

G_GCTA1 = G_GCTA2 = NULL

# Save GRMs
save(G_VR, G_GCTA1, G_GCTA2, delta, pi,
     file = paste0("Runs/GRMs/GRM_rioFIXED", loter_run, 
                   "_", group, ".RData"))

