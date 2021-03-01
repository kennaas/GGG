##################### Libraries ##################################
if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(parallel)) install.packages("parallel")
library(parallel)
if (!require(tictoc)) install.packages("tictoc")
library(tictoc)
source("My_R_code/file_backed_mat.R")

##################### Settings ####################################

group = "Other"
loter_run = 1
# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
numCores = ifelse(Sys.info()[['sysname']] == "Windows", 1, detectCores())

##################### Load data ###################################

# Haplotype data
# load.BGData(file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))

# Local ancestry data for the group in question
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/A", group, "/A", group, ".RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, "/WA", group, "/WA", group, ".RData"))

load(file = paste0("Runs/AlleleFreqs/alleleFreqs", loter_run, ".RData"))
p = get(paste0("p", group))
rm(pInner, pOuter, pOther)

numInds = dim(get(paste0("A", group))@geno)[1] / 2
numSNPs = dim(get(paste0("A", group))@geno)[2]

##################### Create temporary V matrix ##################

V = BGData(
  pheno = get(paste0("A", group))@pheno,
  geno = initFileBackedMatrix(2 * numInds, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/V", group),
                              outputType = "double"))

V_func = function(row) {
  V@geno[row, ] = get(paste0("WA", group))@geno[row, ] - get(paste0("A", group))@geno[row, ] * p
  return(NULL)
}

trash = mclapply(1:(2 * numInds), V_func, mc.cores = numCores)
rm(trash)
##################### Find gamma numerator ########################

passes = round(numSNPs / 5000 / numCores)
getGChunks = ceiling(numSNPs / (passes * numCores))
# Compute four matrix products: all permuations of h=1 and h=2
gammaNumerVR11 = getG(V@geno, i = 1:numInds,
                      center = FALSE, scale = FALSE, scaleG = FALSE,
                      nCores = numCores, verbose = TRUE, chunkSize = getGChunks)

gammaNumerVR22 = getG(V@geno, i = (numInds + 1):(2 * numInds),
                      center = FALSE, scale = FALSE, scaleG = FALSE,
                      nCores = numCores, verbose = TRUE, chunkSize = getGChunks)

gammaNumerVR12 = getG(V@geno, i = 1:numInds, i2 = (numInds + 1):(2 * numInds),
                      center = FALSE, scale = FALSE, scaleG = FALSE,
                      nCores = numCores, verbose = TRUE, chunkSize = getGChunks)

gammaNumerVR21 = t(gammaNumerVR12)

# Sum the four matrix products
gammaNumer = gammaNumerVR11 + gammaNumerVR12 + gammaNumerVR21 + gammaNumerVR22

save(gammaNumer, file = paste0("Data/loter/Run ", loter_run, "/gammaNumerator", group, ".RData"))
