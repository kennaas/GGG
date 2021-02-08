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

loter_run = 1
# The parallel implementation in BGData does not work on Windows.
# So if ran on Windows, use only one core (slow), otherwise use all.
numCores = ifelse(Sys.info()[['sysname']] == "Windows", 1, detectCores())

##################### Load data ###################################

# Must have already run W_setup.R and loter_result_conversion.R

# Haplotype data
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))

# Local ancestry data for the group in question
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/AInner/AInner.RData"))
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/AOuter/AOuter.RData"))
load.BGData(file = paste0("Data/loter/Run ", loter_run, "/AOther/AOther.RData"))

#################### Useful vectors ###############################

numInds = dim(W@geno)[1] / 2
numSNPs = dim(W@geno)[2]
inds = as.character(W@pheno$ind[1:numInds])

#################### Finding group-specific allele frequency ######

# See definition of p-estimator
# Temporary matrix to efficiently compute numerator of p-estimator
WAInner = BGData(
  pheno = W@pheno,
  geno = initFileBackedMatrix(2 * numInds, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/WAInner"),
                              outputType = "boolean"))
WAOuter = BGData(
  pheno = W@pheno,
  geno = initFileBackedMatrix(2 * numInds, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/WAOuter"),
                              outputType = "boolean"))
WAOther = BGData(
  pheno = W@pheno,
  geno = initFileBackedMatrix(2 * numInds, numSNPs,
                              folderOut = paste0("Data/loter/Run ", loter_run, "/WAOther"),
                              outputType = "boolean"))

WA_func = function(row) {
  WAInner@geno[row, ] = W@geno[row, ] & AInner@geno[row, ]
  WAOuter@geno[row, ] = W@geno[row, ] & AOuter@geno[row, ]
  WAOther@geno[row, ] = W@geno[row, ] & AOther@geno[row, ]
  return(NULL)
}
trash = mclapply(1:(2 * numInds), WA_func, mc.cores = numCores)
rm(trash)

# Number of individuals with local ancestry in the group for each position in genome
pDenomInner = chunkedApply(AInner@geno, 2, sum, nCores = numCores, verbose = TRUE)
pDenomOuter = chunkedApply(AOuter@geno, 2, sum, nCores = numCores, verbose = TRUE)
pDenomOther = chunkedApply(AOther@geno, 2, sum, nCores = numCores, verbose = TRUE)

# Number of alternate alleles with local ancestry in the group
# for each position in the genome
pNumerInner = chunkedApply(WAInner@geno, 2, sum, nCores = numCores, verbose = TRUE)
pNumerOuter = chunkedApply(WAOuter@geno, 2, sum, nCores = numCores, verbose = TRUE)
pNumerOther = chunkedApply(WAOther@geno, 2, sum, nCores = numCores, verbose = TRUE)

# Group-specific allele frequency vector
pInner = pNumerInner / pDenomInner
pOuter = pNumerOuter / pDenomOuter
pOther = pNumerOther / pDenomOther

save(WAInner, file = paste0("Data/loter/Run ", loter_run, "/WAInner/WAInner.RData"))
save(WAOuter, file = paste0("Data/loter/Run ", loter_run, "/WAOuter/WAOuter.RData"))
save(WAOther, file = paste0("Data/loter/Run ", loter_run, "/WAOther/WAOther.RData"))

# Delete temporary matrices
# unlink(paste0("Data/loter/Run ", loter_run, "/WAInner/geno_1.bin"), recursive = TRUE)
# unlink(paste0("Data/loter/Run ", loter_run, "/WAOuter/geno_1.bin"), recursive = TRUE)
# unlink(paste0("Data/loter/Run ", loter_run, "/WAOther/geno_1.bin"), recursive = TRUE)
# unlink(paste0("Data/loter/Run ", loter_run, "/WAInner"), recursive = TRUE)
# unlink(paste0("Data/loter/Run ", loter_run, "/WAOuter"), recursive = TRUE)
# unlink(paste0("Data/loter/Run ", loter_run, "/WAOther"), recursive = TRUE)

save(pInner, pOuter, pOther, 
     file = paste0("Runs/AlleleFreqs/alleleFreqs", loter_run, ".RData"))
