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
                          "/AInner/AInner.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/AOuter/AOuter.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/AOther/AOther.RData"))

#################### Useful vectors ###############################

numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2] / 2
inds = W@pheno$ID

#################### Finding group-specific allele frequency ######

# See definition of p-estimator

# Temporary matrix to efficiently compute numerator of p-estimator
WAInner = initFileBackedMatrix(numInds, 2 * numSNPs,
                               folderOut = paste0("Data/loter/Run ",
                                                  loter_run, "/WAInner"),
                               outputType = "boolean")

WAOuter = initFileBackedMatrix(numInds, 2 * numSNPs,
                               folderOut = paste0("Data/loter/Run ",
                                                  loter_run, "/WAOuter"),
                               outputType = "boolean")

WAOther = initFileBackedMatrix(numInds, 2 * numSNPs,
                               folderOut = paste0("Data/loter/Run ",
                                                  loter_run, "/WAOther"),
                               outputType = "boolean")

# Would be nice to parallelize this part too...
# Maybe a future BGData release will have a function chunkedOuter()? (or write self)
for (ind in 1:numInds) {
  WAInner[ind, ] = W@geno[ind, ] * AInner@geno[ind, ]
  WAOuter[ind, ] = W@geno[ind, ] * AOuter@geno[ind, ]
  WAOther[ind, ] = W@geno[ind, ] * AOther@geno[ind, ]
  print(paste0("Subject ", ind, " / ", numInds))
}

# Number of individuals with local ancestry in the group
# for each position in genome
pDenomInner = chunkedApply(AInner@geno, 2, sum,
                           nCores = cores, verbose = TRUE)

pDenomOuter = chunkedApply(AOuter@geno, 2, sum,
                           nCores = cores, verbose = TRUE)

pDenomOther = chunkedApply(AOther@geno, 2, sum,
                           nCores = cores, verbose = TRUE)

# Merge the two entries on the same locus
pDenomInner =
  unname(tapply(pDenomInner,
                (seq_along(pDenomInner) - 1) %/% 2, sum))

pDenomOuter =
  unname(tapply(pDenomOuter,
                (seq_along(pDenomOuter) - 1) %/% 2, sum))

pDenomOther =
  unname(tapply(pDenomOther,
                (seq_along(pDenomOther) - 1) %/% 2, sum))

# Number of alternate alleles with local ancestry in the group
# for each position in the genome
pNumerInner = chunkedApply(WAInner, 2, sum,
                           nCores = cores, verbose = TRUE)

pNumerOuter = chunkedApply(WAOuter, 2, sum,
                           nCores = cores, verbose = TRUE)

pNumerOther = chunkedApply(WAOther, 2, sum,
                           nCores = cores, verbose = TRUE)

# Merge the two entries on the same locus
pNumerInner =
  unname(tapply(pNumerInner,
                (seq_along(pNumerInner) - 1) %/% 2, sum))

pNumerOuter =
  unname(tapply(pNumerOuter,
                (seq_along(pNumerOuter) - 1) %/% 2, sum))

pNumerOther =
  unname(tapply(pNumerOther,
                (seq_along(pNumerOther) - 1) %/% 2, sum))

# Group-specific allele frequency vector
pInner = pNumerInner / pDenomInner
pOuter = pNumerOuter / pDenomOuter
pOther = pNumerOther / pDenomOther

# Delete temporary matrix

unlink(paste0("Data/loter/Run ", loter_run, "/WAInner/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WAOuter/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WAOther/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WAInner"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WAOuter"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WAOther"), 
       recursive = TRUE)

save(pInner, pOuter, pOther, 
     file = paste0("Runs/AlleleFreqs/alleleFreqs", 
                   loter_run, ".RData"))

