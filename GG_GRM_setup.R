##################### Libraries ##################################
if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(BGData)) install.packages("BGData")
library(BGData)
source("My_R_code/file_backed_mat.R")

##################### Settings ####################################

group = "Other"
loter_run = 3
cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)

##################### Load data ###################################

# Local ancestry data for the group in question
load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/W/W.RData"))
# Haplotype data
load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/A", group, "/A", group, ".RData"))

#################### Useful vectors ###############################

numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2] / 2
inds = W@pheno$ID

#################### Find pi ######################################

pi = chunkedApply(get(paste0("A", group))@geno, MARGIN = 1, mean,
                  nCores = cores, verbose = TRUE, chunkSize = 250L)

#################### Finding group-specific allele frequency ######

# See definition of p-estimator

# Temporary matrix to efficiently compute numerator of pi-estimator
WA = initFileBackedMatrix(numInds, 2 * numSNPs,
                          folderOut = paste0("Data/loter/Run ", 
                                             loter_run, "/WA",
                                             group),
                          outputType = "boolean")

for (ind in 1:numInds) {
  WA[ind, ] = W@geno[ind, ] * get(paste0("A", group))@geno[ind, ]
  print(paste0("Subject ", ind, " / ", numInds))
}

# Number of individuals with local ancestry in the group
# for each position in genome
p_denominator = chunkedApply(get(paste0("A", group))@geno, 2, sum, 
                             nCores = cores, verbose = TRUE)

# Merge the two entries on the same locus
p_denominator =
  unname(tapply(p_denominator,
                (seq_along(p_denominator) - 1) %/% 2, sum))

# Number of alternate alleles with local ancestry in the group
# for each position in the genome
p_numerator = chunkedApply(WA, 2, sum,
                           nCores = cores, verbose = TRUE)
# Merge the two entries on the same locus
p_numerator =
  unname(tapply(p_numerator,
                (seq_along(p_numerator) - 1) %/% 2, sum))

# Group-specific allele frequency vector
p = p_numerator / p_denominator

# Delete temporary matrix
rm(WA)
unlink(paste0("Data/loter/Run ", loter_run, "/WA", group, 
              "/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WA", group), 
       recursive = TRUE)
rm(p_numerator)
rm(p_denominator)

###################### Find V #####################################

# Temporary haplotype matrix
V = BGData(geno = initFileBackedMatrix(
  numInds, 2 * numSNPs,  
  folderOut = paste0("Data/loter/Run ", loter_run, "/V", group),
  outputType = "double"), pheno = data.frame(ID = inds))

# Center and scale 
for (ind in 1:numInds) {
  V@geno[ind, ] = get(paste0("A", group))@geno[ind, ] * (W@geno[ind, ] - rep(p, each = 2))
  print(paste0("Subject ", ind, " / ", numInds))
}

##################### Find alpha ##################################

# VR-version
s = 1 / sqrt((rep(p, each = 2) * (1 - rep(p, each = 2))))

gamma_denominator_VR = getG(get(paste0("A", group))@geno, 
                            center = FALSE, scale = s, 
                            scaleG = FALSE, nCores = cores, 
                            verbose = TRUE)

gamma_numerator_VR = getG(V@geno,  center = FALSE, scale = FALSE,
                          scaleG = FALSE, nCores = cores, 
                          verbose = TRUE)

gamma_denominator_VR = ifelse(gamma_denominator_VR == 0,
                              1e-12, gamma_denominator_VR)
gamma_VR = gamma_numerator_VR / gamma_denominator_VR


# GCTA - version 1 (local ancestries in num and mult. by theta)

gamma_GCTA1 = getG(V@geno, center = FALSE, scale = 1 / s, 
                   scaleG = TRUE, nCores = cores, verbose = TRUE)

rm(V)
unlink(paste0("Data/loter/Run ", loter_run, "/V", group, 
              "/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/V", group), 
       recursive = TRUE)

# GCTA - version 2 (only mult. by theta)

gamma_GCTA2 = getG(W@geno, center = rep(p, each = 2), 
                   scale = 1 / s, scaleG = TRUE, nCores = cores,
                   verbose = TRUE)

###################### Find theta #################################

theta = getG(get(paste0("A", group))@geno, 
             center = FALSE, scale = FALSE, scaleG = FALSE,
             nCores = cores, verbose = TRUE) 
theta = theta / (2 * numSNPs)

###################### Find Delta #################################

delta = theta - tcrossprod(pi)

##################### Find G ######################################

G_VR = gamma_VR * theta
G_GCTA1 = gamma_GCTA1 * theta
G_GCTA2 = gamma_GCTA2 * theta

dimnames(G_VR)[[1]] = dimnames(G_VR)[[2]] =
  dimnames(G_GCTA1)[[1]] = dimnames(G_GCTA1)[[2]] = 
  dimnames(G_GCTA2)[[1]] = dimnames(G_GCTA2)[[2]] =
  dimnames(delta)[[1]] = dimnames(delta)[[2]] =
  inds

# Check if GRMs are invertible
print("G_VR invertible:")
all(eigen(G_VR + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)
print("G_GCTA1 invertible:")
all(eigen(G_GCTA1 + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)
print("G_GCTA2 invertible:")
all(eigen(G_GCTA2 + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)


# Save GRMs
save(G_VR, G_GCTA1, G_GCTA2, delta, pi,
     file = paste0("Runs/GRMs/GRM_rio", loter_run, 
                   "_", group, ".RData"))

