if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(BGData)) install.packages("BGData")
library(BGData)
source("My_R_code/file_backed_mat.R")

group = "Other"
loter_run = 4

load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/W/W.RData"))
load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/A", group, "/A", group, ".RData"))

cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)

numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2] / 2
inds = W@pheno$ID

######################## Find pi ##################################

pi = chunkedApply(get(paste0("A", group))@geno, MARGIN = 1, mean,
                  nCores = cores, verbose = TRUE, chunkSize = 250L)

######################## Finding f ################################

WA = initFileBackedMatrix(numInds, 2 * numSNPs,
                          folderOut = paste0("Data/loter/Run ", 
                                             loter_run, "/WA",
                                             group),
                          outputType = "boolean")

for (ind in 1:numInds) {
  WA[ind, ] = W@geno[ind, ] * get(paste0("A", group))@geno[ind, ]
  print(paste0("Subject ", ind, " / ", numInds))
}

f_denominator = chunkedApply(get(paste0("A", group))@geno, 2, sum, 
                             nCores = cores, verbose = TRUE)
f_denominator =
  unname(tapply(f_denominator,
                (seq_along(f_denominator) - 1) %/% 2, sum))

f_numerator = chunkedApply(WA, 2, sum,
                           nCores = cores, verbose = TRUE)
f_numerator =
  unname(tapply(f_numerator,
                (seq_along(f_numerator) - 1) %/% 2, sum))

f = f_numerator / f_denominator

rm(WA)
unlink(paste0("Data/loter/Run ", loter_run, "/WA", group, 
              "/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/WA", group), 
       recursive = TRUE)
rm(f_numerator)
rm(f_denominator)

###################### Find V #####################################

V = BGData(geno = initFileBackedMatrix(
  numInds, 2 * numSNPs,  
  folderOut = paste0("Data/loter/Run ", loter_run, "/V", group),
  outputType = "double"), pheno = data.frame(ID = inds))

for (ind in 1:numInds) {
  V@geno[ind, ] = get(paste0("A", group))@geno[ind, ] * (W@geno[ind, ] - rep(f, each = 2))
  print(paste0("Subject ", ind, " / ", numInds))
}

##################### Find alpha ##################################

# VR-version
s = 1 / sqrt((rep(f, each = 2) * (1 - rep(f, each = 2))))

alpha_denominator_VR = getG(get(paste0("A", group))@geno, 
                            center = FALSE, scale = s, 
                            scaleG = FALSE, nCores = cores, 
                            verbose = TRUE)

alpha_num_VR = getG(V@geno, 
                    center = FALSE, scale = FALSE, scaleG = FALSE,
                    nCores = cores, verbose = TRUE)

alpha_denominator_VR = ifelse(alpha_denominator_VR == 0,
                              1e-12, alpha_denominator_VR)
alpha_VR = alpha_num_VR / alpha_denominator_VR


# GCTA - version 1 (local ancestries in num and mult. by theta)

alpha_GCTA1 = getG(V@geno, center = FALSE, scale = 1 / s, 
                   scaleG = TRUE, nCores = cores, verbose = TRUE)

rm(V)
unlink(paste0("Data/loter/Run ", loter_run, "/V", group, 
              "/geno_1.bin"), 
       recursive = TRUE)
unlink(paste0("Data/loter/Run ", loter_run, "/V", group), 
       recursive = TRUE)

# GCTA - version 2 (only mult. by theta)

alpha_GCTA2 = getG(W@geno, 
                   center = rep(f, each = 2), scale = 1 / s, 
                   scaleG = TRUE, nCores = cores, verbose = TRUE)

###################### Find theta #################################

theta = getG(get(paste0("A", group))@geno, 
             center = FALSE, scale = FALSE, scaleG = FALSE,
             nCores = cores, verbose = TRUE) 
theta = theta / (2 * numSNPs)

###################### Find Delta #################################

delta = theta - tcrossprod(pi)

##################### Find G ######################################

G_VR = alpha_VR * theta
G_GCTA1 = alpha_GCTA1 * theta
G_GCTA2 = alpha_GCTA2 * theta

dimnames(G_VR)[[1]] = dimnames(G_VR)[[2]] =
  dimnames(G_GCTA1)[[1]] = dimnames(G_GCTA1)[[2]] = 
  dimnames(G_GCTA2)[[1]] = dimnames(G_GCTA2)[[2]] =
  dimnames(delta)[[1]] = dimnames(delta)[[2]] =
  inds

print("G_VR invertible:")
all(eigen(G_VR + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)
print("G_GCTA1 invertible:")
all(eigen(G_GCTA1 + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)
print("G_GCTA2 invertible:")
all(eigen(G_GCTA2 + diag(0.01, numInds, numInds),
          only.values = TRUE)$values > 0)

save(G_VR, G_GCTA1, G_GCTA2, delta, pi,
     file = paste0("Runs/GRMs/GRM_rio", loter_run, 
                   "_", group, ".RData"))

