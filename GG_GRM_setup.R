library(Matrix)
library(data.table)
library(BGData)
source("My_R_code/file_backed_mat.R")

group = "Other"

load.BGData(file = "Data/loter/Run 3/W/W.RData")
load.BGData(file = paste0("Data/loter/Run 3/A", group, "/A", group, ".RData"))

cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)

numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2] / 2
inds = W@pheno$ID

V = BGData(geno = initFileBackedMatrix(
  numInds, 2 * numSNPs,  
  folderOut = paste0("Data/loter/Run 3/V", group),
  outputType = "double"), pheno = data.frame(ID = inds))

######################## Find pi ##################################

pi = chunkedApply(get(paste0("A", group))@geno, 1, mean, nCores = cores)

######################## Finding f ################################

WA = initFileBackedMatrix(numInds, 2 * numSNPs,
                          folderOut = paste0("Data/loter/Run 3/WA", group),
                          outputType = "boolean")

for (ind in 1:numInds) {
  row = W@geno[ind, ] * get(paste0("A", group))@geno[ind, ]
  WA[ind, ] = row
  print(paste0("Subject ", ind, " / ", numInds, " ",
               all(WA[ind, ] == row)))
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
unlink(paste0("Data/loter/Run 3/WA", group))
rm(f_numerator)
rm(f_denominator)

###################### Find V #####################################

for (ind in 1:numInds) {
  row = get(paste0("A", group))@geno[ind, ] * (W@geno[ind, ] - rep(f, each = 2))
  V@geno[ind, ] = row
  print(paste0("Subject ", ind, " / ", numInds, " ", all(V@geno[ind, ] == row)))
}

##################### Find alpha ##################################

# VR-version
s = 1 / sqrt((rep(f, each = 2) * (1 - rep(f, each = 2))))

alpha_denominator_VR = getG(get(paste0("A", group))@geno, center = FALSE, scale = s, 
                         scaleG = FALSE, nCores = cores, 
                         verbose = TRUE)

alpha_num_VR = getG(V@geno, center = FALSE, scale = FALSE, scaleG = FALSE,
         nCores = cores, verbose = TRUE)

alpha_denominator_VR = ifelse(alpha_denominator_VR == 0,
                           1e-12, alpha_denominator_VR)
alpha_VR = alpha_num_VR / alpha_denominator_VR


# GCTA - version 1 (local ancestries in num and mult. by theta)

alpha_GCTA1 = getG(V@geno, center = FALSE, scale = 1 / s, 
                  scaleG = TRUE, nCores = cores, verbose = TRUE)

rm(V)
unlink(paste0("Data/loter/Run 3/V", group), recursive = TRUE)

# GCTA - version 2 (only mult. by theta)

alpha_GCTA2 = getG(W@geno, center = rep(f, each = 2), scale = 1 / s, 
                   scaleG = TRUE, nCores = cores, verbose = TRUE)

###################### Find theta #################################

theta = getG(get(paste0("A", group))@geno, center = FALSE, scale = FALSE, scaleG = FALSE,
             nCores = cores, verbose = TRUE) 
theta = theta / (2 * numSNPs)

###################### Find Delta #################################

delta = theta - tcrossprod(pi)

##################### Find G ######################################

G_VR = alpha_VR * theta
G_GCTA1 = alpha_GCTA1 * theta
G_GCTA2 = alpha_GCTA2 * theta

solve(G_VR + diag(0.01, numInds, numInds))
solve(G_GCTA1 + diag(0.01, numInds, numInds))
solve(G_GCTA2 + diag(0.01, numInds, numInds))

dimnames(G_VR)[[1]] = dimnames(G_VR)[[2]] = inds
dimnames(G_GCTA1)[[1]] = dimnames(G_GCTA1)[[2]] = inds
dimnames(G_GCTA2)[[1]] = dimnames(G_GCTA2)[[2]] = inds

save(G_VR, G_GCTA1, G_GCTA2, delta, pi,
     file = paste0("Runs/GRMs/GRM_rio3_", group, ".RData"))

