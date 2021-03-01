library(Matrix)
library(data.table)
library(BGData)
library(parallel)
library(tictoc)
source("My_R_code/file_backed_mat.R")

loter_run = 1

load.BGData(file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))

numCores = ifelse(Sys.info()[['sysname']] == "Windows", 1, detectCores())
numInds = dim(W@geno)[1] / 2
numSNPs = dim(W@geno)[2]
inds = as.character(W@pheno$indS)[1:numInds]

passes = round(numSNPs / 5000 / numCores)
getGChunks = ceiling(numSNPs / (passes * numCores))

p = chunkedApply(W@geno, 2, sum, 
                 nCores = numCores, verbose = TRUE, chunkSize = getGChunks) / (2 * numInds)

G_11 = getG(W@geno, i = 1:numInds,
            center = p, scale = FALSE, scaleG = FALSE,
            nCores = numCores, verbose = TRUE, chunkSize = getGChunks)

G_22 = getG(W@geno, i = (numInds + 1):(2 * numInds),
            center = p, scale = FALSE, scaleG = FALSE,
            nCores = numCores, verbose = TRUE, chunkSize = getGChunks)

G_12 = getG(W@geno, i = 1:numInds, i2 = (numInds + 1):(2 * numInds),
            center = p, scale = FALSE, scaleG = FALSE,
            nCores = numCores, verbose = TRUE, chunkSize = getGChunks)

G_21 = t(G_12)

GRM = (G_11 + G_22 + G_12 + G_21) / (2 * sum(p * (1 - p)))
dimnames(GRM) = list(inds, inds)

save(GRM, file = paste0("Runs/GRMs/GRM_vanRaden_imputed&phased", loter_run , ".RData"))
