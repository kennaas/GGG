library(Matrix)
library(data.table)
library(BGData)
library(parallel)
source("My_R_code/file_backed_mat.R")

loter_run = 1

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/W/W.RData"))

cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, detectCores())
numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2]
inds = W@pheno$ID
SNPs = colnames(W@geno)

SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

########################### Finding p #############################

p_numerator = chunkedApply(W@geno, 2, sum,
                           nCores = cores, verbose = TRUE)
p_numerator =
  unname(tapply(p_numerator,
                (seq_along(p_numerator) - 1) %/% 2, sum))

p = p_numerator / (2 * numInds)

###################### Find G #####################################

G_numerator = getG(W@geno, center = rep(p, each = 2), scale = FALSE,
         scaleG = FALSE, nCores = cores, verbose = TRUE)

# G_denominator = sum(p * (1 - p))
G_denominator = 2 * sum(p * (1 - p)) # mutliply by 2 because haplotype-level?

GRM = G_numerator / G_denominator
dimnames(GRM)[[1]] = dimnames(GRM)[[2]] = inds

save(GRM, file = paste0("Runs/GRMs/GRM_vanRaden_imputed&phased",
                        loter_run , ".RData"))
