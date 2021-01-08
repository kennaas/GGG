library(Matrix)
library(data.table)
library(BGData)
source("My_R_code/file_backed_mat.R")

loter_run = 1

load.BGData(file = paste0("Data/loter/Run ", loter_run, 
                          "/W/W.RData"))

cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)
numInds = dim(W@geno)[1]
numSNPs = dim(W@geno)[2]
inds = dimnames(W@geno)[[1]]
SNPs = dimnames(W@geno)[[2]]

SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

########################### Finding F #############################

f_numerator = chunkedApply(W, 2, sum,
                           nCores = cores, verbose = TRUE)
f_numerator =
  unname(tapply(f_numerator,
                (seq_along(f_numerator) - 1) %/% 2, sum))

f = f_numerator / (2 * numInds)

###################### Find G #####################################

G_numerator = getG(W, center = rep(f, each = 2), scale = FALSE,
         scaleG = FALSE, nCores = cores, verbose = TRUE)

G_denominator = sum(rep(p, each = 2) * (1 - rep(p, each = 2)))

GRM = G_numerator / G_denominator
dimnames(GRM)[[1]] = dimnames(GRM)[[2]] = inds

save(GRM, file = paste0("Runs/GRMs/GRM_vanRaden_imputed&phased",
                        loter_run , ".RData"))
