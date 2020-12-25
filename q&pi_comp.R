################## Check group membership proportions #############
# How similar is pi to q?
library(data.table)
library(BGData)
library(knitr)

load("Runs/morphData_pedigree_version.RData")
load("Runs/pedigree.RData")
load.BGData(file = "Data/loter/Run 3/AInner/AInner.RData")
load.BGData(file = "Data/loter/Run 3/AOuter/AOuter.RData")
load.BGData(file = "Data/loter/Run 3/AOther/AOther.RData")

inds = AInner@pheno$ID

z = rep(-1, length(inds))
pi_mat = data.table(pi_outer = z, pi_inner = z, pi_other = z)

num_alleles =  dim(AInner@geno)[2]

for (ind in 1:length(inds)) {
  pi_mat$pi_inner[ind] = mean(AInner@geno[ind, ])
  pi_mat$pi_outer[ind] = mean(AOuter@geno[ind, ])
  pi_mat$pi_other[ind] = mean(AOther@geno[ind, ])
  
  print(paste0("Individual ", as.character(ind), ": \n",
               as.character(pi_mat[ind, ])))
}

pi_mat$ID = inds

head(pedigreeData[pedigreeData$ringnr %in% inds, c(1, 8, 9, 10)])

head(pedigreeData[match(inds, pedigreeData$ringnr), c(1, 8, 9, 10)])
head(pi_mat)

comp = cbind(
  pedigreeData[match(inds, pedigreeData$ringnr), c(1, 8, 9, 10)], 
  pi_mat[, c(2, 1, 3)])

save(comp, file = "Runs/q&pi_comp3.RData")