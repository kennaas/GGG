################## Check group membership proportions #############

# How similar is pi to q?
library(data.table)
library(BGData)
library(knitr)

###################### Load data ##################################

load("Runs/morphData_pedigree_version.RData")
load("Runs/pedigree.RData")
load("Runs/GRMs/GRM_rio3_Inner.RData")


pi_mat = data.table(ID = rownames(G_VR), pi_inner = pi)
load("Runs/GRMs/GRM_rio3_Outer.RData")
pi_mat$pi_outer = pi
load("Runs/GRMs/GRM_rio3_Other.RData")
pi_mat$pi_other = pi

#head(pedigreeData[pedigreeData$ringnr %in% pi_mat$ID, c(1, 10, 11, 12)])

#head(pedigreeData[match(pi_mat$ID, pedigreeData$ringnr), c(1, 10, 11, 12)])
#head(pi_mat)

# Create data.frame containing group membership proportions based 
# on pedigree and genome.
comp = cbind(pedigreeData[match(pi_mat$ID, pedigreeData$ringnr), c(1, 10:12)], pi_mat[, -1])

save(comp, file = "Runs/q&pi_comp3.RData")