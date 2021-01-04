library(BGData)
library(INLA)

# Pick response and local ancestry run
response = "wing"
loter_run = 3

################### Load data #####################################

# load GRMs
load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Inner.RData"))
innerGRM = G_VR
load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Outer.RData"))
outerGRM = G_VR
load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Other.RData"))
otherGRM = G_VR
rm(G_VR, G_GCTA, pi, delta)

# Load phenotypic data
load("Runs/morphData_pedigree_version.RData")

# Load individuals in each ref or admixed populations
load(paste0("Data/Founders", loter_run, "/founders.RData"))
innerInds = na.omit(unique(innerInds))
outerInds = na.omit(unique(outerInds))
otherInds = na.omit(unique(otherInds))
admixedInds = admixedInds[admixedInds %in% morphData$ringnr]

# Load INLA results
load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, ".RData"))
modelRio = model
load(paste0("Runs/", response, "/inla_result_hetGG_pedigree_", response, ".RData"))
modelPed = model
rm(model)
#Load partial relatedness matrices
load("Runs/A_hetped.RData")

######################### Legarra scaling #######################

innerBasePop = c(innerInds, admixedInds)
outerBasePop = c(outerInds, admixedInds)
otherBasePop = c(otherInds, admixedInds)

innerGRM = innerGRM[colnames(innerGRM) %in% innerBasePop,
                    colnames(innerGRM) %in% innerBasePop]
outerGRM = outerGRM[colnames(outerGRM) %in% outerBasePop,
                    colnames(outerGRM) %in% outerBasePop]
otherGRM = otherGRM[colnames(otherGRM) %in% otherBasePop,
                    colnames(otherGRM) %in% otherBasePop]

 
kRioInner = mean(diag(innerGRM)) - mean(innerGRM)
kRioOuter = mean(diag(outerGRM)) - mean(outerGRM)
kRioOther = mean(diag(otherGRM)) - mean(otherGRM)
 
kRio = rbind(kRioInner, kRioOuter, kRioOther)
 
A_inner = as.matrix(A_inner[colnames(A_inner) %in% innerBasePop,
                            colnames(A_inner) %in% innerBasePop])

A_outer = as.matrix(A_outer[colnames(A_outer) %in% outerBasePop,
                            colnames(A_outer) %in% outerBasePop])

A_other = as.matrix(A_other[colnames(A_other) %in% otherBasePop,
                            colnames(A_other) %in% otherBasePop])

kPedInner = mean(diag(A_inner)) - mean(A_inner)
kPedOuter = mean(diag(A_outer)) - mean(A_outer)
kPedOther = mean(diag(A_other)) - mean(A_other)
kPed = rbind(kPedInner, kPedOuter, kPedOther)

getLegarraPosterior = function(model, precision, k) {
  return(inla.zmarginal(inla.tmarginal(
    function(x) k / x,
    model$marginals.hyperpar[[precision]])))
}

innerRioLegarrad = getLegarraPosterior(modelRio, 5, kRioInner)
innerPedLegarrad = getLegarraPosterior(modelPed, 5, kPedInner)

outerRioLegarrad = getLegarraPosterior(modelRio, 6, kRioOuter)
outerPedLegarrad = getLegarraPosterior(modelPed, 6, kPedOuter)

otherRioLegarrad = getLegarraPosterior(modelRio, 7, kRioOther)
otherPedLegarrad = getLegarraPosterior(modelPed, 7, kPedOther)

save(innerRioLegarrad, innerPedLegarrad,
     outerRioLegarrad, outerPedLegarrad,
     otherRioLegarrad, otherPedLegarrad,
     file = paste0("Runs/legarra/", response, loter_run, ".RData"))