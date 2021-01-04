library(BGData)
library(INLA)

loter_run = 3
cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)

load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Inner.RData"))
innerGRM = G_VR
load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Outer.RData"))
outerGRM = G_VR
load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Other.RData"))
otherGRM = G_VR
rm(G_VR, G_GCTA, pi, delta)

response = "wing"

load("Runs/morphData_pedigree_version.RData")

load(paste0("Data/Founders", loter_run, "/founders.RData"))

innerInds = na.omit(unique(innerInds))
outerInds = na.omit(unique(outerInds))
otherInds = na.omit(unique(otherInds))
admixedInds = admixedInds[admixedInds %in% morphData$ringnr]

load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, ".RData"))
modelRio = model
load(paste0("Runs/", response, "/inla_result_hetGG_pedigree_", response, ".RData"))
modelPed = model
rm(model)
######################### Legarra scaling #######################

#Load partial relatedness matrices
load("Runs/A_hetped.RData")

innerBasePop = c(innerInds)#, admixedInds)
outerBasePop = c(outerInds)#, admixedInds)
otherBasePop = c(otherInds)#, admixedInds)

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

######################## One time: find gamma matrices ############

# load.BGData(file = paste0("Data/loter/Run ", loter_run, 
#                           "/AInner/AInner.RData"))
# 
# load.BGData(file = paste0("Data/loter/Run ", loter_run, 
#                           "/AOuter/AOuter.RData"))
# 
# load.BGData(file = paste0("Data/loter/Run ", loter_run, 
#                           "/AOther/AOther.RData"))

# thetaInner = getG(AInner@geno, 
#              center = FALSE, scale = FALSE, scaleG = FALSE,
#              nCores = cores, verbose = TRUE) 
# thetaInner = thetaInner / dim(geno(AInner))[2]
# gammaInner = innerGRM / thetaInner
# save(gammaInner, file = "Runs/gammaInner.Rdata")
# 
# thetaOuter = getG(AOuter@geno, 
#                   center = FALSE, scale = FALSE, scaleG = FALSE,
#                   nCores = cores, verbose = TRUE) 
# thetaOuter = thetaOuter / dim(geno(AOuter))[2]
# gammaOuter = outerGRM / thetaOuter
# save(gammaOuter, file = "Runs/gammaOuter.Rdata")
# 
# thetaOther = getG(AOther@geno, 
#                   center = FALSE, scale = FALSE, scaleG = FALSE,
#                   nCores = cores, verbose = TRUE) 
# thetaOther = thetaOther / dim(geno(AOther))[2]
# gammaOther = otherGRM / thetaOther
# save(gammaOther, file = "Runs/gammaMOther.Rdata")