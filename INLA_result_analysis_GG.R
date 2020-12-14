######################## Libraries ################################

if (!require(INLA)) 
  install.packages(
    "INLA", 
    repos=c(getOption("repos"), 
            INLA="https://inla.r-inla-download.org/R/stable"), 
    dep=TRUE)
library(INLA)

######################## Load results #############################

# load("Runs/mass/inla_result_mass_generalized_yang_alpha_-1 (different prior!).RData")
# model1 = model
# load("Runs/mass/inla_result_mass_generalized_yang_alpha_-1.RData")
# model2 = model

load("Runs/mass/inla_result_hetGG_rio3_VR_mass.RData")
modelrio = model
# load("Runs/wing/inla_result_hetGG_rio1__full_wing.RData")
# modelriofull = model
load("Runs/mass/inla_result_hetGG_pedigree_mass.RData")
modelped = model
# load("Runs/wing/inla_result_hetGG_pedigree__full_wing.RData")
# modelpedfull = model

# Load rio Gs
# load("Runs/GRMs/GRM_rio_Inner.RData")
# innerGRM = G
# load("Runs/GRMs/GRM_rio_Outer.RData")
# outerGRM = G
# load("Runs/GRMs/GRM_rio_Other.RData")
# otherGRM = G

#Load partial relatedness matrices
# load("Runs/morphData_pedigree_het_version.RData")

######################## INLA result analysis #####################

# inlaPostVariances = function(marginal){
#   sigmaMarg = inla.tmarginal(function(x) 1 / x, marginal)
#   summaryStats =
#     inla.zmarginal(sigmaMarg, silent = TRUE)[c(1, 5, 3, 7)]
#   summaryStats = round(as.numeric(summaryStats), digits = 2)
#   names(summaryStats) = c("mean", "mode", "2.5%", "97.5%")
#   return(summaryStats)
# }
# 
# variances = do.call(
#   "rbind", lapply(modelrio$marginals.hyperpar, inlaPostVariances))
# rownames(variances) =
#   gsub("Precision", "Variance", rownames(variances))
# modelrio$variances = variances


modelrio$summary.fixed
# modelriofull$summary.fixed
modelped$summary.fixed
# modelpedfull$summary.fixed
modelrio$variances
# modelriofull$variances
modelped$variances
# modelpedfull$variances

# # Legarra scaling
# innerGRM_reduced = innerGRM[colnames(innerGRM) %in% morphData$ringnr,
#                             colnames(innerGRM) %in% morphData$ringnr]
# 
# outerGRM_reduced = outerGRM[colnames(outerGRM) %in% morphData$ringnr,
#                             colnames(outerGRM) %in% morphData$ringnr]
# 
# otherGRM_reduced = otherGRM[colnames(otherGRM) %in% morphData$ringnr,
#                             colnames(otherGRM) %in% morphData$ringnr]
# 
# innerGRM_reduced = innerGRM_reduced[diag(innerGRM_reduced) != 0,
#                                     diag(innerGRM_reduced) != 0]
# outerGRM_reduced = outerGRM_reduced[diag(outerGRM_reduced) != 0,
#                                     diag(outerGRM_reduced) != 0]
# otherGRM_reduced = otherGRM_reduced[diag(otherGRM_reduced) != 0,
#                                     diag(otherGRM_reduced) != 0]
# 
# k_rio_inner = mean(diag(innerGRM_reduced)) - mean(innerGRM_reduced)
# k_rio_outer = mean(diag(outerGRM_reduced)) - mean(outerGRM_reduced)
# k_rio_other = mean(diag(otherGRM_reduced)) - mean(otherGRM_reduced)
# 
# k_rio = rbind(k_rio_inner, k_rio_outer, k_rio_other)
# 
# A_inner_reduced = A_inner[colnames(A_inner) %in% morphData$id,
#                           colnames(A_inner) %in% morphData$id]
# 
# A_outer_reduced = A_outer[colnames(A_outer) %in% morphData$id,
#                           colnames(A_outer) %in% morphData$id]
# 
# A_other_reduced = A_other[colnames(A_other) %in% morphData$id,
#                           colnames(A_other) %in% morphData$id]
# 
# A_inner_reduced = A_inner_reduced[diag(A_inner_reduced) != 0,
#                                   diag(A_inner_reduced) != 0]
# A_outer_reduced = A_outer_reduced[diag(A_outer_reduced) != 0,
#                                   diag(A_outer_reduced) != 0]
# A_other_reduced = A_other_reduced[diag(A_other_reduced) != 0,
#                                   diag(A_other_reduced) != 0]
# 
# k_ped_inner = mean(diag(A_inner_reduced)) - mean(A_inner_reduced)
# k_ped_outer = mean(diag(A_outer_reduced)) - mean(A_outer_reduced)
# k_ped_other = mean(diag(A_other_reduced)) - mean(A_other_reduced)
# k_ped = rbind(k_ped_inner, k_ped_outer, k_ped_other)
# 
# 
# 
# apply(modelrio$variances[4:6,], 2, function(x) x * k_rio)
# apply(modelped$variances[4:6,], 2, function(x) x * k_ped)
