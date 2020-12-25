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


load("Runs/wing/inla_result_wing_vanRaden.RData")
model1 = model
load("Runs/wing/inla_result_wing_vanRaden_impute.RData")
model2 = model
load("Runs/wing/inla_result_wing_pedigree.RData")
model3 = model

# Load rio Gs
# load("Runs/GRMs/GRM_rio_Inner.RData")
# innerGRM = G
# load("Runs/GRMs/GRM_rio_Outer.RData")
# outerGRM = G
# load("Runs/GRMs/GRM_rio_Other.RData")
# otherGRM = G

######################## INLA result analysis #####################

model1$summary.fixed
model2$summary.fixed
model3$summary.fixed
model1$variances
model2$variances
model3$variances
