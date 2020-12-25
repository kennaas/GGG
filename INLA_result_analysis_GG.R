######################## Libraries ################################

if (!require(INLA)) 
  install.packages(
    "INLA", 
    repos=c(getOption("repos"), 
            INLA="https://inla.r-inla-download.org/R/stable"), 
    dep=TRUE)
library(INLA)

######################## Prior sensitivity analysis ###############
response = "wing"

load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, ".RData"))
modelmoderate = model
load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, "_strictprior.RData"))
modelstrict = model
load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, "_kindprior.RData"))
modelkind = model

load(paste0("Runs/", response, "/inla_result_hetGG_pedigree_", response, ".RData"))
modelpedmoderate = model
load(paste0("Runs/", response, "/inla_result_hetGG_pedigree_", response, "_strictprior.RData"))
modelpedstrict = model
load(paste0("Runs/", response, "/inla_result_hetGG_pedigree_", response, "_kindprior.RData"))
modelpedkind = model


modelmoderate$summary.fixed
modelkind$summary.fixed
modelstrict$summary.fixed
modelmoderate$variance
modelkind$variance
modelstrict$variance

modelpedmoderate$summary.fixed
modelpedkind$summary.fixed
modelpedstrict$summary.fixed
modelpedmoderate$variance
modelpedkind$variance
modelpedstrict$variance

rm(list = ls())

######################### Run times ###############################
response = "wing"

load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, ".RData"))
summary(model01)
summary(model02)
summary(model)

riomodel = model

load(paste0("Runs/", response, "/inla_result_hetGG_pedigree_", response, ".RData"))
summary(model01)
summary(model02)
summary(model)

pedmodel = model

rm(model)
rm(model01)
rm(model02)

############# Using all genotyped (3) vs. only phenotyped (4) #####

response = "tarsus"

load(paste0("Runs/", response, "/inla_result_hetGG_rio3_VR_", response, ".RData"))

model3 = model

load(paste0("Runs/", response, "/inla_result_hetGG_rio4_VR_", response, ".RData"))

model4 = model

model3$summary.fixed - model4$summary.fixed

model3$variances - model4$variances

# all responses: barely effects results -> safe to use smaller data (faster) data set?
