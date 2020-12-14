if (!require(devtools)) install.packages("devtools")
library(devtools)
if (!require(roxygen2)) install.packages("roxygen2")
library(roxygen2)
if (!require(BONE)) install("BONE")
library(BONE)
if (!require(glmnet)) install.packages("glmnet")
library(glmnet)
if (!require(data.table)) install.packages("data.table")
library(data.table)


load("Runs/morphData_pedigree_version.RData")
load("Runs/pedigree.RData")
# Reference populations consist of all purebreeds here:
load("Data/referenceIDs1_purebreed.RData")

fullPed = fread(file = "Data/Data_from_henrik/Helgeland_01_2018.ped")

#cutFullPed = fullPed[, c(2, 7:100)]

admixed = (fullPed$V2 %in% admixedRef)

admixedPed = fullPed[admixed, ]
referencePed = fullPed[!admixed, ]

#rm(cutFullPed)
rm(fullPed)

admixedPed = cbind(data.frame(
  sample_type = rep("mixture", sum(admixed)),
  repunit = rep(100, sum(admixed)),
  collection = rep(100)), admixedPed)
colnames(admixedPed)[4] = "indiv"

referencePed = cbind(data.frame(
  sample_type = rep("reference", sum(!admixed)),
  repunit = rep(-1, sum(!admixed)),
  collection = rep(-1, sum(!admixed))), referencePed)
colnames(referencePed)[4] = "indiv"
referencePed[referencePed$indiv %in% innerRef, 2:3] = c(0, 0)
referencePed[referencePed$indiv %in% outerRef, 2:3] = c(1, 1)
referencePed[referencePed$indiv %in% otherRef, 2:3] = c(2, 2)

# admixedCutPed = ifelse(admixedCutPed == "A", 1, 
#                        ifelse(admixedCutPed == "C", 2,
#                               ifelse(admixedCutPed == "G", 3,
#                                      ifelse(admixedCutPed == "T", 4,
#                                             admixedCutPed))))

MixtureY = impute(admixedPed, genotyping = "Axiom")
#MixtureY = impute(admixedCutPed, genotyping = "Axiom")
#MixtureY = MixtureY$Y
MixtureY[1:5,1:4]
ReferenceY = impute(referencePed, genotyping = "Axiom")
#ReferenceY = ReferenceY$Y
