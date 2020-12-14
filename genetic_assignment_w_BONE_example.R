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

################## Example ########################################

Classes = c(rep("character", 4), rep("integer", 2000))
MixtureData = read.table("Bone/Data/MixtureData.txt", 
                         colClasses = Classes, header = TRUE)
ReferenceData = read.table("Bone/Data/ReferenceData.txt",
                           colClasses = Classes, header = TRUE)

str(MixtureData[, 1:6])
str(ReferenceData[, 1:6])

MixtureY = impute(MixtureData, genotyping = "Axiom")
#MixtureY = impute(MixtureData, genotyping = "Axiom")
MixtureY = MixtureY$Y
MixtureY[1:5,1:4]
ReferenceY = impute(ReferenceData, genotyping = "Axiom")
ReferenceY = ReferenceY$Y

setdiff(rownames(MixtureY),rownames(ReferenceY))

Markers = intersect(rownames(MixtureY),rownames(ReferenceY))
MixtureY = MixtureY[Markers,]
ReferenceY = ReferenceY[Markers,]

SampleSizeBaselinePop = ncol(ReferenceY)
Y = cbind(ReferenceY,MixtureY)

lambda = seq(0.4,0.02,length.out=40)

MBapprox = LASSOSolPath(Y, lambda, intercept = TRUE,
                        SampleSizeBaselinePop, Baseline = TRUE)

SampleNames=colnames(Y)

NetworkResultsSolpath = SolPathInference(MBapprox, ReferenceData,
                                         SampleNames = SampleNames,
                                         SampleSizeBaselinePop,
                                         alpha = 0.05)
NetworkResultsWTA = WTAInference(MBapprox, ReferenceData,
                                 SampleNames = SampleNames, 
                                 SampleSizeBaselinePop)


NetworkResultsWTA$ProbofOrigin[1:5,]
NetworkResultsWTA$MixtureProp
round(NetworkResultsSolpath$ProbofOrigin[1:5,], 4)
NetworkResultsSolpath$MixtureProp
