######################## Libraries ################################
if (!require(Matrix)) install.packages("Matrix")
library(Matrix)
if (!require(INLA)) 
  install.packages(
    "INLA", 
    repos=c(getOption("repos"), 
            INLA="https://inla.r-inla-download.org/R/stable"), 
    dep=TRUE)
library(INLA)
if (!require(optparse)) install.packages("optparse")
library(optparse)

######################## Options ##################################

# Arguments: Response, GRMVariant, alpha, usePedigree
option_list = list(
  make_option(c("-r", "--response"), type = "character",
              default = "wing", help = "animal model response",
              metavar = "character"),
  make_option(c("-g", "--GRMVariant"), type = "character",
              default = "GCTA", help = "which definition of G",
              metavar = "character"),
  make_option(c("-a", "--alpha"), type = "integer", default = -1,
              help = "alpha parameter in generalized yang",
              metavar = "character"),
  make_option(c("-p", "--usePedigree"), type = "logical", 
              default = FALSE, 
              help = "use pedigree or genome-based matrix", 
              metavar = "character"),
  make_option(c("-i", "--iterate"), type = "logical", 
              default = TRUE, 
              help = "rerun inla model? (improves stability)", 
              metavar = "character"),
  make_option(c("-s", "--shrink"), type = "logical", 
              default = TRUE, 
              help = "rerun inla model? (improves stability)", 
              metavar = "character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Possibilities: wing, mass, more? tarsus?
response = opt$response
# BGData_default, allele_sharing, vanRaden, GCTA
GRMVariant = opt$GRMVariant
# Yang has alpha = -1
alpha = opt$alpha
# Use genomic or pedigree relatedness
usePed = opt$usePedigree
# Rerun inla model or not
iterate = opt$iterate
# Whether to only keep relevant relatednesses
shrink = opt$shrink


######################## Load Data ################################

# Load morphological data
load(file = "Runs/morphData_pedigree_version.RData")
#load(file = "Runs/morphData.RData")

# Load pedigree data
load(file = "Runs/pedigree.RData")

switch (GRMVariant,
        BGData_default = {
          load("Runs/GRMs/GRM_BGData_default.RData")
        },
        allele_sharing = {
          load("Runs/GRMs/GRM_allele_sharing.RData")
        },
        vanRaden = {
          load("Runs/GRMs/GRM_vanRaden.RData") 
        },
        vanRaden_impute = {
          load("Runs/GRMs/GRM_vanRaden_imputed&phased.RData")
        },
        GCTA = {
          GRMVariant = paste0(GRMVariant, "_alpha_", alpha)
          load(paste0("Runs/GRMs/GRM_GCTA_alpha_", alpha, ".RData"))
        }
        
)

######################## Prepare data #############################

# All pedigreed morphdata individuals have been SNPed
# (not the extra non-ped ones, which is clear from using non-ped
# version of morphData) 
# End up with N = 1984 individuls
N = length(unique(morphData$ringnr))
#(N == sum(colnames(GRM) %in% morphData$ringnr))
# No duplicates also
#any(duplicated(colnames(GRM)))

morphData$ID2 = morphData$ID1 = morphData$ID

if (usePed) {
  if (shrink) {
    # Keep only relatednesses for individuals with data
    A = A[colnames(A) %in% morphData$ringnr,
          colnames(A) %in% morphData$ringnr]
    dimnames(A)[[1]] = dimnames(A)[[2]] =
      morphData[match(colnames(A), morphData$ringnr), "ID"]
  } else {
    pedigreeMap$ID2 = 1:nrow(pedigreeMap)
    numRel = dim(A)[1]
    morphData$ID2 =
      pedigreeMap[match(morphData$ringnr, pedigreeMap$ringnr),
                  "ID2"]
  }
  Ainv = solve(A)
  dimnames(Ainv) = dimnames(A)
} else {
  if (shrink) {
    # Keep only relatednesses for individuals with data
    GRM = GRM[colnames(GRM) %in% morphData$ringnr,
              colnames(GRM) %in% morphData$ringnr]
    dimnames(GRM)[[1]] = dimnames(GRM)[[2]] =
      morphData[match(colnames(GRM), morphData$ringnr), "ID"]
  } else {
    pedigreeMap$ID2 = 1:nrow(pedigreeMap)
    numRel = dim(GRM)[1]
    morphData$ID2 =
      pedigreeMap[match(morphData$ringnr, pedigreeMap$ringnr),
                  "ID2"]
  }
}

###################### Make Pos def ###############################

# The function getG() usually gives a pos.def result. So the
# tricks are not needed here? But it would be nice to see what
# trick it uses.
if (usePed) {
  stopifnot(all(eigen(A, only.values = TRUE)$values > 0))
} else {
  stopifnot(all(eigen(GRM, only.values = TRUE)$values > 0))
} 

##################### INLA Animal model ###########################
if (shrink) {
  val = "as.numeric(colnames(Cmatrix))"
} else {
  val = "1:numRel"
}

# Formula
s = paste0(response, " ~ sex + FGRM + month + age + outer + other +
  f(hatchYear, model = \"iid\", 
    hyper = list(prec = list(initial = log(1),
                             prior = \"pc.prec\", 
                             param = c(1, 0.05)))) + 
  f(ID1, model = \"iid\", 
    hyper = list(prec = list(initial = log(1),
                             prior = \"pc.prec\",
                             param = c(2, 0.05)))) + 
  f(ID2, values = ", val, ", 
    model = \"generic0\", Cmatrix = Cmatrix, constr = TRUE,
    hyper = list(prec = list(initial = log(0.69),
                             prior = \"pc.prec\",
                             param = c(3, 0.05)))
    )")
formula = as.formula(s)

# Find inverse of G (or A)
if (usePed) {
  Cmatrix = round(Ainv, 5)
  Cmatrix[Cmatrix == 0] = 0
} else {
  Cmatrix = solve(GRM)
  if (!isSymmetric(Cmatrix)){
    Cmatrix = forceSymmetric(Cmatrix)
  }
}

model = inla(formula = formula, family = "gaussian",
                   data = morphData, verbose = TRUE,
                   control.family = list(
                     hyper = list(theta = list(initial = log(1),
                                               prior = "pc.prec",
                                               param=c(1, 0.05)))))
if (iterate) {
  model = inla.rerun(model)
  model = inla.rerun(model)
}
  

inlaPostVariances = function(marginal){
  sigmaMarg = inla.tmarginal(function(x) 1 / x, marginal)
  summaryStats =
    inla.zmarginal(sigmaMarg, silent = TRUE)[c(1, 5, 3, 7)]
  summaryStats = round(as.numeric(summaryStats), digits = 2)
  names(summaryStats) = c("mean", "mode", "2.5%", "97.5%")
  return(summaryStats)
}

variances = do.call(
  "rbind", lapply(model$marginals.hyperpar, inlaPostVariances))
rownames(variances) =
  gsub("Precision", "Variance", rownames(variances))
model$variances = variances


path = paste0("Runs/inla_result_", response)

if (usePed) {
  path = paste0(path, "_pedigree")
} else {
  path = paste0(path, "_", GRMVariant)
}
if (!shrink)
  path = paste0(path, "_full")

save(model, file = paste0(path, ".RData"))
