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

# List of options for running from command line
# Response: "wing", "mass" or "tarsus"
# GRMVariant: "VR" or "GCTA". VR is recommended, GCTA is not tested
# UsePedigree: TRUE/FALSE on whether to use pedigree-based model. 
  # FALSE (default) means to use genome-based model.
# Iterate: TRUE/FALSE on whether to rerun the model. 
  # TRUE recommended for better stability
# Shrink: TRUE/FALSE on whether to trim kinship matrices to 
  # only phenotyped individuals. TRUE recommended, as computations
  # are sped up and results are not notably affected.
# Segregation: TRUE/FALSE on whether to attempt to fit segregation
  # variances as random effect. FALSE recommended for convergence.
# h-value: choose custom values for the h-parameter in INLA.
# local_ancestry: which run of loter to apply (3 was used in analysis)
option_list = list(
  make_option(c("-r", "--response"), type = "character",
              default = "wing", help = "animal model response",
              metavar = "character"),
  make_option(c("-g", "--GRMVariant"), type = "character",
              default = "VR", help = "which definition of G",
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
              help = "only keep relevant relatedness entries?", 
              metavar = "character"),
  make_option(c("-z", "--zegregation", type = "logical"),
              default = FALSE,
              help = "include segregation variances?",
              metavar = "character"),
  make_option(c("-h", "--h_value", type = "numeric"),
              default = inla.set.control.inla.default()$h,
              help = "value of h in INLA",
              metavar = "character"),
  make_option(c("-l", "--local_ancestry", type = "numeric"),
              default = 3,
              help = "which run of loter local ancestry inference",
              metavar = "character")
)
opt_parser = OptionParser(option_list = option_list, 
                          add_help_option = FALSE)
opt = parse_args(opt_parser)

# Possibilities: wing, mass, more? tarsus?
response = opt$response

#Which GRM to use
GRMVariant = opt$GRMVariant

# Use genomic or pedigree relatedness
usePed = opt$usePedigree

# Rerun inla model or not
iterate = opt$iterate

# Whether to only keep relevant relatednesses
shrink = opt$shrink

# Whether to include segregation variances
segregation = opt$zegregation

h_val = opt$h_value

loter_run = opt$local_ancestry

######################## Load Data ################################
load(file = "Runs/morphData_pedigree_version.RData")
load(file = "Runs/pedigree.RData")

if (usePed) {
  # Load data
  load(file = "Runs/A_hetped.RData")
} else {
  
  # Choose type of GRM
  load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Inner.RData"))
  if (GRMVariant == "GCTA1") {
    innerGRM = G_GCTA1
  } else if (GRMVariant == "GCTA2") {
    innerGRM = G_GCTA2
  } else if(GRMVariant == "VR") {
    innerGRM = G_VR
  }
  morphData$innerPi = pi[match(morphData$ringnr, colnames(innerGRM))]
  innerDelta = delta
  dimnames(innerDelta) = dimnames(innerGRM)
  
  load(paste0("Runs/GRMs/GRM_rio", loter_run ,"_Outer.RData"))
  if (GRMVariant == "GCTA1") {
    outerGRM = G_GCTA1
  } else if (GRMVariant == "GCTA2") {
    outerGRM = G_GCTA2
  } else if(GRMVariant == "VR") {
    outerGRM = G_VR
  }
  morphData$outerPi = pi[match(morphData$ringnr, colnames(outerGRM))]
  outerDelta = delta
  dimnames(outerDelta) = dimnames(outerGRM)
  
  load(paste0("Runs/GRMs/GRM_rio", loter_run, "_Other.RData"))
  if (GRMVariant == "GCTA1") {
    otherGRM = G_GCTA1
  } else if (GRMVariant == "GCTA2") {
    otherGRM = G_GCTA2
  } else if(GRMVariant == "VR") {
    otherGRM = G_VR
  }
  morphData$otherPi = pi[match(morphData$ringnr, colnames(otherGRM))]
  otherDelta = delta
  dimnames(otherDelta) = dimnames(otherGRM)
  
  rm(G_VR)
  rm(G_GCTA1)
  rm(G_GCTA2)
  rm(pi)
  rm(delta)
}

######################## Prepare data #############################

# All pedigreed morphdata individuals have been SNPed
# (not the extra non-ped ones, which is clear from using non-ped
# version of morphData) 
# End up with N = 1984 individuls
length(unique(morphData$ringnr))

# (N == sum(colnames(otherGRM) %in% morphData$ringnr))
# No duplicates also
#any(duplicated(colnames(GRM)))

morphData$ID4 = morphData$ID3 =
  morphData$ID2 = morphData$ID1 = morphData$ID

if (usePed) {
  if (shrink) {
    # Only kinship for individuals with data
    A_inner = A_inner[colnames(A_inner) %in% morphData$ringnr,
                      colnames(A_inner) %in% morphData$ringnr]
    # Rename dimnames for use in INLA
    dimnames(A_inner)[[1]] = dimnames(A_inner)[[2]] =
      morphData[match(colnames(A_inner), morphData$ringnr), "ID"]
    
    A_outer = A_outer[colnames(A_outer) %in% morphData$ringnr,
                      colnames(A_outer) %in% morphData$ringnr]
    dimnames(A_outer)[[1]] = dimnames(A_outer)[[2]] =
      morphData[match(colnames(A_outer), morphData$ringnr), "ID"]
    
    A_other = A_other[colnames(A_other) %in% morphData$ringnr,
                      colnames(A_other) %in% morphData$ringnr]
    dimnames(A_other)[[1]] = dimnames(A_other)[[2]] =
      morphData[match(colnames(A_other), morphData$ringnr), "ID"]
    
    # IDs used to fit random effects in INLA
    morphData$ID2 = ifelse(morphData$inner > 0,
                               morphData$ID2, NA)
    
    morphData$ID3 = ifelse(morphData$outer > 0,
                               morphData$ID3, NA)
    
    morphData$ID4 = ifelse(morphData$other > 0,
                               morphData$ID4, NA)
  } else {
    pedigreeMap$ID2 = 1:nrow(pedigreeMap)
    
    numRel = dim(A)[1]
    morphData$ID2 = morphData$ID3 = morphData$ID4 =
      pedigreeMap[match(morphData$ringnr, pedigreeMap$ringnr),
                  "ID2"]
  }
} else {
  if (shrink) {
    # Only keep kinship for individuals with data and rename dims
    
    # Inner
    innerGRM = innerGRM[colnames(innerGRM) %in% morphData$ringnr,
                        colnames(innerGRM) %in% morphData$ringnr]
    dimnames(innerGRM)[[1]] = dimnames(innerGRM)[[2]] =
      morphData[match(colnames(innerGRM), morphData$ringnr), "ID"]
    
    # Outer
    outerGRM = outerGRM[colnames(outerGRM) %in% morphData$ringnr,
                        colnames(outerGRM) %in% morphData$ringnr]
    dimnames(outerGRM)[[1]] = dimnames(outerGRM)[[2]] =
      morphData[match(colnames(outerGRM), morphData$ringnr), "ID"]
    
    # Other
    otherGRM = otherGRM[colnames(otherGRM) %in% morphData$ringnr,
                        colnames(otherGRM) %in% morphData$ringnr]
    dimnames(otherGRM)[[1]] = dimnames(otherGRM)[[2]] =
      morphData[match(colnames(otherGRM), morphData$ringnr), "ID"]
    
    morphData$ID2 = ifelse(morphData$innerPi > 0,
                           morphData$ID2, NA)
    
    morphData$ID3 = ifelse(morphData$outerPi > 0,
                           morphData$ID3, NA)
    
    morphData$ID4 = ifelse(morphData$otherPi > 0,
                           morphData$ID4, NA)
  } else {
    pedigreeMap$ID2 = 1:nrow(pedigreeMap)
    numRel = dim(innerGRM)[1]
    morphData$ID2 = morphData$ID3 = morphData$ID4 =
      pedigreeMap[match(morphData$ringnr, pedigreeMap$ringnr),
                  "ID2"]
  }
}

###################### Make Pos def ###############################

# Add small value to diagonal of genome-based GRMs
# and check if matrices are positive semi definite

if (!usePed) {
  N = dim(innerGRM)[1]
  innerGRM = innerGRM + diag(1e-12, N, N)
  outerGRM = outerGRM + diag(1e-12, N, N)
  otherGRM = otherGRM + diag(1e-12, N, N)
  
  stopifnot(all(eigen(innerGRM, only.values = TRUE)$values > 0))
  stopifnot(all(eigen(outerGRM, only.values = TRUE)$values > 0))
  stopifnot(all(eigen(otherGRM, only.values = TRUE)$values > 0))
} else {
  N = dim(A_inner)[1]
  # A_inner = A_inner + diag(1e-12, N, N)
  # A_outer = A_outer + diag(1e-12, N, N)
  # A_other = A_other + diag(1e-12, N, N)
  
  stopifnot(all(eigen(A_inner, only.values = TRUE)$values > 0))
  stopifnot(all(eigen(A_outer, only.values = TRUE)$values > 0))
  stopifnot(all(eigen(A_other, only.values = TRUE)$values > 0))
}



##################### INLA Animal model ###########################
if (shrink) {
  val = "as.numeric(colnames(innerCmatrix))"
} else {
  val = "1:numRel"
}

# Define INLA formula

# Fixed effects
effects = c("sex", "FGRM", "month", "age")

# Pedigree-based or genome-based genetic effects?
if (usePed) {
  effects = c(effects, "outer", "other")
} else {
  effects = c(effects, "outerPi", "otherPi")
}

# Random effects
randomEffects = c("f(hatchYear, model = \"iid\", 
                  hyper = list(prec = list(initial = log(10), 
                  prior = \"pc.prec\", param = c(1, 0.05))))",
                  
                  "f(islandCurrent, model = \"iid\", 
                  hyper = list(prec = list(initial = log(10), 
                  prior = \"pc.prec\", param = c(1, 0.05))))",
                  
                  "f(ID1, model = \"iid\", 
                  hyper = list(prec = list(initial = log(1),
                  prior = \"pc.prec\", param = c(1, 0.05))))",
                  
                  paste0("f(ID2,values =", val, ", 
                  model = \"generic0\", Cmatrix = innerCmatrix,
                  constr = TRUE, hyper = list(
                  prec = list(initial = log(0.50), 
                  prior = \"pc.prec\",
                  param = c(1, 0.05))))"),
                  
                  paste0("f(ID3,values = ", val, ",
                  model = \"generic0\", Cmatrix = outerCmatrix,
                  constr = TRUE, hyper = list(
                  prec = list(initial = log(0.50),
                  prior = \"pc.prec\",
                  param = c(1, 0.05))))"),
                  
                  paste0("f(ID4, values = ", val, ",
                  model = \"generic0\", Cmatrix = otherCmatrix,
                  constr = TRUE, hyper = list(
                  prec = list(initial = log(0.50),
                  prior = \"pc.prec\",
                  param = c(1, 0.05))))"))

# Only when segregation variances are fit:
if (segregation) {
  morphData$ID7 = morphData$ID6 = morphData$ID5 = morphData$ID4 
  
  # Definition of S
  innerSeg = (-innerDelta + outerDelta + otherDelta) / 2
  outerSeg = (-outerDelta + innerDelta + otherDelta) / 2
  otherSeg = (-otherDelta + innerDelta + outerDelta) / 2
  
  
  # Shrink matrices
  innerSeg = innerSeg[colnames(innerSeg) %in% morphData$ringnr,
                      colnames(innerSeg) %in% morphData$ringnr]
  dimnames(innerSeg)[[1]] = dimnames(innerSeg)[[2]] =
    morphData[match(colnames(innerSeg), morphData$ringnr), "ID"]

  outerSeg = outerSeg[colnames(outerSeg) %in% morphData$ringnr,
                      colnames(outerSeg) %in% morphData$ringnr]
  dimnames(outerSeg)[[1]] = dimnames(outerSeg)[[2]] =
    morphData[match(colnames(outerSeg), morphData$ringnr), "ID"]

  otherSeg = otherSeg[colnames(otherSeg) %in% morphData$ringnr,
                      colnames(otherSeg) %in% morphData$ringnr]
  dimnames(otherSeg)[[1]] = dimnames(otherSeg)[[2]] =
    morphData[match(colnames(otherSeg), morphData$ringnr), "ID"]
  
  #innerSeg = innerSeg + diag(1e-12, N, N)
  #outerSeg = outerSeg + diag(1e-12, N, N)
  #otherSeg = otherSeg + diag(1e-12, N, N)
  
  # Check if matrices are semi def (they are not)
  stopifnot(all(eigen(innerSeg, only.values = TRUE)$values > 0))
  stopifnot(all(eigen(outerSeg, only.values = TRUE)$values > 0))
  stopifnot(all(eigen(otherSeg, only.values = TRUE)$values > 0))
  
  innerSegInv = solve(innerSeg)
  outerSegInv = solve(outerSeg)
  otherSegInv = solve(otherSeg)
  
  # Add more random effects
  randomEffects = c(randomEffects, 
                    "f(ID5,values=as.numeric(colnames(innerSegInv)), 
                  model = \"generic0\", Cmatrix = innerSegInv,
                  constr = TRUE, hyper = list(
                  prec = list(initial = log(1.2), 
                  prior = \"pc.prec\",
                  param = c(2.5, 0.05))))",
                    
                    "f(ID6,values=as.numeric(colnames(outerSegInv)),
                  model = \"generic0\", Cmatrix = outerSegInv,
                  constr = TRUE, hyper = list(
                  prec = list(initial = log(1.2),
                  prior = \"pc.prec\",
                  param = c(2.5, 0.05))))",
                    
                    "f(ID7,values=as.numeric(colnames(otherSegInv)),
                  model = \"generic0\", Cmatrix = otherSegInv,
                  constr = TRUE, hyper = list(
                  prec = list(initial = log(1.2),
                  prior = \"pc.prec\",
                  param = c(2.5, 0.05))))")
}

# Add random effects to fixed effects
effects = c(effects, randomEffects)

# Make formula
s = reformulate(effects, response = response)

# Find inverse of G (or A)
if (usePed) {
  # Round for INLA-reasons...
  innerCmatrix = round(solve(A_inner), 5) # round(AInv_inner, 5)
  outerCmatrix = round(solve(A_outer), 5) # round(AInv_outer, 5)
  otherCmatrix = round(solve(A_other), 5) # round(AInv_other, 5)
  
  dimnames(innerCmatrix) = dimnames(A_inner)
  dimnames(outerCmatrix) = dimnames(A_outer)
  dimnames(otherCmatrix) = dimnames(A_other)
  
  innerCmatrix[innerCmatrix == 0] = 0
  outerCmatrix[outerCmatrix == 0] = 0
  otherCmatrix[otherCmatrix == 0] = 0
} else {
  
  innerCmatrix = solve(innerGRM)
  if (!isSymmetric(innerCmatrix)){
    innerCmatrix = forceSymmetric(innerCmatrix)
  }
  
  outerCmatrix = solve(outerGRM)
  if (!isSymmetric(outerCmatrix)){
    outerCmatrix = forceSymmetric(outerCmatrix)
  }
  
  otherCmatrix = solve(otherGRM)
  if (!isSymmetric(otherCmatrix)){
    otherCmatrix = forceSymmetric(otherCmatrix)
  }
}

# First INLA run
model = inla(formula = s, family = "gaussian",
             data = morphData, verbose = TRUE,
             control.inla = list(h = h_val),
             control.family = list(
               hyper = list(theta = list(initial = log(1),
                                         prior = "pc.prec",
                                         param=c(1, 0.05)))))#,
#             control.compute=list(dic = T, config = TRUE))
if (iterate) {
  model01 = model
  model02 = inla.rerun(model01)
  model = inla.rerun(model02)
}
  
# Transform precisions to variances and store for easy access
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

# Save model results with an informative name
path = "Runs/inla_result_hetGG"

if (usePed) {
  path = paste0(path, "_pedigree_")
} else {
  path = paste0(path, "_rio", loter_run,
                "_", as.character(GRMVariant), "_")
}

if (segregation) {
  path = paste0(path, "_segregation_")
}

if (!shrink)
  path = paste0(path, "_full_")

if (iterate) {
  save(model01, model02, model, 
       file = paste0(path, response, "h=", h_val, ".RData"))
} else {
  save(model, 
       file = paste0(path, response, "h=", h_val, ".RData"))
}
