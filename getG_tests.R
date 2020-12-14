compute_maf = function(genotype) {
  # Frequency of the counted allele, i.e. the allele appended 
  # after the SNP name. 
  allele_freq = summarize(genotype,
                          nCores = nCores, verbose = TRUE)[, 2]
  
  # Minor allele frequency
  return(ifelse(allele_freq < 0.5, allele_freq, 1 - allele_freq))
}


# Number of individuals
N = dim(geno(BGData))[1]
# Number of SNPs
m = dim(geno(BGData))[2]
nCores = 1

load.BGData(file = "BGData_Helgeland_01_2018_Dosage/BGData.RData")

inds = sample(1:N, 6, replace = FALSE)
SNPs = sample(1:m, 6, replace = FALSE)
genotype = geno(BGData)[inds,SNPs]


# default
getG(X = genotype, verbose = TRUE, nCores = nCores)

# as. ok. That is, m inside or outside is the same
# (1 + getG(X = genotype, verbose = TRUE, 
#           scaleG = FALSE, center = rep(1, 6), 
#           scale = rep(sqrt(6), 6), nCores = nCores)) / 2
# 
# (1 + getG(X = genotype, verbose = TRUE, 
#           scaleG = FALSE, center = rep(1, 6), 
#           scale = FALSE, nCores = nCores) / 6) / 2

#vr
maf = compute_maf(genotype)
# What to divide each entry by
scaling = 2 * sum(maf * (1 - maf))
getG(X = genotype, verbose = TRUE, scaleG = FALSE,
           center = 2 * maf, 
           scale = FALSE,
           nCores = nCores) / scaling

G = getG(X = genotype, verbose = TRUE, scaleG = FALSE,
           center = 2 * maf, 
           scale = rep(sqrt(scaling), 6),
           nCores = nCores)
