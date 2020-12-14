##################### Libraries ###################################
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(optparse)) install.packages("optparse")
library(optparse)

##################### Options #####################################


option_list = list(
  make_option(c("-g", "--GRMVariant"), type="character", default="vanRaden", 
              help="which definition of G", metavar="character"),
  make_option(c("-a", "--alpha"), type = "integer", default = -1,
              help="alpha parameter in GCTA", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# BGData_default, allele_sharing, vanRaden, GCTA
GRMVariant = opt$GRMVariant
# Yang has alpha = -1
alpha = opt$alpha

# Windows does not support multiple cores 
OS = Sys.info()[['sysname']]
nCores = ifelse(OS == "Windows", 1, 8)

##################### Data import #################################

# Instead:
# The dosage-file is in .raw format from PLINK.
# See https://www.cog-genomics.org/plink/1.9/formats#raw

# First time: load using .raw file to a BGData object.
#dosageData = readRAW(fileIn = "Data/Helgeland_01_2018_Dosage.txt",
#                     idCol = 2)

# After first time, load the data now in BGData format.
load.BGData(file = "Data/BGData_Helgeland_01_2018_Dosage/BGData.RData")

################### Minor Allele Frequency Function ###############

# Returns the minor allele frequency based on a genotype matrix
compute_maf = function(genotype) {
  # Frequency of the counted allele, i.e. the allele appended 
  # after the SNP name. 
  allele_freq = summarize(genotype,
                          nCores = nCores, verbose = TRUE)[, 2]
  
  # Minor allele frequency
  return(ifelse(allele_freq < 0.5, allele_freq, 1 - allele_freq))
}

################### Computing GRM #################################
# Can do custom centerings and scalings for different version of G. 

# Number of individuals
N = dim(geno(BGData))[1]
# Number of SNPs
m = dim(geno(BGData))[2]

switch (GRMVariant,
  BGData_default = {
    # Default scaling and centering in BGData
    GRM = getG(X = geno(BGData), verbose = TRUE, nCores = nCores)
  },
  allele_sharing = {
    # Allele-sharing coefficient, eq. (6) in Speed & Balding 2014
    # Gives strange (?) result.
    GRM =
      (1 + getG(X = geno(BGData), verbose = TRUE, scaleG = FALSE,
                center = rep(1, m), scale = FALSE, nCores = nCores)
       / m) / 2
  },
  vanRaden = {
    # maf = compute_maf(geno(BGData))
    allele_freq = summarize(geno(BGData),
                            nCores = nCores, verbose = TRUE)[, 2]
    
    # What to divide each entry by
    scaling = 2 * sum(allele_freq * (1 - allele_freq))
    
    GRM = getG(X = geno(BGData), verbose = TRUE, scaleG = FALSE,
               center = 2 * allele_freq, 
               scale = FALSE,
               nCores = nCores) / scaling
  },
  GCTA = {
    # Mark the file with the alpha value
    GRMVariant = paste0(GRMVariant, "_alpha_", as.character(alpha))
    
    allele_freq = summarize(geno(BGData),
                            nCores = nCores, verbose = TRUE)[, 2]
    
    #maf = compute_maf(geno(BGData))
    # What to divide each entry by
    scaling = (2 * allele_freq * (1 - allele_freq)) ^
      (-as.numeric(alpha) / 2)
    
    GRM = getG(X = geno(BGData), verbose = TRUE, scaleG = FALSE,
               center = 2 * allele_freq,
               scale = scaling,
               nCores = nCores) / m
  }
)

################### Save Result ###################################

rm(BGData)
save(GRM, file = paste0("Runs/GRMs/GRM_", GRMVariant, ".RData"))