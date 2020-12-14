##################### Libraries ###################################
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(optparse)) install.packages("optparse")
library(optparse)
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(R.utils)) install.packages("R.utils")
library(R.utils)

##################### Options #####################################

# Windows does not support multiple cores 
OS = Sys.info()[['sysname']]
nCores = ifelse(OS == "Windows", 1, 8)

##################### Data import #################################

# Instead:
# The dosage-file is in .raw format from PLINK.
# See https://www.cog-genomics.org/plink/1.9/formats#raw


admixed_vcf = fread(file = "Data/loter/Loter input 1/admixed_phased.vcf.gz")
# outer_vcf = fread(file = "Data/loter/Loter input 1/outer_phased.vcf.gz")
# inner_vcf = fread(file = "Data/loter/Loter input 1/inner_phased.vcf.gz")
# other_vcf = fread(file = "Data/loter/Loter input 1/other_phased.vcf.gz")

inner_w = fread(file = "Data/loter/Loter input 1/inner_phased_dosage/inner_loter_input_recoded.raw")
other_w = fread(file = "Data/loter/Loter input 1/other_phased_dosage/other_loter_input_recoded.raw")
outer_w = fread(file = "Data/loter/Loter input 1/outer_phased_dosage/outer_loter_input_recoded.raw")
admixed_w = fread(file = "Data/loter/Loter input 1/admixed_phased_dosage/admixed_loter_input_recoded.raw")

W = rbind(inner_w, other_w, outer_w, admixed_w)

refLAOuter = fread(file = "Data/loter/refLAOuter.txt")



# After first time, load the data now in BGData format.
load.BGData(file = "BGData_refLAOuter/BGData.RData")
load.BGData(file = "Data/BGData/BGData_Helgeland_01_2018_Dosage/BGData.RData")

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
          maf = compute_maf(geno(BGData))
          # What to divide each entry by
          scaling = 2 * sum(maf * (1 - maf))
          
          GRM = getG(X = geno(BGData), verbose = TRUE, scaleG = FALSE,
                     center = 2 * maf, 
                     scale = FALSE,
                     nCores = nCores) / scaling
        },
        vanRaden2 = {
          maf = compute_maf(geno(BGData))
          # What to divide each entry by
          scaling = 2 * sum(maf * (1 - maf))
          
          GRM = getG(X = geno(BGData), verbose = TRUE, scaleG = FALSE,
                     center = 1 + 2 * maf, 
                     scale = FALSE,
                     nCores = nCores) / scaling
        },
        generalized_yang = {
          # Mark the file with the alpha value
          GRMVariant = paste0(GRMVariant, "_alpha_", as.character(alpha))
          
          maf = compute_maf(geno(BGData))
          # What to divide each entry by
          scaling = (2 * maf * (1 - maf)) ^ (-as.numeric(alpha) / 2)
          
          GRM = getG(X = geno(BGData), verbose = TRUE, scaleG = FALSE,
                     center = 2 * maf,
                     scale = scaling,
                     nCores = nCores) / m
        }, 
        RAE = {
          
        }
        
)

################### Save Result ###################################

rm(BGData)
save(GRM, file = paste0("Runs/GRM_", GRMVariant, ".RData"))