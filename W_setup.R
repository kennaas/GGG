if (!require(R.utils)) install.packages("R.utils")
library(R.utils)
if (!require(data.table)) install.packages("data.table")
library(data.table)

source("My_R_code/file_backed_mat.R")

loter_run = 1

# Load phased haplotype data from different populations

print("Loading data")
admixedVCF = 
  fread(file = paste0("Data/loter/Run ", loter_run,
               "/Loter input ", loter_run,
               "/admixed_phased", loter_run, ".vcf.gz"),
        key = "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT")
innerVCF = 
  fread(file = paste0("Data/loter/Run ", loter_run,
               "/Loter input ", loter_run, 
               "/inner_phased", loter_run, ".vcf.gz"),
        key = "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT")
outerVCF = 
  fread(file = paste0("Data/loter/Run ", loter_run,
                      "/Loter input ", loter_run,
                      "/outer_phased", loter_run, ".vcf.gz"),
        key = "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT")
otherVCF = 
  fread(file = paste0("Data/loter/Run ", loter_run,
                      "/Loter input ", loter_run, 
                      "/other_phased", loter_run, ".vcf.gz"),
        key = "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT")

# Merge into one data.table (still on vcf format)
fullVCF = 
  merge(innerVCF, merge(outerVCF, merge(otherVCF, admixedVCF)))
print("Merge complete.")

inds = colnames(fullVCF)[-(1:9)]
numSNPs = length(fullVCF$ID)

print(numSNPs)

rm(innerVCF)
rm(outerVCF)
rm(otherVCF)
rm(admixedVCF)

# Indices
SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

# Create file-backed matrix to hold haplotypes
W = BGData(pheno = data.frame(ID = inds),
           geno = initFileBackedMatrix(length(inds),
                                       2 * numSNPs,
                                       outputType = "logical",
                                       folderOut = paste0(
                                         "Data/loter/Run ", loter_run, "/W")))

# Fill matrix with haplotype data from the VCF format
for (ind in 10:(9 + length(inds))) {
  W@geno[ind - 9, SNPs1] = sapply(fullVCF[, ..ind],
                             function(x) as.numeric(grepl("^1", x)))
  W@geno[ind - 9, SNPs2] = sapply(fullVCF[, ..ind],
                             function(x) as.numeric(grepl("1$", x)))
}

# Give columns reasonable names
colnames(W@geno) = do.call("paste0",
                           c(list(rep(fullVCF$ID, each = 2)),
                             "_",
                             list(1:2)))

# Save
save(W, file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))

