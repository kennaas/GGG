if (!require(R.utils)) install.packages("R.utils")
library(R.utils)
if (!require(data.table)) install.packages("data.table")
library(data.table)

source("My_R_code/file_backed_mat.R")

loter_run = 4

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

SNPs1 = seq(1, by = 2, to = 2 * numSNPs)
SNPs2 = seq(2, by = 2, to = 2 * numSNPs)

W = BGData(pheno = data.frame(ID = inds),
           geno = initFileBackedMatrix(length(inds),
                                       2 * numSNPs,
                                       outputType = "logical",
                                       folderOut = paste0(
                                         "Data/loter/Run ", loter_run, "/W")))

for (ind in 10:(9 + length(inds))) {
  W@geno[ind - 9, SNPs1] = sapply(fullVCF[, ..ind],
                             function(x) as.numeric(grepl("^1", x)))
  W@geno[ind - 9, SNPs2] = sapply(fullVCF[, ..ind],
                             function(x) as.numeric(grepl("1$", x)))
}

colnames(W@geno) = do.call("paste0",
                           c(list(rep(fullVCF$ID, each = 2)),
                             "_",
                             list(1:2)))

save(W, file = paste0("Data/loter/Run ", loter_run, "/W/W.RData"))

#################### Old way, writing to txt files ################

# print("Making W1 matrix")
# W1 = apply(fullVCF[, -(1:9)], c(1,2), function(x) as.numeric(grepl("^1", x)))
# print("W1 created.")
# rownames(W1) = fullVCF$ID
# W1 = t(W1)
# print("W1 transposed. Adding ID column.")
# W1 = cbind(data.frame(ID = inds), as.data.frame(W1))
# print("Writing W1 table.")
# write.table(W1, file = "Data/W1_beagle3.txt",
#             quote = FALSE, row.names = FALSE)
# rm(W1)
# print("W1 done")
# 
# print("Making W2 matrix")
# W2 = apply(fullVCF[, -(1:9)], c(1,2), function(x) as.numeric(grepl("1$", x)))
# print("W2 created.")
# rownames(W2) = fullVCF$ID
# rm(fullVCF)
# W2 = t(W2)
# print("W2 transposed. Adding ID column.")
# W2 = cbind(data.frame(ID = inds), as.data.frame(W2))
# print("Writing W2 table.")
# write.table(W2, file = "Data/W2_beagle3.txt",
#             quote = FALSE, row.names = FALSE)
# print("W2 done")

# save(W1, W2, file = "Data/GG_W_matrices.RData")

######################## Packaged version #########################

# Won't work, wrong coding. 2 is both 1/0 and 0/1

# library(VariantAnnotation)
# 
# inner = readVcf(file = "Data/loter/Loter input 1/inner_phased.vcf.gz")
# outer = readVcf(file = "Data/loter/Loter input 1/outer_phased.vcf.gz")
# other = readVcf(file = "Data/loter/Loter input 1/other_phased.vcf.gz")
# admixed = readVcf(file = "Data/loter/Loter input 1/admixed_phased.vcf.gz")
# 
# inner = genotypeToSnpMatrix(inner)
# outer = genotypeToSnpMatrix(outer)
# other = genotypeToSnpMatrix(other)
# admixed = genotypeToSnpMatrix(admixed)
# 
# full = rbind2(inner, outer, other, admixed)
# 
# W1 = ifelse(full$genotypes == 1, 1, NA)


