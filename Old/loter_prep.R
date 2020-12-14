library(data.table)

load(file = "Runs/morphData_pedigree_version.RData")
load(file = "Runs/pedigree.RData")
map = fread(file = "Data/Helgeland_01_2018.map")
# load.BGData(file = "Data/BGData_Helgeland_01_2018_Dosage/BGData.RData")

outerRef = pedigreeData$ringnr[pedigreeData$outer == 1]
innerRef = pedigreeData$ringnr[pedigreeData$inner == 1]
otherRef = pedigreeData$ringnr[pedigreeData$other == 1]
admixedRef = pedigreeData$ringnr[
  !pedigreeData$ringnr %in% c(otherRef, outerRef, innerRef)]

# fullVCF = fread(file = "Data/beagle_out/beagle_out.vcf")
# 
# admixedIndices = c(1:9, which(colnames(fullVCF) %in% admixedRef))
# admixedVCF = fullVCF[, ..admixedIndices]
# 
# outerIndices = c(1:9, which(colnames(fullVCF) %in% outerRef))
# outerVCF = fullVCF[, ..outerIndices]
# 
# innerIndices = c(1:9, which(colnames(fullVCF) %in% innerRef))
# innerVCF = fullVCF[, ..innerIndices]
# 
# otherIndices = c(1:9, which(colnames(fullVCF) %in% otherRef))
# otherVCF = fullVCF[, ..otherIndices]
# 
# write.table(admixedVCF, file = "Data/beagle_out/admixed.vcf", 
#             quote = FALSE, row.names = FALSE)
# write.table(outerVCF, file = "Data/beagle_out/outer.vcf", 
#             quote = FALSE, row.names = FALSE)
# write.table(innerVCF, file = "Data/beagle_out/inner.vcf", 
#             quote = FALSE, row.names = FALSE)
# write.table(otherVCF, file = "Data/beagle_out/other.vcf", 
#             quote = FALSE, row.names = FALSE)
# 
# admixed = fread(file = "Data/beagle_out/admixed.vcf")
# outer = fread(file = "Data/beagle_out/outer.vcf")
# inner = fread(file = "Data/beagle_out/inner.vcf")
# other = fread(file = "Data/beagle_out/other.vcf")
# full = fread(file = "Data/beagle_out/beagle_out.vcf")

unAdmixed = fread(file = "Beagle/Data_from_henrik/NOT_admixed_inds.txt")
uninner = fread(file = "Beagle/Data_from_henrik/NOT_inner_inds.txt")
unOuter = fread(file = "Beagle/Data_from_henrik/NOT_outer_inds.txt")
unOther = fread(file = "Beagle/Data_from_henrik/NOT_other_inds.txt")