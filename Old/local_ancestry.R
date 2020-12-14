######################## Libraries ################################
if (!require(EILA)) install.packages("EILA")
library(EILA)
if (!require(BGData)) install.packages("BGData")
library(BGData)
if (!require(data.table)) install.packages("data.table")
library(data.table)

######################## Load Data ################################
load(file = "Runs/morphData_pedigree_version.RData")
load(file = "Runs/pedigree.RData")
map = fread(file = "Data/Helgeland_01_2018.map")
# load.BGData(file = "Data/BGData_Helgeland_01_2018_Dosage/BGData.RData")


outerRef = pedigreeData$ringnr[pedigreeData$outer == 1]
innerRef = pedigreeData$ringnr[pedigreeData$inner == 1]
otherRef = pedigreeData$ringnr[pedigreeData$other == 1]
admixedRef = pedigreeData$ringnr[
  !pedigreeData$ringnr %in% c(otherRef, outerRef, innerRef)]

###################### Unimputed ##################################

# One time thing to get BGData objects for each subpopulation
# unimputedRaw = fread(file = "Data/Helgeland_01_2018_Dosage.txt")
# write.table(fullRaw[IID %in% outerRef, ], 
#             quote = FALSE, row.names = FALSE, 
#             file = "Data/Helgeland_01_2018_Dosage_outer.txt")
# write.table(fullRaw[IID %in% innerRef, ], 
#             quote = FALSE, row.names = FALSE, 
#             file = "Data/Helgeland_01_2018_Dosage_inner.txt")
# write.table(fullRaw[IID %in% otherRef, ], 
#             quote = FALSE, row.names = FALSE, 
#             file = "Data/Helgeland_01_2018_Dosage_other.txt")
# write.table(fullRaw[IID %in% admixedRef, ], 
#             quote = FALSE, row.names = FALSE, 
#             file = "Data/Helgeland_01_2018_Dosage_admixed.txt")
# readRAW(fileIn = "Data/Helgeland_01_2018_Dosage_outer.txt", 
#         idCol = 2, verbose = TRUE)
# readRAW(fileIn = "Data/Helgeland_01_2018_Dosage_inner.txt", 
#         idCol = 2, verbose = TRUE)
# readRAW(fileIn = "Data/Helgeland_01_2018_Dosage_other.txt", 
#         idCol = 2, verbose = TRUE)
# readRAW(fileIn = "Data/Helgeland_01_2018_Dosage_admixed.txt", 
#         idCol = 2, verbose = TRUE)

# load.BGData(
#   file = "Data/BGData_Helgeland_01_2018_Dosage_outer/BGData.RData")
# outerBGData = BGData
# load.BGData( # this one must be reread from raw, oops
#   file = "Data/BGData_Helgeland_01_2018_Dosage_inner/BGData.RData")
# innerBGData = BGData
# load.BGData(
#   file = "Data/BGData_Helgeland_01_2018_Dosage_other/BGData.RData")
# otherBGData = BGData
# load.BGData(
#   file = "Data/BGData_Helgeland_01_2018_Dosage_admixed/BGData.RData")
# admixedBGData = BGData
# rm(BGData)


######################### Beagle imputed data #####################

# imputedRaw = fread(file = "Data/beagle_out_raw/beagle_out.raw")

# readRAW(file = "Data/beagle_out_raw/beagle_out.raw",
#         idCol = 2, verbose = TRUE)

# write.table(imputedRaw[IID %in% outerRef, ], 
#             quote = FALSE, row.names = FALSE, 
#             file = "Data/imputed_dosage_outer.txt")
# write.table(imputedRaw[IID %in% innerRef, ],
#             quote = FALSE, row.names = FALSE,
#             file = "Data/imputed_dosage_inner.txt")
# write.table(imputedRaw[IID %in% otherRef, ],
#             quote = FALSE, row.names = FALSE,
#             file = "Data/imputed_dosage_other.txt")
# write.table(imputedRaw[IID %in% admixedRef, ],
#             quote = FALSE, row.names = FALSE,
#             file = "Data/imputed_dosage_admixed.txt")
# readRAW(fileIn = "Data/imputed_dosage_outer.txt", 
#         idCol = 2, verbose = TRUE)
# readRAW(fileIn = "Data/imputed_dosage_inner.txt",
#         idCol = 2, verbose = TRUE)
# readRAW(fileIn = "Data/imputed_dosage_other.txt",
#         idCol = 2, verbose = TRUE)
# readRAW(fileIn = "Data/imputed_dosage_admixed.txt",
#         idCol = 2, verbose = TRUE)



# load.BGData(
#   file = "Data/BGData_imputed_dosage_outer/BGData.RData")
# outerBGData = BGData
# load.BGData( # this one must be reread from raw, oops
#   file = "Data/BGData_imputed_dosage_inner/BGData.RData")
# innerBGData = BGData
# load.BGData(
#   file = "Data/BGData_imputed_dosage_other/BGData.RData")
# otherBGData = BGData
# load.BGData(
#   file = "Data/BGData_imputed_dosage_admixed/BGData.RData")
# admixedBGData = BGData
load.BGData(file = "Data/BGData_beagle_out/BGData.RData")

unwanted_markers = read.table("Data/markers_on_0_or_30.txt")$V1
map_reduced = map[!map$V2 %in% unwanted_markers, ]

######################### EILA ####################################

# Chromosome-wise LA

chromosomes = unique(map_reduced$V1)

all_markers = dimnames(geno(BGData))[[1]]

for (chr in chromosomes) {
  whichChr = which(map_reduced$V1 == chr)
  whichAdmixed = which(dimnames(geno(BGData))[[1]] %in% admixedRef)
  whichInner = which(dimnames(geno(BGData))[[1]] %in% innerRef)
  whichOuter = which(dimnames(geno(BGData))[[1]] %in% outerRef)
  whichOther = which(dimnames(geno(BGData))[[1]] %in% otherRef)

  la = eila(admixed = t(geno(BGData)[whichAdmixed, whichChr]),
            anc1 = t(geno(BGData)[whichInner, whichChr]),
            anc2 = t(geno(BGData)[whichOuter, whichChr]),
            anc3 = t(geno(BGData)[whichOther, whichChr]),
            position = map_reduced$V4[whichChr])
  print(chr)
  save(la, file = paste0("Runs/la/la_chr,", chr, ".RData"))
}

