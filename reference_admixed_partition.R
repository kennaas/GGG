# Creating a new founder pop:
# Use only natalData to assign individuals to reference or admixed

library(BGData)

load(file = "Runs/morphData_pedigree_version.RData")
load(file = "Runs/pedigree.RData")
load(file = "Runs/natalData_ped.RData")

head(natalData$ringnr)
#founders = with(morphData, ringnr[hatchYear == 2007])

phenotypedNatalData = natalData[match(morphData$ringnr, natalData$ringnr), ]

innerInds = with(phenotypedNatalData, na.omit(unique(ringnr[islandIO == "inner"])))
outerInds = with(phenotypedNatalData, na.omit(unique(ringnr[islandIO == "outer"])))
otherInds = with(phenotypedNatalData, na.omit(unique(ringnr[islandIO == "other"])))

admixedInds = unique(morphData$ringnr[!morphData$ringnr %in%
                                        c(innerInds, outerInds, otherInds)])

load.BGData(file =
              "Data/BGData/BGData_Helgeland_01_2018_Dosage/BGData.RData")
genotypedInds = dimnames(BGData@geno)[[1]]
# All natal inds were genotyped
all(natalData$ringnr %in% genotypedInds)

## Founders4: only phenotyped and genotyped, based only on natal##


write.table(unique(genotypedInds[!genotypedInds %in% innerInds]),
            file = "Data/Founders4/NOT_inner_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(genotypedInds[!genotypedInds %in% outerInds]),
            file = "Data/Founders4/NOT_outer_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(genotypedInds[!genotypedInds %in% otherInds]),
            file = "Data/Founders4/NOT_other_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(genotypedInds[!genotypedInds %in% admixedInds]),
            file = "Data/Founders4/NOT_admixed_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

save(admixedInds, innerInds, otherInds, outerInds,
     file = "Data/Founders4/founders.RData")


################### Founders 3: only based on natal, include unphenotyped ########

# write.table(genotypedInds[!genotypedInds %in% innerInds],
#             file = "Data/Founders3/NOT_inner_inds.txt",
#             row.names = FALSE, col.names = FALSE,
#             quote = FALSE)
# 
# write.table(genotypedInds[!genotypedInds %in% outerInds],
#             file = "Data/Founders3/NOT_outer_inds.txt",
#             row.names = FALSE, col.names = FALSE,
#             quote = FALSE)
# 
# write.table(genotypedInds[!genotypedInds %in% otherInds],
#             file = "Data/Founders3/NOT_other_inds.txt",
#             row.names = FALSE, col.names = FALSE,
#             quote = FALSE)
# 
# admixedInds = genotypedInds[
#   !genotypedInds %in% c(innerInds, outerInds, otherInds)]
# 
# write.table(genotypedInds[!genotypedInds %in% admixedInds],
#             file = "Data/Founders3/NOT_admixed_inds.txt",
#             row.names = FALSE, col.names = FALSE,
#             quote = FALSE)
# 
# save(admixedInds, unique(innerInds), unique(otherInds), unique(outerInds),
#      file = "Data/Founders3/founders.RData")
