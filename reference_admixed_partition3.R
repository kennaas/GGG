# Use only natalData to assign individuals to reference or admixed

library(BGData)

# Load data
load(file = "Runs/morphData_pedigree_version.RData")
load(file = "Runs/pedigree.RData")
load(file = "Runs/natalData_ped.RData")

# Natal data for phenotyped individuals
phenotypedNatalData = natalData[match(morphData$ringnr, natalData$ringnr), ]

# Reference populations
innerInds = with(phenotypedNatalData, 
                 na.omit(unique(ringnr[islandIO == "inner"])))
outerInds = with(phenotypedNatalData, 
                 na.omit(unique(ringnr[islandIO == "outer"])))
otherInds = with(phenotypedNatalData, 
                 na.omit(unique(ringnr[islandIO == "other"])))

# Get genotype data
load.BGData(
  file="Data/BGData/BGData_Helgeland_01_2018_Dosage/BGData.RData")
# Genotyped individuals
genotypedInds = dimnames(BGData@geno)[[1]]
# All natal inds were genotyped
all(natalData$ringnr %in% genotypedInds)

# Admixed population
admixedInds = genotypedInds[!genotypedInds %in% c(innerInds, outerInds, otherInds)]

# Write files excluding each population for use in Beagle
write.table(unique(genotypedInds[!genotypedInds %in% innerInds]),
            file = "Data/Founders3/NOT_inner_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(genotypedInds[!genotypedInds %in% outerInds]),
            file = "Data/Founders3/NOT_outer_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(genotypedInds[!genotypedInds %in% otherInds]),
            file = "Data/Founders3/NOT_other_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(genotypedInds[!genotypedInds %in% c(innerInds,
                                                       outerInds, 
                                                       otherInds)]),
            file = "Data/Founders3/NOT_admixed_inds.txt",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# Save partition for future reference.
save(admixedInds, unique(innerInds), unique(otherInds),
     unique(outerInds), file = "Data/Founders3/founders.RData")
