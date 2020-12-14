####################### Libraries #################################
if (!require(nadiv)) install.packages('nadiv')
  library(nadiv)

####################### Options ###################################

usePed = TRUE
center = TRUE

####################### Data Importing ############################

# We load the sparrow data.

# The data table `yearData` contains the year data for individuals 
# (hatchyear and last year observed), 
# as well as last and first island observed.
yearData =
  read.csv(
    file =
      "data/spuRvebasen-hatchyear-island_1993-2016_february2017_forSteffi.csv")
names(yearData)[2:5] = c("hatchYear", "lastYear", 
                      "firstIsland", "lastIsland")

# The data table `morphData`contains measured
# (including repeated measurements) morphological data of 
# Helgeland sparrows in the period 1993-2016.
morphData = read.csv(
  file =
    "data/RepeatedAdMorph-200K-Helgeland_20170623_extrainfo.csv")
names(morphData)[c(5, 7:8, 18:19)] =
  c("islandCurrent", "hatchIsland", "hatchYear",
    "totalBadge", "visibleBadge")

# Natal island data
natalData = read.table("data/Dispersaldata_20181121.txt", 
                  header = TRUE, sep = " ") 
names(natalData)[2:3] = c("natalIsland", "adultIsland")


# Load the morph data, where inbreeding coefficients F (genomic?)
# should be ok (as per 7.11.17)
FData = read.table(
  file = "data/steffie_inbreeding_morph_data_071117.txt",
  header = TRUE, stringsAsFactors = FALSE)

# A map of Ringnr and ID for the d.morph file
mapRingnrID = read.csv(file =
                         "data/Ringnr_ID_Link_N3147_SNP183384.csv")
names(mapRingnrID)[1] = "ringnr"

if (usePed) {
  # Import the pedigree
  pedigreeData = read.table(
    file = "data/SNP_pedigree_Helgeland_05122017.txt",
    header = T, sep=" ")
  names(pedigreeData) = c("ringnr", "mother", "father")
}


####################### Data Prep #################################

# The following ringumbers must be removed,
# because they were duplicated:
duplicatedInds = read.table(file =
                         "data/Duplicates_to_remove_200kSNP.txt",
                         header = FALSE)
morphData = morphData[!(morphData$ringnr %in% duplicatedInds[,1]),]

# When pedigree is used
if (usePed)
  morphData = morphData[morphData$ringnr %in% pedigreeData$ringnr,]

# Set missing values to NA
morphData[morphData$tarsus < 0, "tarsus"] = NA
morphData[morphData$wing < 0, "wing"] = NA
morphData[morphData$billD < 0, "billD"] = NA
morphData[morphData$billL < 0, "billL"] = NA
morphData[morphData$mass < 0, "mass"] = NA
morphData[morphData$totalBadge < 0, "totalBadge"] = NA
morphData[morphData$visibleBadge < 0, "visibleBadge"] = NA

# Recast parts of the data into `R` factors, where it makes sense 

# Add inbreeding to morph data frame
morphData$FGRM = FData[match(morphData$ringnr, FData$IID), "FGRM"]

# Recast parts of the data into `R` factors, where it makes sense
morphData$islandCurrent = as.factor(morphData$islandCurrent)
morphData$hatchYear = as.factor(morphData$hatchYear)
morphData$year = as.factor(morphData$year)

# Encode sex as 0/1 (not 1/2)
morphData$sex = ifelse(morphData$sex == 2, 1, 0)

if (center) {
  # Center continuous data
  morphData$age = scale(morphData$age, scale = FALSE)
  morphData$month = scale(morphData$month, scale = FALSE)
  morphData$FGRM = scale(morphData$FGRM, scale = FALSE)
}


# Pick only relevant columns
#d.morph = d.morph[,c("ID","ringnr","indnr","sex","island_current",
                      # "adultSNPisland","hatchisland","hatchyear",
                      # "year","month","age","wing","mass",
                      # "island","FGRM","id")]


####################### Pedigree Setup ############################

# Assign groups based on birth island
# Islands 22, 23, 24, 88 are the "outer" islands
# Islands 20, 26, 27, 28, 88 are the "inner" islands
# The rest are "other" islands

if (usePed) {
  # Organizing the pedigree -- estimate and store relatedness.
  # Need an ordered Pedigree for the inverseA() function:
  
  pedigreeData = prepPed(pedigreeData)
  
  # Also replace the ringnumber by an id (1:#animals)
  pedigreeData$id = 1:(nrow(pedigreeData))
  pedigreeMap = pedigreeData[, c("ringnr", "id")]
  pedigreeMap$ID = 
    mapRingnrID[match(pedigreeMap$ringnr,mapRingnrID$ringnr), "ID"]
  pedigreeMap = pedigreeMap[order(pedigreeMap$ID), ]
  
  # Give mother and father the id
  pedigreeData$motherId =
    pedigreeMap[match(pedigreeData$mother, pedigreeMap$ringnr),
                "id"]
  pedigreeData$fatherId =
    pedigreeMap[match(pedigreeData$father, pedigreeMap$ringnr),
                "id"]
  
  # and again the same to keep track:
  pedigreeData$motherIdOriginal = pedigreeData$motherId
    # mapRingnrID[match(pedigreeData$mother, pedigreeMap$ringnr),
    #               "id"]
  pedigreeData$fatherIdOriginal = pedigreeData$fatherId
    # mapRingnrID[match(pedigreeData$father, pedigreeMap$ringnr), 
    #               "id"]
  
  # Use full pedigree for inversion and
  # calculation of inbreeding coefficients
  A = makeA(pedigreeData[, 1:3])
  
  # Add id to the morph and year data 
  # (used for the Ainv and A matrices from pedigree)
  morphData$id =
    pedigreeMap[match(morphData$ringnr, pedigreeMap$ringnr), "id"]
  yearData$id =
    pedigreeMap[match(yearData$ringnr, pedigreeMap$ringnr), "id"]
}

####################### Genetic Group Setup #######################

if (usePed) {
  # Categorize natal islands
  natalData$islandIO =
    ifelse(natalData$natalIsland %in% c(22:24, 88),
           "outer", 
           ifelse(natalData$natalIsland %in% c(20, 26:28, 38),
                  "inner",
                  "other"))
  
  # Categorize adult islands
  natalData$adultIslandIO =
    ifelse(natalData$adultIsland %in% c(22:24, 88),
           "outer",
           ifelse(natalData$adultIsland %in% c(20, 26:28, 38),
                  "inner",
                  "other"))
  
  # Dispersal rate between island systems:
  sum(natalData$adultIslandIO != natalData$islandIO) /
    nrow(natalData)
  # Dispersal rate between islands
  sum(natalData$disp == 1) / nrow(natalData)
  
  # Assign best knowledge of natal island. The islandIO in yearData
  # should be the genetic natal island from natalData (if 
  # available), and otherwise the `firstIsland` from yearData.
  # Create therefore a new column `natalIsland` that is initiated
  # with `firstIsland`.
  yearData$natalIsland = yearData$firstIsland
  yearData$natalIsland = 
    ifelse(yearData$ringnr %in% natalData$ringnr, 
           natalData[match(yearData$ringnr, natalData$ringnr),
                     "natalIsland"], 
           yearData$firstIsland)
  yearData$islandIO =
    ifelse(yearData$natalIsland %in% c(22:24, 88),
           "outer",
           ifelse(yearData$natalIsland %in% c(20, 26:28, 38),
                  "inner",
                  "other"))
  
  # Introduce a group variable
  pedigreeData$GG = rep(NA, nrow(pedigreeData))
  
  # Set genetic group to island where animal lived at first year
  # that it was observed. This corresponds to hatchIsland in
  # the `morphData` file but there are less individuals than in
  # the pedigree. If animal is a dummy animal, i.e. has 
  # name length 5, assign it the island of its mother or father.
  
  pedigreeData$GG =
    ifelse( # if ringnr is a real animal with 7 digits ringnr
      nchar(as.character(pedigreeData$ringnr)) == 7,
      yearData[match(pedigreeData$ringnr, yearData$ringnr),
               "islandIO"],
      # Otherwise if mother is a real animal
      ifelse(nchar(as.character(pedigreeData$mother)) == 7,
             yearData[match(pedigreeData$mother, yearData$ringnr),
                      "islandIO"],
             # Else if father is a real animal
             ifelse(nchar(as.character(pedigreeData$father)) == 7,
                    yearData[match(pedigreeData$father,
                                   yearData$ringnr), "islandIO"],
                    NA)))
  
  # Finally, there are 271 dummy animals left without parents.
  # Assign them genetic group of their offspring:
  for (i in 1:nrow(pedigreeData)) {
    if (is.na(pedigreeData$GG[i])) {
      pedigreeData$GG[i] =
        pedigreeData[match(pedigreeData$ringnr[i],
                           pedigreeData$father), "GG"]
    }
    if (is.na(pedigreeData$GG[i])) {
      pedigreeData$GG[i] =
        pedigreeData[match(pedigreeData$ringnr[i],
                           pedigreeData$mother), "GG"]
    }
  }
  
  # Three individuals have been missed. Give them the GG of their
  # fathers, as was supposed to happen earlier.
  pedigreeData[is.na(pedigreeData$GG), "GG"] =
    yearData[match(pedigreeData[is.na(pedigreeData$GG), "father"],
                   yearData$ringnr), 
             "islandIO"]
  
  # If mother or father is missing, replace by GG (phantom parents)
  pedigreeData[is.na(pedigreeData$motherId), "motherId"] = 
    pedigreeData[is.na(pedigreeData$motherId), "GG"]
  pedigreeData[is.na(pedigreeData$fatherId), "fatherId"] =
    pedigreeData[is.na(pedigreeData$fatherId), "GG"]
  
  # Need q_{ij} values for each individual, use ggcontrib() from
  # nadiv package.
  Q = ggcontrib(pedigreeData[, 4:6], 
                ggroups = c("inner", "outer", "other"))
  pedigreeData = cbind(pedigreeData, Q)
  
  morphData$inner =
    pedigreeData[match(morphData$id, pedigreeData$id), "inner"]
  morphData$outer =
    pedigreeData[match(morphData$id, pedigreeData$id), "outer"]
  morphData$other =
    pedigreeData[match(morphData$id, pedigreeData$id), "other"]
}

# # to find out how many animals had inner/outer/other natal island
# d.tmp = d.natal[d.natal$ringnr %in% d.morph$ringnr,]
# table(d.tmp$islandIO )

####################### Save Data #################################

if (center) {
  f = "Runs/morphData"
} else {
  f = "Runs/morphData_uncentered"
}

if (usePed) {
  f = paste0(f, "_pedigree_version")
}

save(morphData, file = paste0(f, ".RData"))

if (usePed) {
  save(A, Q, pedigreeMap, pedigreeData,
       file = "Runs/pedigree.RData")
}

save(natalData, file = "Runs/natalData_ped.RData")
save(yearData, file = "Runs/yearData_ped.RData")
