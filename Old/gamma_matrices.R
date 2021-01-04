library(BGData)

load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/AInner/AInner.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/AOuter/AOuter.RData"))

load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/AOther/AOther.RData"))

cores = ifelse(Sys.info()[['sysname']] == "Windows", 1, 8)

thetaInner = getG(AInner@geno,
             center = FALSE, scale = FALSE, scaleG = FALSE,
             nCores = cores, verbose = TRUE)
thetaInner = thetaInner / dim(geno(AInner))[2]
gammaInner = innerGRM / thetaInner
save(gammaInner, file = "Runs/gammaInner.Rdata")

thetaOuter = getG(AOuter@geno,
                  center = FALSE, scale = FALSE, scaleG = FALSE,
                  nCores = cores, verbose = TRUE)
thetaOuter = thetaOuter / dim(geno(AOuter))[2]
gammaOuter = outerGRM / thetaOuter
save(gammaOuter, file = "Runs/gammaOuter.Rdata")

thetaOther = getG(AOther@geno,
                  center = FALSE, scale = FALSE, scaleG = FALSE,
                  nCores = cores, verbose = TRUE)
thetaOther = thetaOther / dim(geno(AOther))[2]
gammaOther = otherGRM / thetaOther
save(gammaOther, file = "Runs/gammaMOther.Rdata")