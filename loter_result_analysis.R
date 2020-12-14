library(data.table)
# 
# load(file = "Runs/morphData_pedigree_version.RData")
# load(file = "Runs/pedigree.RData")
# 
loter =  fread(file = "Data/loter_result.txt")#"Data/loter_test.txt")
# 
# full_vcf = fread(file = "Beagle/Data_from_henrik/Helgeland_01_2018.vcf")
# 
# admixed_in = fread(file = "Data/Loter input 1/admixed_phased.vcf/admixed_phased.vcf")
# other_in = fread(file = "Data/Loter input 1/other_phased.vcf/other_phased.vcf")
# outer_in = fread(file = "Data/Loter input 1/outer_phased.vcf/outer_phased.vcf")
# inner_in = fread(file = "Data/Loter input 1/inner_phased.vcf/inner_phased.vcf")
# 
# outerRef = pedigreeData$ringnr[pedigreeData$outer == 1]
# innerRef = pedigreeData$ringnr[pedigreeData$inner == 1]
# otherRef = pedigreeData$ringnr[pedigreeData$other == 1]
# admixedRef = pedigreeData$ringnr[
#   !pedigreeData$ringnr %in% c(otherRef, outerRef, innerRef)]
# 
# unAdmixed = fread(file = "Beagle/Data_from_henrik/NOT_admixed_inds.txt")
# uninner = fread(file = "Beagle/Data_from_henrik/NOT_inner_inds.txt")
# unOuter = fread(file = "Beagle/Data_from_henrik/NOT_outer_inds.txt")
# unOther = fread(file = "Beagle/Data_from_henrik/NOT_other_inds.txt")

# # The right individuals, i.e. the admixed, were LA'd
# all_beagle_inds = colnames(full_vcf[, -(1:9)])
# (length(all_beagle_inds) - length(unAdmixed$V1)) == dim(loter)[1] / 2

#rownames(loter) = 

# library(lattice)
# 
# #Build the data
# nrowcol <- 1000
# dat <- matrix(ifelse(runif(nrowcol*nrowcol) > 0.5, 1, 0),
#               nrow=nrowcol)
# 
# #Build the palette and plot it
# pal = colorRampPalette(c("red", "yellow"), space = "rgb")
# levelplot(as.matrix(loter), main="", xlab="", ylab="", 
#           cuts=5, at=seq(-0.1, 2, 1))#, useRaster = TRUE)
# 
# install.packages("heatmap3")  
# library(heatmap3)
# heatmap3(loter[,], useRaster = TRUE)
  
  