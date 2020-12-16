library(data.table)
library(BGData)
source("My_R_code/file_backed_mat.R")

loter_run = 4

innerFile = file(paste0("Data/loter/Run ", loter_run,
                        "/Loter input ", loter_run,
                        "/inner_phased", loter_run, ".vcf.gz"), 
                 open = "r")
innerLine  = scan(innerFile, nlines = 1, skip = 8,
                  what = character(), quiet = TRUE)
close(innerFile)
innerInds = innerLine[-(1:9)]

outerFile = file(paste0("Data/loter/Run ", loter_run,
                        "/Loter input ", loter_run,
                        "/outer_phased", loter_run, ".vcf.gz"), 
                 open = "r")
outerLine  = scan(outerFile, nlines = 1, skip = 8,
                  what = character(), quiet = TRUE)
close(outerFile)
outerInds = outerLine[-(1:9)]

otherFile = file(paste0("Data/loter/Run ", loter_run,
                        "/Loter input ", loter_run,
                        "/other_phased", loter_run, ".vcf.gz"), 
                 open = "r")
otherLine  = scan(otherFile, nlines = 1, skip = 8,
                  what = character(), quiet = TRUE)
close(otherFile)
otherInds = otherLine[-(1:9)]

admixedFile = file(paste0("Data/loter/Run ", loter_run,
                          "/Loter input ", loter_run,
                          "/admixed_phased", loter_run, ".vcf.gz"), 
                   open = "r")
admixedLine  = scan(admixedFile, nlines = 1, skip = 8,
                  what = character(), quiet = TRUE)
close(admixedFile)
admixedInds = admixedLine[-(1:9)]
SNPs = fread(file = paste0("Data/loter/Run ", loter_run,
                           "/Loter input ", loter_run,
                           "/admixed_phased", loter_run,".vcf.gz"),
              select = 3)$ID


inds = c(innerInds, outerInds, otherInds, admixedInds)

AOuter = BGData(
  geno = initFileBackedMatrix(length(inds), 2 * length(SNPs),
                              folderOut = paste0("Data/loter/Run ",
                                                 loter_run,
                                                 "/AOuter"),
                              outputType = "boolean"),
  pheno = data.frame(ID = inds))
AInner = BGData(
  geno = initFileBackedMatrix(length(inds), 2 * length(SNPs),
                              folderOut = paste0("Data/loter/Run ",
                                                 loter_run,
                                                 "/AInner"),
                              outputType = "boolean"),
  pheno = data.frame(ID = inds))
AOther = BGData(
  geno = initFileBackedMatrix(length(inds), 2 * length(SNPs),
                              folderOut = paste0("Data/loter/Run ",
                                                 loter_run,
                                                 "/AOther"),
                              outputType = "boolean"),
  pheno = data.frame(ID = inds))

i = c(1, length(innerInds))
for (row in i[1]:i[2]) {
  AOuter@geno[row,] = AOther@geno[row,] = rep(0, 2 * length(SNPs))
  AInner@geno[row,] = rep(1, 2 * length(SNPs))
  print(row)
}

i = i + c(length(innerInds), length(outerInds))

for (row in i[1]:i[2]) {
  AInner@geno[row, ] = AOther@geno[row,] = rep(0, 2 * length(SNPs))
  AOuter@geno[row, ] = rep(1, 2 * length(SNPs))
  print(row)
}

i = i + c(length(outerInds), length(otherInds))

for (row in i[1]:i[2]) {
  AInner@geno[row, ] = AOuter@geno[row,] = rep(0, 2 * length(SNPs))
  AOther@geno[row, ] = rep(1, 2 * length(SNPs))
  print(row)
}
i = i + c(length(otherInds), length(admixedInds))
# Order in loter: Outer == 0, Inner == 1, Other == 2
SNPs1 = seq(1, by = 2, to = 2 * length(SNPs))
SNPs2 = seq(2, by = 2, to = 2 * length(SNPs))
admixedFile = file(paste0("Data/loter/Run ", loter_run,
                          "/loter_result", loter_run, ".txt"), 
                   open = "r")
for (row in i[1]:i[2]) {
  line1  = scan(admixedFile, nlines = 1, 
                what = integer(), quiet = TRUE)
  AOuter@geno[row, SNPs1] = ifelse(line1 == 0, 1, 0)
  AInner@geno[row, SNPs1] = ifelse(line1 == 1, 1, 0)
  AOther@geno[row, SNPs1] = ifelse(line1 == 2, 1, 0)
  line2 = scan(admixedFile, nlines = 1, 
               what = integer(), quiet = TRUE)
  AOuter@geno[row, SNPs2] = ifelse(line2 == 0, 1, 0)
  AInner@geno[row, SNPs2] = ifelse(line2 == 1, 1, 0)
  AOther@geno[row, SNPs2] = ifelse(line2 == 2, 1, 0)
  print(row)
}
close(admixedFile)

load.BGData(file = paste0("Data/loter/Run ", loter_run,
                          "/W/W.RData"))

dimnames(AInner@geno) = dimnames(W@geno)
dimnames(AOuter@geno) = dimnames(W@geno)
dimnames(AOther@geno) = dimnames(W@geno)


save(AInner, file = paste0("Data/loter/Run ", loter_run,
                           "/AInner/AInner.RData"))
save(AOuter, file = paste0("Data/loter/Run ", loter_run,
                           "/AOuter/AOuter.RData"))
save(AOther, file = paste0("Data/loter/Run ", loter_run,
                           "/AOther/AOther.RData"))


# write.table(
#   cbind(
#     data.table(ID = c(innerInds, outerInds, otherInds, admixedInds)),
#     rbind(matrix(0, nrow = length(innerInds), ncol = dim(admixedLA)[2]),
#           matrix(1, nrow = length(outerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(otherInds), ncol = dim(admixedLA)[2]),
#           ifelse(loterAllele1 == 0, 1, 0))),
#   file = "Data/loter/Run 2/LA_allele1Outer.txt", quote = FALSE,
#   row.names = FALSE)
# print("Writing inner1")
# write.table(
#   cbind(
#     data.table(ID = c(innerInds, outerInds, otherInds, admixedInds)),
#     rbind(matrix(1, nrow = length(innerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(outerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(otherInds), ncol = dim(admixedLA)[2]),
#           ifelse(loterAllele1 == 1, 1, 0))),
#   file = "Data/loter/Run 2/LA_allele1Inner.txt", quote = FALSE,
#   row.names = FALSE)
# print("Writing other1")
# write.table(
#   cbind(
#     data.table(ID = c(innerInds, outerInds, otherInds, admixedInds)),
#     rbind(matrix(0, nrow = length(innerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(outerInds), ncol = dim(admixedLA)[2]),
#           matrix(1, nrow = length(otherInds), ncol = dim(admixedLA)[2]),
#           ifelse(loterAllele1 == 2, 1, 0))),
#   file = "Data/loter/Run 2/LA_allele1Other.txt", quote = FALSE,
#   row.names = FALSE)
# 
# rm(loterAllele1)
# loterAllele2 = admixedLA[(1:(dim(admixedLA)[1] / 2)) * 2, ]
# 
# print("Writing outer2")
# write.table(
#   cbind(
#     data.table(ID = c(innerInds, outerInds, otherInds, admixedInds)),
#     rbind(matrix(0, nrow = length(innerInds), ncol = dim(admixedLA)[2]),
#           matrix(1, nrow = length(outerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(otherInds), ncol = dim(admixedLA)[2]),
#           ifelse(loterAllele2 == 0, 1, 0))),
#   file = "Data/loter/Run 2/LA_allele2Outer.txt", quote = FALSE,
#   row.names = FALSE)
# print("Writing inner2")
# write.table(
#   cbind(
#     data.table(ID = c(innerInds, outerInds, otherInds, admixedInds)),
#     rbind(matrix(1, nrow = length(innerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(outerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(otherInds), ncol = dim(admixedLA)[2]),
#           ifelse(loterAllele2 == 1, 1, 0))),
#   file = "Data/loter/Run 2/LA_allele2Inner.txt", quote = FALSE,
#   row.names = FALSE)
# print("Writing other2")
# write.table(
#   cbind(
#     data.table(ID = c(innerInds, outerInds, otherInds, admixedInds)),
#     rbind(matrix(0, nrow = length(innerInds), ncol = dim(admixedLA)[2]),
#           matrix(0, nrow = length(outerInds), ncol = dim(admixedLA)[2]),
#           matrix(1, nrow = length(otherInds), ncol = dim(admixedLA)[2]),
#           ifelse(loterAllele2 == 2, 1, 0))),
#   file = "Data/loter/Run 2/LA_allele2Other.txt", quote = FALSE,
#   row.names = FALSE)
