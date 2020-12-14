library(data.table)

orig = fread(file = "Data/Helgeland_01_2018_Dosage.txt")
new2 = fread(file = "Data/plink2.raw")

x = which(colnames(orig) != colnames(new2))

rm(orig)
rm(new2)

ped = fread(file = "Data/Data_from_henrik/Helgeland_01_2018.ped", select = 1:50)

orig = fread(file = "Data/Helgeland_01_2018_Dosage.txt")
new2 = fread(file = "Data/plink2.raw")

origmap = fread(file = "Data/Helgeland_01_2018.map")
newmap = fread(file = "Data/Helgeland_01_2018new.map")
new3 = fread(file = "Data/plink3.raw")

