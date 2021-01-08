library(BGData)

load.BGData("Data/loter/Run 3/AInner/AInner.RData")

proportion_different = numeric(3116)

for (row_num in 1798:3116) {
  x = as.numeric(AInner@geno[row_num, ])
  
  locus_LA_sum = unname(tapply(x, (seq_along(x) - 1) %/% 2, sum ))
  proportion_different[row_num] = sum(locus_LA_sum == 1) / length(locus_LA_sum)
  print(row_num)
}

