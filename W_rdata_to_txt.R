load("Data/GG_W_matrices.RData")
rm(W2)
# Add ID
W1 = cbind(data.frame(ID = dimnames(W1)[[1]]), as.data.frame(W1))
# Write table
write.table(W1, file = "Data/W1_beagle1.txt",
            quote = FALSE, row.names = FALSE)
rm(W1)

# W2
load("Data/GG_W_matrices.RData")
rm(W1)
# Add ID
W2 = cbind(data.frame(ID = dimnames(W2)[[1]]), as.data.frame(W2))
# Write table
write.table(W2, file = "Data/W2_beagle1.txt",
            quote = FALSE, row.names = FALSE)