load("Runs/pedigree.RData")
load("Runs/morphData_pedigree_version.RData")

############################## Decomposition of A #################

# # Only consider relevant individuals
# Acut = A[(sparrow.ped$ID %in% gg.data$ID), (sparrow.ped$ID %in% gg.data$ID)]

# Do Choleksy decomposition
T_chol = t(chol(A)) 
D = diag(diag(T_chol)) ^ 2
T_mat = T_chol %*% sqrt(solve(D))
TInv = solve(T_mat)

# Checks that matrices are correct
# sum(A - (T_mat %*% D %*% t(T_mat)))
# AInv = solve(A)
# sum(AInv - (t(TInv) %*% solve(D) %*% TInv))
# Is ok 

###########################################################

pedigreeData$scalingInnerEq10 = pedigreeData$inner
pedigreeData$scalingOuterEq10 = pedigreeData$outer
pedigreeData$scalingOtherEq10 = pedigreeData$other

pedigreeData[, c("scalingInnerEq10")] = 
  ifelse(!is.na(pedigreeData$motherIdOriginal) &
          is.na(pedigreeData$fatherIdOriginal) &
         !is.na(pedigreeData[
           match(pedigreeData$motherIdOriginal, pedigreeData$id),
           c("inner")]), # if only mother is known
  pedigreeData[
    match(pedigreeData$motherIdOriginal, pedigreeData$id),
    c("inner")], # then give it genetic group of mother
  ifelse(
    is.na(pedigreeData$motherIdOriginal) & 
   !is.na(pedigreeData$fatherIdOriginal) & 
   !is.na(pedigreeData[
     match(pedigreeData$fatherIdOriginal, pedigreeData$id), 
     c("inner")]), # if only father is known
    pedigreeData[
      match(pedigreeData$fatherIdOriginal, pedigreeData$id), 
      c("inner")], # then give it genetic group of father
    pedigreeData[ ,c("inner")])
)

pedigreeData[,c("scalingOuterEq10")] = 
  ifelse(!is.na(pedigreeData$motherIdOriginal) & 
          is.na(pedigreeData$fatherIdOriginal) & 
         !is.na(pedigreeData[
           match(pedigreeData$motherIdOriginal, pedigreeData$id), 
           c("outer")]), # if only mother is known
  pedigreeData[
    match(pedigreeData$motherIdOriginal, pedigreeData$id), 
    c("outer")],
  ifelse(
    is.na(pedigreeData$motherIdOriginal) &
   !is.na(pedigreeData$fatherIdOriginal) &
   !is.na(pedigreeData[
     match(pedigreeData$fatherIdOriginal, pedigreeData$id), 
     c("outer")]),
    pedigreeData[
      match(pedigreeData$fatherIdOriginal, pedigreeData$id), 
      c("outer")],
    pedigreeData[, c("outer")])
)

pedigreeData[,c("scalingOtherEq10")] =
  ifelse(!is.na(pedigreeData$motherIdOriginal) & 
          is.na(pedigreeData$fatherIdOriginal) & 
         !is.na(pedigreeData[
           match(pedigreeData$motherIdOriginal, pedigreeData$id), 
           c("other")]), # if only mother is known
  pedigreeData[
    match(pedigreeData$motherIdOriginal, pedigreeData$id), 
    c("other")],
  ifelse(is.na(pedigreeData$motherIdOriginal) & 
        !is.na(pedigreeData$fatherIdOriginal) & 
        !is.na(pedigreeData[
          match(pedigreeData$fatherIdOriginal, pedigreeData$id), 
          c("other")]),
    pedigreeData[
      match(pedigreeData$fatherIdOriginal, pedigreeData$id),
      c("other")],
    pedigreeData[, c("other")])
)

####################### Group Specific Ainv's #####################

# Prepare numerical IDs
# gg.ped$IDnum => pedigreeData$id
# gg.data$IDnum => morphData$id

# Prepare Inner-specific relatedness matrix

# Replace divisions by zero with large numbers
scalingInnerEq12 = ifelse(pedigreeData$inner > 0,
                       1 / (pedigreeData$inner) ^ 2, 1e12) 

# Eq. 10 in Muff et al 2019
d_inner = 1 - pedigreeData$scalingInnerEq10 * (1 - diag(D)) 

# As defined in eq. 12 in Muff et al 2019
D_tilde_inner = Diagonal(x = d_inner * 1 / scalingInnerEq12)
DInv_tilde_inner = Diagonal(x = 1 / (d_inner) * scalingInnerEq12) 

# Equation 13 in Muff et al 2019
A_inner = T_mat %*% D_tilde_inner %*% t(T_mat)
AInv_inner = t(TInv) %*% DInv_tilde_inner %*% TInv 

# Recast for reasons related to implementation of MCMCglmm
AInv_inner = as(AInv_inner, "dgCMatrix") 

# Give appropriate dimnames as required by MCMCglmm
AInv_inner@Dimnames = A_inner@Dimnames = A@Dimnames

# Repeat procedure for Vega specific relatedness matrix
scalingOuterEq12 = ifelse(pedigreeData$outer > 0, 
                       1 / (pedigreeData$outer) ^ 2, 1e12)
d_outer = 1 - pedigreeData$scalingOuterEq10 * (1- diag(D))
D_tilde_outer = Diagonal(x = d_outer * 1 / scalingOuterEq12)
DInv_tilde_outer = Diagonal(x = 1 / (d_outer) * scalingOuterEq12)
A_outer = T_mat %*% D_tilde_outer %*% t(T_mat)
AInv_outer = t(TInv) %*% DInv_tilde_outer %*% TInv
AInv_outer = as(AInv_outer, "dgCMatrix")
AInv_outer@Dimnames = A_outer@Dimnames = A@Dimnames

# Repeat procedure for "other" specific relatedness matrix
scalingOtherEq12 = ifelse(pedigreeData$other > 0,
                       1 / pedigreeData$other ^ 2, 1e12)  
d_other = 1 - pedigreeData$scalingOtherEq10 * (1 - diag(D))
D_tilde_other = Diagonal(x = d_other * 1 / scalingOtherEq12)
DInv_tilde_other = Diagonal(x = 1 / (d_other) * scalingOtherEq12)
A_other = T_mat %*% D_tilde_other %*% t(T_mat)
AInv_other = t(TInv) %*% DInv_tilde_other %*% TInv
AInv_other = as(AInv_other, "dgCMatrix")
AInv_other@Dimnames = A_other@Dimnames = A@Dimnames


###################################################################

save(A_inner,# AInv_inner,
     A_outer,# AInv_outer, 
     A_other,# AInv_other, 
     file = "Runs/A_hetped.RData")


