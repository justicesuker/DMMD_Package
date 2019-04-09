require(matlib)

#' Calculate exact column joint and individual structure of two matrices
#' 
#' Detailed description...What is this?
#' @param A First matrix
#' @param B Second matrix
#' @param tol Tolerence, default is the square root of machine precision.
#'
#' @return A list that contains joint and individual signals. 
#' \item{J1}{Joint signal of first matrix}
#' \item{J2}{Joint signal of second matrix}
#' \item{I1}{Individual signal of the first matrix}
#' \item{I2}{Individual signal of the second matrix}
#'          
#' @examples No example is a good example...

col_jive = function(A, B, tol = .Machine$double.eps^0.5)
{
#get basis for column spaces of A and B
result1 = svd(A, nu = nrow(A), nv = ncol(A))
result2 = svd(B, nu = nrow(B), nv = ncol(B))
U1 = result1$u
S1 = result1$d
U2 = result2$u
S2 = result2$d
m = 1
n = 1
index1 = 0
index2 = 0
while (m > tol) {
  index1 = index1 + 1
  if (index1 > length(S1))
    {break}
  m = S1[index1]
}
rank1 = index1 - 1
while (n > tol){
  index2 = index2 + 1
  if (index2 > length(S2))
    {break}
  n = S2[index2]
}
rank2 = index2 - 1
basismat1 = U1[,1:rank1]
basismat1 = cbind(basismat1)
basismat2 = U2[,1:rank2]

#Stack those basis together to be X and find the null space of X
X = cbind(basismat1,basismat2)
result = svd(X, nu = nrow(X), nv = ncol(X))
S = result$d
V = result$v
num = 1
index = 0
while (num > tol){
  index = index + 1
  if (index > length(S))
    {break}
  num = S[index]
}
p = dim(V)[2]

# First consider there is no non-zero intersection of column space
if (index == p + 1){
  J1 = matrix(rep(0),nrow = nrow(A), ncol = ncol(A))
  J2 = matrix(rep(0),nrow = nrow(A), ncol = ncol(A))
  I1 = A
  I2 = B
  newList <- list("J1" = J1, "J2" = J2, "I1" = I1, "I2" = I2)
}

# Non-trivial cases
else {
  N = V[,index:p]
  N = cbind(N)
  #find the basis of the intersection of column space
  tempmat = N[1:rank1,]
  if (rank1 == 1){
    tempmat = rbind(tempmat)
  }
  J = basismat1%*%tempmat
  #Construct Joint and Individual structure
  G = t(J)%*%J
  if (dim(G)[1] == 1 && dim(G)[2] == 1){
    P = J%*%(1/G)%*%t(J)
  }
  else{
    P = J%*%inv(G)%*%t(J)
  }
  J1 = P%*%A
  J2 = P%*%B
  I1 = A - J1
  I2 = B - J2
  newList <- list("J1" = J1, "J2" = J2, "I1" = I1, "I2" = I2)
}
return(newList)
}