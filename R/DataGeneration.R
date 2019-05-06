#' Function to generate simulation data
#' @description Function that generates a dataset using provided true joint rank and individual ranks. Number of matrices in one dataset is restricted to two.
#' @param n The number of rows of matrices that will be generated
#' @param p The number of columns of matrices that will be generated
#' @param joint_rank The true joint rank of matrices to be generated 
#' @param individual_rank A vector of length two. The true individual ranks of matrices that will be generated 
#' @param nrep Number of datasets to be generated. A pair of matrices is called a dataset.
#' @param std The standard deviation of normal distribution used to generatd noise matrices. Default is 0.1
#' 
#' @return A list that contains the following:
#' \item{X1_list}{A list of first matrices generated.}
#' \item{X2_list}{A list of second matrices generated}
#' \item{trueJ1_list}{A list of true joint signal matrices of the \code{X1_list}}
#' \item{trueJ2_list}{A list of true joint signal matrices of the \code{X2_list}}
#' \item{trueI1_list}{A list of true individual signal matrices of the \code{X1_list}}
#' \item{trueI2_list}{A list of true individual signal matrices of the \code{X2_list}}
#' 
#' @examples
#' datagen(n = 100, p = 20, joint_rank = 5, individual_rank = c(5,7), nrep = 20)
#' 
datagen <- function(n = 1000, p = 300, joint_rank = 30, individual_rank = c(100,120), nrep = 20, std = 0.1){
  if (n %% 4 != 0){stop("4 must be a factor of n.")}
  if (max(joint_rank + individual_rank) > min(n, p)){stop("Total rank exceeds the size of matrices.")}
  X1_list = list()
  X2_list = list()
  trueJ1_list = list()
  trueJ2_list = list()
  trueI1_list = list()
  trueI2_list = list()
  for (i in 1:nrep){
    joint_index = sample(1:(n/2),joint_rank)
    individual_index1 = sample(1:(n/4),individual_rank[1])
    individual_index2 = sample(1:(n/4),individual_rank[2])
    tempv = rep(0,n)
    tempv[joint_index[1]] = 1
    joint_matrix = tempv
    tempv21 = rep(0,n)
    tempv21[individual_index1[1]+(n/2)] = 1
    individual_matrix1 = tempv21
    tempv22 = rep(0,n)
    tempv22[individual_index2[1]+ (3/4*n)] = 1
    individual_matrix2 = tempv22
    for (i in 2:joint_rank){
      tempv = rep(0,n)
      tempv[joint_index[i]] = 1
      joint_matrix = cbind(joint_matrix,tempv)
    }
    for (i in 2:individual_rank[1]){
      tempv21 = rep(0,n)
      tempv21[individual_index1[i]+(n/2)] = 1
      individual_matrix1 = cbind(individual_matrix1,tempv21)
    }
    for (i in 2:individual_rank[2]){
      tempv22 = rep(0,n)
      tempv22[individual_index2[i]+(3/4*n)] = 1
      individual_matrix2 = cbind(individual_matrix2,tempv22)
    }
    temp_loading_joint_1 = matrix(rt(joint_rank * p ,df = 5), 
                                  nrow = joint_rank, ncol = p)
    temp_loading_joint_2 = matrix(rt(joint_rank * p ,df = 10), 
                                  nrow = joint_rank, ncol = p)
    temp_loading_ind_1 = matrix(rnorm(individual_rank[1] * p,sd = 2),
                                nrow = individual_rank[1], ncol = p)
    temp_loading_ind_2 = matrix(rnorm(individual_rank[2] * p,sd = 2.5),
                                nrow = individual_rank[2], ncol = p)
    
    Joint_1 = joint_matrix %*% temp_loading_joint_1
    Joint_2 = joint_matrix %*% temp_loading_joint_2
    Ind_1 = individual_matrix1 %*% temp_loading_ind_1
    Ind_2 = individual_matrix2 %*% temp_loading_ind_2
    
    E1 = matrix(rnorm(n*p, mean = 0, sd = std), nrow = n)
    E2 = matrix(rnorm(n*p, mean = 0, sd = std), nrow = n)
    X1 = Joint_1 + Ind_1 + E1
    X2 = Joint_2 + Ind_2 + E2
    
    X1_list = append(X1_list, list(X1))
    X2_list = append(X2_list, list(X2))
    trueJ1_list = append(trueJ1_list, list(Joint_1))
    trueJ2_list = append(trueJ2_list, list(Joint_2))
    trueI1_list = append(trueI1_list, list(Ind_1))
    trueI2_list = append(trueI2_list, list(Ind_2))
  }
  return(list("X1_list" = X1_list, "X2_list" = X2_list, "trueJ1_list" = trueJ1_list, "trueJ2_list" = trueJ2_list,
         "trueI1_list" = trueI1_list, "trueI2_list" = trueI2_list))
}