#' Function that uses edge distribution (ED) method to calculate the rank i.e., the number of factors.
#' @description This function is a direct translation from Python function 'select_K_by_ED' in dcca.py.
#' @references Onatski, Alexei. “DETERMINING THE NUMBER OF FACTORS FROM EMPIRICAL DISTRIBUTION OF EIGENVALUES.” The Review of Economics and Statistics 92, no. 4 (2010): 1004–16. http://www.jstor.org/stable/40985808.
#' @importFrom stats lm
#' @param x The input vector for selecting the rank
#' @param k_max The maximum rank to be considered, restricted to less than or equal to half of the length of the vector. If not given, one tenth of the length of the vector is considered in general
#' @param k_min The minimum rank to be considered. If not given, it is treated as 1
#' @param maxiter The maximum iterations for the algorithm. Default 100
#'
#' @return The rank estimated by ED method
#' @export
#'
#' @examples
#' x = c(20,9.5,9,8,6,5,4.5,3)
#' Select_ED_Rank(x)
Select_ED_Rank <- function(x, k_max = NULL, k_min = NULL, maxiter = 100){
  m = length(x)
  if (is.null(k_min)){
    k_min = 1
  }
  if (is.null(k_max)){
    k_max_star = sum(x >= mean(x))
    k_max = ceiling(min(k_max_star, 0.1*m))
  }
  if (k_max > 0.5*m){
    stop('k_max cannot be greater than half of the length of the vector, please try another number.')
  }
  if (k_max < k_min){
    stop('k_min is greater than K_max, please try another number.')
  }

  j = k_max + 1
  x_diff = x[1:k_max] - x[2:j]
  k_pre = -1
  for (t in 1:maxiter){
    y = x[j:(j+4)]
    j_vec = k_max + seq(0,4,1)
    temp_x = (j_vec)^(2/3)
    coeff = lm(y~temp_x)$coefficients["temp_x"]
    names(coeff) = NULL
    delta = 2 * abs(coeff)
    index = which(x_diff >= delta)
    if (length(index) == 0){
      k = 0
    }
    else{
      k = max(index)
    }
    if (k_pre == k){
      break
    }
    k_pre = k
    j = k + 1
  }
  return(max(k,k_min))
}