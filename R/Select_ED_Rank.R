# Function that uses edge distribution (ED) method to calculate the rank i.e., the number of factors
# Refer to Select_ED_Rank_K in D-CCA.py in EDRank_FromDCCA folder.
# Onatski, Alexei. "Determining the number of factors from empirical distribution of eigenvalues." 
# The Review of Economics and Statistics 92, no. 4 (2010): 1004-1016.
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