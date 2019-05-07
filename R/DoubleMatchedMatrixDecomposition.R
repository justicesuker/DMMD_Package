# This algorithm is still in development. 
# Do not use it until another update. Commented on May 6, 2019.

DoubleMatchedMatrixDecomposition <- function(X1, X2, angle_threshold = 90 * pi/180, variance1 = "unequal", variance2 = "unequal",  tol = .Machine$double.eps^0.5){
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)
  r1 = ProfileLikCluster(svd_x1$d, variance = variance1)$index
  r2 = ProfileLikCluster(svd_x2$d, variance = variance1)$index
  #Calculate JIVE in column space
  angle_result_c = angle_cal(X1, X2, tol = tol)
  principal_angle_c = angle_result_c$angle
  pv1_c = angle_result_c$principal_vector1
  pv2_c = angle_result_c$principal_vector2
  joint_rank_c = joint_angle_cluster(
    principal_angle_c, angle_threshold = angle_threshold, variance = variance2)$joint_rank
  joint_space_c = (pv1_c[,1:joint_rank_c] + pv2_c[,1:joint_rank_c])/2
  # J1 = projection(joint_space_c) %*% X1
  # J2 = projection(joint_space_c) %*% X2
  P_c = projection(joint_space_c)
  X1_per_c =  X1 - P_c %*% X1
  X2_per_c =  X2 - P_c %*% X2
  svd_x1_per_c = svd(X1_per_c)
  svd_x2_per_c = svd(X2_per_c)
  individual_rank_c = c(r1 - joint_rank_c, r2 - joint_rank_c)
  individual_space_c1 = svd_x1_per_c$u[,1:individual_rank_c[1]]
  individual_space_c2 = svd_x2_per_c$u[,1:individual_rank_c[2]]
  column_space_1 = cbind(joint_space_c,individual_space_c1)
  column_space_2 = cbind(joint_space_c,individual_space_c2)
  
  # Calculate JIVE in row space
  angle_result_r = angle_cal(t(X1), t(X2), tol = tol)
  principal_angle_r = angle_result_r$angle
  pv1_r = angle_result_r$principal_vector1
  pv2_r = angle_result_r$principal_vector2
  joint_rank_r = joint_angle_cluster(
    principal_angle_r, angle_threshold = angle_threshold, variance = variance2)$joint_rank
  joint_space_r = (pv1_r[,1:joint_rank_r] + pv2_r[,1:joint_rank_r])/2
  # J1 = projection(joint_space_c) %*% X1
  # J2 = projection(joint_space_c) %*% X2
  P_r = projection(joint_space_r)
  X1_per_r =  X1 - X1 %*% P_r
  X2_per_r =  X2 - X2 %*% P_r
  svd_x1_per_r = svd(X1_per_r)
  svd_x2_per_r = svd(X2_per_r)
  individual_rank_r = c(r1 - joint_rank_r, r2 - joint_rank_r)
  individual_space_r1 = svd_x1_per_r$v[,1:individual_rank_r[1]]
  individual_space_r2 = svd_x2_per_r$v[,1:individual_rank_r[2]]
  row_space_1 = cbind(joint_space_r,individual_space_r1)
  row_space_2 = cbind(joint_space_r,individual_space_r2) 
  signal_mat1 = projection(column_space_1) %*% X1 %*% projection(row_space_1)
  signal_mat2 = projection(column_space_2) %*% X2 %*% projection(row_space_2)
  result_c = col_jive(signal_mat1,signal_mat2)
  result_r = col_jive(t(signal_mat1,signal_mat2))
  return(list("column signals" = result_c, "row signals" = result_r))
}