# DMMD_Package
R package that performs Double-Matched Matrix Decomposition (DMMD). 

Description
-------

The goal of DMMD package is to extract joint and individual signals from double-matched multi-view matrices. The main function is DMMD_Fit, which can do double-matched matrix factorization efficiently. DMMD_i is the updated version of DMMD_Fit which also updates the joint structure but is less efficient. The rank selection is based on the method of profile likelihood or edge distribution. There is a function of \code{DoubleDataGen} which can generate double-matched matrices that can be used to test functions.

Reference
-------
[Dongbang Yuan & Irina Gaynanova (2022) Double-Matched Matrix Decomposition for Multi-View Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2022.2067860](https://www.tandfonline.com/doi/abs/10.1080/10618600.2022.2067860)

See also
-------
[Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.](https://arxiv.org/pdf/1704.02060.pdf).

[Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.](https://arxiv.org/pdf/1102.4110.pdf).

[Zhu, Mu, and Ali Ghodsi. "Automatic dimensionality selection from the scree plot via the use of profile likelihood." Computational Statistics & Data Analysis 51.2 (2006): 918-930.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf).

Installation
-------
`devtools::install_github("justicesuker/DMMD_package")`

Example
-------
```{r}
X = matrix(c(1,1,1,1,1,0),nrow = 3, ncol = 2)
Y = matrix(c(1,1,1,2,1,0),nrow = 3, ncol = 2)
angle_cal(X,Y)
x = c(20,9.5,9,8,6,5,4.5,3)
ProfileLikCluster(x)

data = DoubleDataGen(n = 20, p = 16, rank = c(4, 3), rc = 2, rr = 1, nrep = 1)
X1 = data$X1[[1]]
X2 = data$X2[[1]]

svd_x1 = svd(X1)
svd_x2 = svd(X2)
ProfileLikCluster(svd_x1$d)$index
ProfileLikCluster(svd_x2$d)$index
r1 = 4
r2 = 3

# Get the estimated column/row space of X1 and X2
X1_est_c = as.matrix(svd_x1$u[,1:r1])
X2_est_c = as.matrix(svd_x2$u[,1:r2])
X1_est_r = as.matrix(svd_x1$v[,1:r1])
X2_est_r = as.matrix(svd_x2$v[,1:r2])
  
# Get the principal angles
angle_result_c = angle_cal(X1_est_c, X2_est_c, tol = tol)
principal_angle_c = angle_result_c$angle
joint_rank_c = joint_angle_cluster(principal_angle_c)$joint_rank

result_DMMD = DMMD_Fit(X1,X2)
result_DMMD$rc
result_DMMDi = DMMD_i(X1,X2, verbose = TRUE)
```
