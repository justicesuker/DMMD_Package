# DMMD_Package
A R package that performs Double-Matched Matrix Decomposition (DMMD). 

Description 
-------

The goal of DMMD package is to extract joint and individual signals from double matched multi-view matrices. For example, if we have one matrice that contains gene data collected in normal tissue and one matrice collected in cancer tissue, it is helpful to use my package to find out genes that are common to both normal and cancer tissues and genes that are unique. Another example is about recommendation system. If we have two data sets that contain rankings of movies from various people in different time, by using the package, researchers can know signals (movies) that are not changed with time. There is a R package called *r.jive* which has similar functionality, but is slow, and can only do one-way matched matrix analysis. 

In this package, the main function is DMMD, which can do single-matched matrix factorization efficiently. The rank selection is based on the method of profile likelihood. There is a function of *DoubleDataGen* which can generate double-matched matrices. Another main function is DMMD_i which also updates joint structure. It is generally more accurate compared to DMMD but computationally less efficient.

References
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
