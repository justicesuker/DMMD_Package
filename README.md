# DMMD_Package
R package that does decomposition on double-matched multi-view matrices

##Description 
The goal of DMMD package is to extract joint and individual signals from double matched multi-view matrices. For example, if we have one matrice that contains gene data collected in normal tissue and one matrice collected in cancer tissue, it is helpful to use my package to find out genes that are common to both normal and cancer tissues and genes that are unique. Another example is about recommendation system. If we have two data sets that contain rankings of movies from various people in different time, by using the package, researchers can know signals (movies) that are not changed with time. There is a R package called ``r.jive" which has similar functionality, but is slow, and can only do one-way matched matrix analysis. 
In this package, the main function is SMMD, which can do single-matched matrix factorization efficiently. The rank selection is based on the method of profile likelihood. There is a function of \code{datagen} which can generate single-matched matrices that can be used to test functions. There are also some supporting functions that can double standardize a matrix, calculate principal angles, find joint structures and so on.
Although the name of my package is double-matched matrix decomposition, as of now, the final DMMD function is still under investgation.

##Installation
`devtools::install_github("justicesuker/DMMD_package")`

##Examples
```{r}
X = matrix(c(1,1,1,1,1,0),nrow = 3, ncol = 2)
Y = matrix(c(1,1,1,2,1,0),nrow = 3, ncol = 2)
angle_cal(X,Y)
x = c(20,9.5,9,8,6,5,4.5,3)
ProfileLikCluster(x)

data = datagen(n = 100, p = 20, joint_rank = 5, individual_rank = c(5,7), nrep = 1)
X1 = data$X1_list[[1]]
X2 = data$X2_list[[1]]
angle_vec = angle_cal(X1,X2)$angle
joint_angle_cluster(angle_vec)
result_SMMD = SMMD(X1,X2)
result_SMMD$joint_rank
result_SMMD$individual_rank
```
