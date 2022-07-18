\name{DMMD-package}
\alias{DMMD-package}
\alias{DMMD}
\docType{package}
\title{
\packageTitle{DMMD}
}
\description{
\packageDescription{DMMD}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{DMMD}
\packageIndices{DMMD}
The goal of DMMD package is to extract joint and individual signals from double matched multi-view matrices. For example, if we have one matrice that contains gene data collected in normal tissue and one matrice collected in cancer tissue, it is helpful to use my package to find out genes that are common to both normal and cancer tissues and genes that are unique. Another example is about recommendation system. If we have two data sets that contain rankings of movies from various people in different time, by using the package, researchers can know signals (movies) that are not changed with time. There is a R package called ``r.jive" which has similar functionality, but is slow, and can only do one-way matched matrix analysis. 

In this package, the main function is SMMD, which can do single-matched matrix factorization efficiently. The rank selection is based on the method of profile likelihood. There is a function of \code{datagen} which can generate single-matched matrices that can be used to test functions. There are also some supporting functions that can double standardize a matrix, calculate principal angles, find joint structures and so on.

Although the name of my package is double-matched matrix decomposition, as of now, the final DMMD function is still under investgation.
}
\author{
\packageAuthor{DMMD}

Maintainer: \packageMaintainer{DMMD}
}
\references{
Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.

Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.

Zhu, Mu, and Ali Ghodsi. "Automatic dimensionality selection from the scree plot via the use of profile likelihood." Computational Statistics & Data Analysis 51.2 (2006): 918-930.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
...
}
\examples{
data = DoubleDataGen(n = 20, p = 16, rank = c(4, 3), rc = 2, rr = 1, nrep = 1)
X1 = data$X1[[1]]
X2 = data$X2[[1]]

svd_x1 = svd(X1)
svd_x2 = svd(X2)
r1 = ProfileLikCluster(svd_x1$d)$index
r2 = ProfileLikCluster(svd_x2$d)$index

result_DMMD = DMMD(X1,X2)
result_DMMD$rc
result_DMMDi = DMMD_i(X1,X2)
}
