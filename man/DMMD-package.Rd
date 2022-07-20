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
The goal of DMMD package is to extract joint and individual signals from double-matched multi-view matrices. The main function is DMMD_Fit, which can do double-matched matrix factorization efficiently. DMMD_i is the updated version of DMMD_Fit which also updates the joint structure but is less efficient. The rank selection is based on the method of profile likelihood or edge distribution. There is a function of \code{DoubleDataGen} which can generate double-matched matrices that can be used to test functions.
}

\author{
Dongbang Yuan and Irina Gaynanova \cr
Maintainer: \packageMaintainer{DMMD}
}

\keyword{data integration}
\keyword{dimension reduction}
\keyword{matrix factorization}
\keyword{singular value decomposition}

\references{
Dongbang Yuan & Irina Gaynanova (2022) Double-Matched Matrix Decomposition for Multi-View Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2022.2067860
}
\seealso{
Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.

Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.

Zhu, Mu, and Ali Ghodsi. "Automatic dimensionality selection from the scree plot via the use of profile likelihood." Computational Statistics & Data Analysis 51.2 (2006): 918-930.
}
\examples{
data = DoubleDataGen(n = 20, p = 16, rank = c(4, 3), rc = 2, rr = 1, nrep = 1)
X1 = data$X1[[1]]
X2 = data$X2[[1]]

svd_x1 = svd(X1)
svd_x2 = svd(X2)
r1 = ProfileLikCluster(svd_x1$d)$index
r2 = ProfileLikCluster(svd_x2$d)$index

result_DMMD = DMMD_Fit(X1,X2)
result_DMMD$rc
result_DMMDi = DMMD_i(X1,X2)
}