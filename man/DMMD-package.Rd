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

}
\author{
\packageAuthor{DMMD}

Maintainer: \packageMaintainer{DMMD}
}
\references{
Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.

Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
r.jive: Perform JIVE Decomposition for Multi-Source Data. It is created by Michael J. O'Connell and Eric F. Lock, which is also available in R.
\code{\link[<r.jive>:<r.jive>-package]{https://cran.r-project.org/web/packages/r.jive/}}
}
\examples{
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
}
