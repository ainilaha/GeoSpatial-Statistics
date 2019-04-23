# list
coords<-c(2,3,9,6,5,2,7,9,5,3)
coords<-matrix(coords,nrow=5,ncol=2)
coords
dist(coords)
z<-c(3,4,2,4,6)
z<-matrix(z,nrow=5,ncol=1)
z
# install.packages(geoR)
library(geoR)
covmat<-varcov.spatial(coords = coords,cov.model = "spherical", nugget = 2.5,cov.pars = c(7.5,10))
covmat

# setup gamma
coords_all<-c(2,3,9,6,5,5,2,7,9,5,3,5)
coords_all<-matrix(coords_all,nrow=6,ncol=2)
coords_all
# setup covariance matrix with point for prediction
gamma1<-varcov.spatial(coords = coords_all,cov.model = "spherical", nugget = 2.5,cov.pars = c(7.5,10))
gamma1

# we are interested in last column
gamma<-gamma1$varcov[1:5,6]
gamma


# Now perform kriging using eq. 2.19 Book Hierarchical methods


# assuming a mean of 0
y_pred<-t(gamma)%*%solve(covmat$varcov)%*%z
y_pred

# using a mean

m<-c(3.8,3.8,3.8,3.8,3.8)
m<-matrix(m,nrow=5,ncol=1)
m

y_pred<-3.8+t(gamma)%*%solve(covmat$varcov)%*%(z-m)
y_pred




# now compute variance
var_pred<-10-t(gamma)%*%solve(covmat$varcov)%*%gamma
var_pred