# GENERATE CORRELATED DATA AND INFERENCE FOR VARIOGRAM PARAMETERS
# This procedure computes least squares mean
# Load needed libraries
library(geoR)
library(sp)

# set.seed(2401)
# set.seed(5137)
# set up data matrix size nxn
# R is number of simulations
# count is the number of CI that capture the mean
# TOL is tolerance for minimum eigenvalue
TOL<-1E-8
R<-1
mu <- 5
n<-25
graphs <-TRUE

count<-0
for( r in 1:R){
  X<-matrix(mu,nrow=n,ncol=n, byrow= TRUE)
  # Gaussian semivariogram
  # Gaussian variogram parameters
  sigma2<-1
  theta1<-0
  theta2<-1
  theta3<-3
  
  # Coordinates
  unit<-1
  n2<-n*n
  Coord<-matrix(0,nrow=n2, ncol = 2)
  for(i in 1:n) {
    for(j in 1:n){
      Coord[j+(i-1)*n,1] <-i*unit
      Coord[j+(i-1)*n,2] <-j*unit
    }
  }
  
  # distance matrix
  D<-as.matrix(dist(Coord))
  
  
  # Covariance Matrix
  SGM<-matrix(0,nrow=n2,ncol=n2, byrow= TRUE)
  for(i in 1:n2) {
    for(j in 1:n2){
      if (D[i,j]==0) SGM[i,j] <- sigma2 else SGM[i,j] <- sigma2-(theta1+theta2*(1-exp(-(D[i,j]/theta3)^2)))
    }
  }
  ev<-eigen(SGM)
  
  # Substitute small eigenvalues for TOL
  for(i in 2:n2){
    if( ev$values[i] < TOL) ev$values[i] <- TOL 
  }
  SGM<- ev$vectors%*%diag(ev$values)%*%t(ev$vectors)
  
  
  # image graph for variance matrix
  if (graphs){
    x <- seq(1, nrow(SGM), length.out = nrow(SGM))
    y <- seq(1, ncol(SGM), length.out = ncol(SGM))
    # rainbow(n, s = 1, v = 1, start = 0, end = max(1,999)/1000, gamma = 1, alpha = 1)
    image(x,y,SGM, col=heat.colors(10000))
    # image(x,y,SGM, col=gray((0:1000)/1000))
  }
  
  # Choleski decomposition
  Lu<-chol(SGM)
  L<-t(Lu)
  e<-rnorm(n2, mean = 0, sd = 1)
  MU<-mu*c(rep(1,n2))
  Z<-MU+L%*%e
  Z<-matrix(Z,nrow = n ,ncol = n)
  
  # Graphs
  if (graphs){
    x <- seq(0, nrow(Z), length.out = nrow(Z))
    y <- seq(0, ncol(Z), length.out = ncol(Z))
    dev.new()
    image(x,y,Z)
    contour(x,y,Z, add=T)
    dev.new()
    persp(x,y,Z,theta = -70, phi = 45,d=1,shade=0.75,scale=FALSE,ticktype = "detailed", col ="green")
  }
  
  # Arithmetic mean
  Zbar<-mean(Z)
  
  # Estimation of the mean using least squares
  Zv<-as.vector(Z)
  unit<-c(rep(1,n2))
  SGMINV<-chol2inv(Lu,size=n2)
  Mlsq<-solve(t(unit)%*%SGMINV%*%unit)%*%(t(unit)%*%SGMINV%*%Zv)
  
  
  ###################################
  # PARAMETRIC BOOTSTRAP
  # Inference for variogram parameters
  # Compute residuals
  ###################################
  
  Muh<-Mlsq[1,1]
  # Muh<-Zbar
  Zd<-Z-Muh
  Coords<-c(Coord,as.vector(Zd))
  Coords<-matrix(Coords, nrow=n2, ncol=3)
  Coords<-as.geodata(Coords)
  
  vario1<-variog(Coords, max.dist=30 )
  plot(vario1)
  Ini<-matrix(c(.5,.8,1.2,1.8,2.5,3,1.5,2,2.5,3,3.5,4),nrow=6,ncol=2)
  g1 <- variofit(vario1, ini.cov.pars = Ini, fix.nugget = TRUE, nugget = 0, max.dist = 10, cov.model = "gaussian")
  #ml<-likfit(vario1, ini=c(1,1), fix.nug=TRUE, lik.met="REML")
  
  # COMPUTATION OF CONFIDENCE INTERVALS
  # Gaussian semivariogram
  # Gaussian variogram parameters obtained from previous step
  sigma2<-g1$cov.pars[1]
  theta1<-0
  theta2<-g1$cov.pars[1]
  theta3<-g1$cov.pars[2]
  
  # New Covariance Matrix computed from fitted variogram
  SGM<-matrix(0,nrow=n2,ncol=n2, byrow= TRUE)
  for(i in 1:n2) {
    for(j in 1:n2){
      if (D[i,j]==0) SGM[i,j] <- sigma2 else SGM[i,j] <- sigma2-(theta1+theta2*(1-exp(-(D[i,j]/theta3)^2)))
    }
  }
  ev<-eigen(SGM)
  
  # Substitute small eigenvalues for TOL
  for(i in 2:n2){
    if( ev$values[i] < TOL) ev$values[i] <- TOL 
  }
  SGM<- ev$vectors%*%diag(ev$values)%*%t(ev$vectors)
  
}
theta1
theta2
theta3
sigma2
