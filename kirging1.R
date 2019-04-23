# list
coords<-c(2,3,9,6,5,2,7,9,5,3)
coords<-matrix(coords,nrow=5,ncol=2)
coords
dist(coords)
z<-c(3,4,2,4,6)
z<-matrix(z,nrow=5,ncol=1)
z
t(coords)%*%z

data1<-list(xy=coords,dat=z)
data1$dat
data1$xy
library(geoR)
library(akima)
int.data1<-interp(data1$xy[,1],data1$xy[,2],data1$dat,extrap=TRUE)
int.data1

# Graphs
image(int.data1)
contour(int.data1, add=T)
persp(int.data1)
persp(int.data1,theta = -70, phi = 45, d=1)

# variograms
var<-variog(coords = data1$xy, data = data1$dat,max.dist=15,estimator.type="classical")
plot(var)

# kriging
loci<-expand.grid(seq(0,10,l=11),seq(0,10,l=11))
kc<-krige.conv(coords = data1$xy, data = data1$dat,loc=loci,krige=krige.control(cov.pars=c(1,.25)))
image(kc)
contour(kc,add=T)
persp(kc)
persp(kc,theta = -70, phi = 45, d=1)

point<-krige.conv(coords = data1$xy, data = data1$dat,loc=c(5,1),
                  krige=krige.control(cov.pars=c(7.5,10),cov.model="spherical",nugget=2.5))
point

dist<-dist(data1$xy)
dist
