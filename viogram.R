library(geoR)
library(tcltk)
data(s100)
s100
s100$coords[,1]
plot.geodata(s100)
# install.packages('scatterplot3d')
plot.geodata(s100,scatter3d=TRUE)
library(akima)
int.s100<-interp.new(s100$coords[,1],s100$coords[,2],s100$data,extrap=TRUE)
int.s100

# Graphs
image(int.s100)
contour(int.s100, add=T)
persp(int.s100)
persp(int.s100,theta = -70, phi = 45, d=1)
persp(int.s100,theta = -70, phi = 45, d=1,col="blue",shade=0.75,border=NA)

# variograms
var4<-variog4(s100,max.dist=1)
plot(var4)
var<-variog(s100,max.dist=1,estimator.type="classical")
plot(var)
var<-variog(s100,max.dist=1,estimator.type="modulus")
plot(var)
eyefit(var,silent=TRUE)
varf<-variofit(var,ini.cov.pars=c(0.91,0.59),cov.model="exponential",fix.nugget=FALSE, nugget=0.23)
varf
varlik<-likfit(s100,ini.cov.pars=c(0.91,0.59),cov.model="exponential",fix.nugget=FALSE,nugget=0.23, method.lik = "ML")
summary(varlik)

# Classical Kriging
loci<-expand.grid(seq(0,1,l=11),seq(0,1,l=11))
kc<-krige.conv(s100,loc=loci,krige=krige.control(cov.pars=c(1,.25)))
image(kc)
contour(kc,add=T)