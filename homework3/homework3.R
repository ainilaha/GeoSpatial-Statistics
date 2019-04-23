url <- "https://www.counterpointstat.com/uploads/1/1/9/3/119383887/myscallops.txt" 
myscallops <- read.table(url,header = T)
coords <- as.matrix(myscallops[,c("lat","long")]) 
lgcatch<- myscallops$lgcatch

library(geoR)
library(spBayes)
bins = 50
max.dist <- 0.5*max(iDist(coords))

myscps.vario <- variog(coords = coords, data = lgcatch, uvec = (seq(0, max.dist, length = bins)))
plot(myscps.vario)
eyefit(myscps.vario,silent=TRUE)


fit.lgcatch<- variofit(myscps.vario,cov.model="exponential",fix.nugget=FALSE, nugget=0.7)
fit.lgcatch


point<-krige.conv(coords = coords, data = lgcatch,loc=c(length(lgcatch),1),
                  krige=krige.control(cov.pars=c(4.5,0.5),cov.model="exponential",nugget=0.63))
point

pred_low <-point$predict - 2*sqrt(point$krige.var)
pred_high <-point$predict + 2*sqrt(point$krige.var)
print(paste("The PI is between",pred_low,"and",pred_high))



url_coal <- "https://www.counterpointstat.com/uploads/1/1/9/3/119383887/coal.ash.txt" 
coalash <- read.table(url_coal,header = T)


coords_coal <- as.matrix(coalash[,c("x","y")]) 
coal <- coalash$coal
vario.coal <- variog(coords = coords_coal, data = coal, uvec = (seq(0, length = bins)))

plot(vario.coal)

eyefit(vario.coal,silent=TRUE)

fit.coal<- variofit(vario.coal,cov.model="exponential",
                    fix.nugget=FALSE,max.dist = 1/17, nugget=1.2)

fit.coal


coal_point<-krige.conv(coords = coords_coal, data = coal,loc=c(length(coal),1),
                  krige=krige.control(cov.pars=c(1.14,1/17),cov.model="exponential",nugget=1.2))
coal_point

pred_low <-coal_point$predict - 2*sqrt(coal_point$krige.var)
pred_high <-coal_point$predict + 2*sqrt(coal_point$krige.var)
print(paste("The PI is between",pred_low,"and",pred_high))



given_coords <- coords[1:dim(coords)[1]-1,]
myscallops.covmat<-varcov.spatial(coords = given_coords,
                                  cov.model = "exponential", 
                                  nugget = 1.67,
                                  cov.pars = c(5.34,0.88))
# just print myscallops.covmat will be too long!
# myscallops.covmat
# setup covariance matrix with point for prediction
gamma_all<-varcov.spatial(coords = coords,
                          cov.model = "exponential", 
                          nugget = 1.67,
                          cov.pars = c(5.34,0.88))
# we are interested in last column
gamma <- gamma_all$varcov[,ncol(gamma_all$varcov)]
gamma <- gamma[1:length(gamma)-1]

z <- lgcatch[1:length(lgcatch-1)]

mu <- mean(z)
m <- rep(mu, length(z))
m <- matrix(m,nrow=length(z),ncol=1)
z <- matrix(z,nrow=length(z),ncol=1)

y_pred<- mu +t(gamma)%*%solve(myscallops.covmat$varcov)%*%(z-m)
y_pred
