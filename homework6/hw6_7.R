#load libraries
library(spBayes)
library(MBA)
library(geoR)
library(fields)
library(sp)
library(maptools)
library(rgdal)
library(classInt)
library(lattice)
library(dplyr)


bat_df <- read.table("homework6/batonrouge.dat",header = T)

## Extract the coordinates
coords <- as.matrix(bat_df[,c("Longitude","Latitude")])

bins = 50
max.dist <- 0.7*max(iDist(coords))
log.selling.vario <- variog(coords = coords, data = bat_df$logSellingPr, 
                          uvec = (seq(0, max.dist, length = bins)))
plot(log.selling.vario)
eyefit(log.selling.vario,silent=TRUE)


## OLS 
lm.log.selling = lm(logSellingPr~LivingArea + Age + OtherArea + Bathrooms, data=bat_df)


## Obtain OLS residuals
selling.resid = resid(lm.log.selling)


selling.resid.vario <- variog(coords = coords, data = selling.resid, 
                            uvec = (seq(0, max.dist, length = bins)))
plot(selling.resid.vario)
eyefit(selling.resid.vario,silent=TRUE)

# the values of sigma square, tau sqaure  and nugget are from the above ressults
point<-krige.conv(coords = coords, data = bat_df$logSellingPr,loc=c(length(bat_df$logSellingPr),1),
                  krige=krige.control(cov.pars=c(0.22,0.17),cov.model="exponential",nugget=0.05))
point
pred_low <-point$predict - 2*sqrt(point$krige.var)
pred_high <-point$predict + 2*sqrt(point$krige.var)
print(paste("The 95% confident PI is between",pred_low,"and",pred_high))




n.samples = 1000

log.selling.mcmc <- spLM(logSellingPr~LivingArea + Age + OtherArea + Bathrooms, 
                         data=bat_df,
                         coords=coords, 
                         starting=list("phi"=0.17,"sigma.sq"=0.22,
                                       "tau.sq"=0.05), 
                         tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
                         priors=list("phi.Unif"=c(3/1500, 10),
                                     "sigma.sq.IG"=c(0.1, 0.1),
                                     "tau.sq.IG"=c(0.1, 0.1)), 
                         cov.model="exponential",n.samples=n.samples)


round(summary(mcmc(log.selling.mcmc$p.theta.samples))$quantiles,3)

burn.in <- floor(0.75*n.samples)
log.selling.mcmc <- spRecover(log.selling.mcmc, start=burn.in, thin=2)
Dic1 = spDiag(log.selling.mcmc,start=burn.in,verbose=FALSE)
Dic1


beta.samples = log.selling.mcmc$p.beta.recover.samples
w.samples = log.selling.mcmc$p.w.recover.samples

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)
