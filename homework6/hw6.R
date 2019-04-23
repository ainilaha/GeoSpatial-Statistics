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




lithology_df <- read.csv("homework6/lithology.csv",header = T)
# drop the NA or miss data
lithology_df <- lithology_df %>% 
  filter( !is.na(Thickness_ft) & 
          !is.na(Surf_Elevation_ft_amsl) &
         !is.na(A_B_Elevation_ft_amsl)) %>%
  distinct(Easting_ft,Northing_ft,.keep_all = TRUE)
  
lithology_df$Surf_Elevation_ft_amsl <- as.numeric(lithology_df$Surf_Elevation_ft_amsl)
lithology_df$A_B_Elevation_ft_amsl <- as.numeric(lithology_df$A_B_Elevation_ft_amsl)

## Extract the coordinates
coords <- as.matrix(lithology_df[,c("Easting_ft","Northing_ft")])

x.res <- 200; y.res <- 200

surf <- mba.surf(cbind(coords, 
                       lithology_df$Thickness_ft), 
                 no.X=x.res, no.Y=y.res, h=5, 
                 m=2, extend=FALSE)$xyz.est
zlim.Thickness_ft <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r",
           yaxs = "r", xlab="Easting (ft)", 
           ylab="Northing (ft)",
           main="Thickness(ft)" )

points(coords)
contour(surf,add = T)

surf <- mba.surf(cbind(coords, 
                       lithology_df$Surf_Elevation_ft_amsl), 
                 no.X=x.res, no.Y=y.res, h=5, 
                 m=2, extend=FALSE)$xyz.est
zlim.Surf_Elevation_ft_amsl <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r",
           yaxs = "r", xlab="Easting (ft)", 
           ylab="Northing (ft)",
           main="Surf Elevation (ft amsl)")

points(coords)
contour(surf,add = T)


surf <- mba.surf(cbind(coords, 
                       lithology_df$A_B_Elevation_ft_amsl), 
                 no.X=x.res, no.Y=y.res, h=5, 
                 m=2, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r",
           yaxs = "r", xlab="Easting (ft)", 
           ylab="Northing (ft)",
           main="A-B Elevation(ft amsl)")
points(coords)
contour(surf,add = T)



log.thickness <- log(lithology_df$Thickness_ft)
bins = 50
max.dist <- 0.1*max(iDist(coords))
log.thick.vario <- variog(coords = coords, data = log.thickness, 
                       uvec = (seq(0, max.dist, length = bins)))
plot(log.thick.vario)
eyefit(log.thick.vario,silent=TRUE)





p <- 3 ## This is the number of columns in the design matrix
## Set the prior mean and precision for the regression
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)
sigma.sq.prior.shape <- 0.1 ## Set IG shape for sigma.sq (partial sill)
sigma.sq.prior.rate <- 0.1 ## Set IG scale for sigma.sq (partial sill)

## For use with bayesGeostatExact, do the following
phi <- 1/1100 ## Set the spatial range (from the variogram)
alpha <- 0.094/0.58 ## Set the nugget/partial-sill ratio

sp.exact <- bayesGeostatExact(
  log.thickness ~ Surf_Elevation_ft_amsl + A_B_Elevation_ft_amsl,
  data=lithology_df, coords=coords, n.samples=1000,
  beta.prior.mean=beta.prior.mean,
  beta.prior.precision=beta.prior.precision,
  cov.model="exponential",
  phi=phi, alpha=alpha,
  sigma.sq.prior.shape=sigma.sq.prior.shape,
  sigma.sq.prior.rate=sigma.sq.prior.rate,
  sp.effects=FALSE)
round(summary(sp.exact$p.samples)$quantiles,3)


phi <- 1/1000 ## Set the spatial range (from the variogram)
alpha <- 0.094/0.58 ## Set the nugget/partial-sill ratio
nu <- 0.5

## Run bayesGeostatExact to deliver exact posterior samples
sp.exact <- bayesGeostatExact(
  log.thickness ~ Surf_Elevation_ft_amsl + A_B_Elevation_ft_amsl,
  data=lithology_df, coords=coords, n.samples=1000,
  beta.prior.mean=beta.prior.mean,
  beta.prior.precision=beta.prior.precision,
  cov.model="matern",
  phi=phi, alpha=alpha,
  nu=nu,
  sigma.sq.prior.shape=sigma.sq.prior.shape,
  sigma.sq.prior.rate=sigma.sq.prior.rate,
  sp.effects=FALSE)
round(summary(sp.exact$p.samples)$quantiles,3)





## Run spLM to deliver MCMC samples from marginal posterior distributions
n.samples <- 1000
thick.sp <- spLM(log.thickness ~ Surf_Elevation_ft_amsl + A_B_Elevation_ft_amsl,
               data=lithology_df, coords=coords, 
               starting=list("phi"=3/1100,"sigma.sq"=0.08,"tau.sq"=1E-10), 
               tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=1E-10),
               priors=list("phi.Unif"=c(3/1500, 3/50), 
                           "sigma.sq.IG"=c(0.1,0.1),
                           "tau.sq.IG"=c(0.1, 0.1)), 
               cov.model="exponential",n.samples=n.samples)

round(summary(mcmc(thick.sp$p.theta.samples))$quantiles,3)



thick.sp2 <- spLM(log.thickness ~ Surf_Elevation_ft_amsl + A_B_Elevation_ft_amsl,
               data=lithology_df, coords=coords, 
               starting=list("phi"=3/1100,"sigma.sq"=0.08,"tau.sq"=0), 
               tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0),
               priors=list("phi.Unif"=c(3/1500, 3/50), 
                           "sigma.sq.IG"=c(0.1,0.1),
                           "tau.sq.IG"=c(0.1, 0.1)), 
               cov.model="exponential",n.samples=n.samples)

round(summary(mcmc(thick.sp2$p.theta.samples))$quantiles,3)
thick.sp <- spRecover(thick.sp, start=burn.in, thin=2)
Dic1 = spDiag(thick.sp,start=burn.in,verbose=FALSE)



## Recover spatial residuals using spRecover
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start=burn.in, thin=2)

## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples = bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

## Obtain trace plots for regression coefficients
par(mfrow=c(2,2))
plot(beta.samples, auto.layout=TRUE, density=FALSE)

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)
par(mfrow=c(1,1))
surf <- mba.surf(cbind(coords, w.hat.mu), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", zlim= main="log.thickness Mean Spatial Effects")
contour(surf,add = T)


## Run spLM to deliver MCMC samples from marginal posterior distributions
n.samples <- 1000
bef.sp <- spLM(log.thickness ~ Surf_Elevation_ft_amsl + A_B_Elevation_ft_amsl,
               data=lithology_df, coords=coords, 
               starting=list("phi"=3/1100,"sigma.sq"=0.08,"tau.sq"=0.02), 
               tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
               priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(0.1,0.1),"tau.sq.IG"=c(0.1, 0.1)), 
               cov.model="exponential",n.samples=n.samples)

round(summary(mcmc(bef.sp$p.theta.samples))$quantiles,3)


## Recover spatial residuals using spRecover
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start=burn.in, thin=2)

## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples = bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

## Obtain trace plots for regression coefficients
par(mfrow=c(2,2))
plot(beta.samples, auto.layout=TRUE, density=FALSE)

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)


surf <- mba.surf(cbind(coords, w.hat.mu),
                 no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", 
           main="Mean Spatial Effects")
contour(surf,add = T)

surf <- mba.surf(cbind(coords, w.hat.mu),
                 no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", 
           zlim = zlim.Thickness_ft,
           main="log(Thickness) Mean Spatial Effects Over Thickness(ft)")
contour(surf,add = T)


surf <- mba.surf(cbind(coords, w.hat.mu),
                 no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", 
           zlim = zlim.Surf_Elevation_ft_amsl,
           main="log(Thickness) Mean Spatial Effects Over Surf Elevation (ft amsl)")
contour(surf,add = T)
surf <- mba.surf(cbind(coords, w.hat.mu),
                 no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", 
           zlim = zlim.Surf_Elevation_ft_amsl,
           main="log(Thickness) Mean Spatial Effects Over A-B Elevation(ft amsl) ")
contour(surf,add = T)




