### R code from BCG Second Edition, Section 9.8, pp. 297--301. 
### Encoding: UTF-8

###################################################
### code chunk number 1: libs
###################################################
options(width = 67) 
library(spBayes)
library(MBA)
library(geoR)
library(fields)


###################################################
### code chunk number 2: getData
###################################################
dat <- read.table("CostaRica/T4.csv", header=T, sep=",")
coords <- as.matrix(dat[,c("X","Y")])
nut.names <- c("Ca","K","Mg")
log.nut <- log(dat[,nut.names])


###################################################
### code chunk number 3: nutSurf
###################################################
par(mfrow=c(2,2))
for(i in 1:length(nut.names)){
  surf <- mba.surf(cbind(coords,data=log.nut[,i]), no.X=100, no.Y=100)$xyz.est
  image.plot(surf, main=paste("Log ",nut.names[i],sep=""))
  points(coords)
}


###################################################
### code chunk number 4: geoRVariog
###################################################
max <- 0.25*max(as.matrix(dist(dat[,c("X","Y")])))
bins <- c(9,8,9)

par(mfrow=c(3,1))
for(i in 1:length(nut.names)){
  vario <- variog(coords=coords,data=log.nut[,i], uvec=(seq(0,max, length=bins[i])))

  fit <-variofit(vario, ini.cov.pars=c(0.3, 20/-log(0.05)), ##sigma^2 and 1/phi 
                     cov.model="exponential", minimisation.function="nls", weights="equal")

  plot(vario, pch=19, main=paste("Log ",nut.names[i],sep=""))
  lines(fit)
  abline(h=fit$nugget, col="blue")
  abline(h=fit$cov.pars[1]+fit$nugget, col="green")
  abline(v=-log(0.05)*fit$cov.pars[2], col="red3")
}


###################################################
### code chunk number 5: spMvLM
###################################################
q <- 3
##Specify starting values and collect samples
A.starting <- diag(0.1,q)[lower.tri(diag(1,q), TRUE)]
Psi.starting <- rep(0.5, q)

A.tuning <- rep(0.0005,length(A.starting))
Psi.tuning <- rep(0.0005, q)

n.samples <- 10000

starting <- list("phi"=rep(3/20,q),"A"=A.starting,"Psi"=Psi.starting)
tuning <- list("phi"=rep(0.1,q),"A"=A.tuning, "Psi"=Psi.tuning)
priors <- list("phi.Unif"=list(rep(3/60,q), rep(3/10,q)), "K.IW"=list(q+1, diag(0.001,q)), "Psi.IG"=list(rep(2,q), c(0.05,0.08,0.05)))

m.1 <- spMvLM(list(Ca~1,K~1,Mg~1), coords=coords, data=log.nut,
              starting=starting, tuning=tuning, priors=priors,
              cov.model="exponential",  n.samples=n.samples, n.report=2000)

##other starting
A.starting <- diag(0.05,q)[lower.tri(diag(1,q), TRUE)]
Psi.starting <- rep(0.1, q)

m.2 <- spMvLM(list(Ca~1,K~1,Mg~1), coords=coords, data=log.nut,
              starting=starting, tuning=tuning, priors=priors,
              cov.model="exponential",  n.samples=n.samples, n.report=2000)

round(summary(m.1$p.theta.samples)$quantiles[,c(1,3,5)],3)

par(mfrow=c(1,2))
plot(mcmc.list(m.1$p.theta.samples, m.2$p.theta.samples), density=FALSE)


###################################################
### code chunk number 6: spMvLM-range
###################################################
burn.in <- 0.75*n.samples
p.theta.samples <- as.matrix(window(mcmc.list(m.1$p.theta.samples, m.2$p.theta.samples), window=burn.in))
n.samples <- nrow(p.theta.samples)

fn <- function(d, a, phi){
  0.05 - sum(a^2*exp(-phi*d))/sum(a^2)
}

get.A <- function(K, q){
  A <- matrix(0, q, q)
  A[lower.tri(A,TRUE)] <- K
  A[upper.tri(A,FALSE)] <- t(A)[upper.tri(A,FALSE)]
  t(chol(A))
}

eff.range <- matrix(0, 3, n.samples)

for(s in 1:n.samples){
  A <- get.A(p.theta.samples[s, c("K[1,1]","K[2,1]","K[3,1]","K[2,2]","K[3,2]","K[3,3]")], q)
  phi <- p.theta.samples[s, c("phi[1]", "phi[2]", "phi[3]")] 
  
  for(r in 1:q){
    eff.range[r,s] <- uniroot(fn, lower=0, upper=1000, tol=1e-15, a=A[i,1:r], phi=phi[1:r])$root
  }
}

rownames(eff.range) <- paste("Eff. range ",1:q,sep="")

round(summary(mcmc(cbind(p.theta.samples, t(eff.range))))$quantiles[,c(1,3,5)],3)


###################################################
### code chunk number 7: spLMW
###################################################
Ca.resids <- resid(lm(Ca~1, data=log.nut))
K.resids <- resid(lm(K~1, data=log.nut))
Mg.resids <- resid(lm(Mg~1, data=log.nut))

m.1 <- spRecover(m.1, start=burn.in, thin=10)

w <- rowMeans(m.1$p.w.recover.samples)
w.Ca <- w[seq(1,length(w),q)]
w.K <- w[seq(2,length(w),q)]
w.Mg <- w[seq(3,length(w),q)]

res <- 100
par(mfrow=c(3,2))
surf.r <- mba.surf(cbind(coords,Ca.resids), no.X=res, no.Y=res, extend=FALSE)$xyz.est
surf.w <- mba.surf(cbind(coords,w.Ca), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(c(surf.r[["z"]], surf.w[["z"]]), na.rm=TRUE)

image.plot(surf.r, zlim=z.lim, main="Ca lm residuals"); points(coords)
image.plot(surf.w, zlim=z.lim, main="Ca spatial effects"); points(coords)

surf.r <- mba.surf(cbind(coords,K.resids), no.X=res, no.Y=res, extend=FALSE)$xyz.est
surf.w <- mba.surf(cbind(coords,w.K), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(c(surf.r[["z"]], surf.w[["z"]]), na.rm=TRUE)

image.plot(surf.r, zlim=z.lim, main="K lm residuals");points(coords)
image.plot(surf.w, zlim=z.lim, main="K spatial effects");points(coords)

surf.r <- mba.surf(cbind(coords,K.resids), no.X=res, no.Y=res, extend=FALSE)$xyz.est
surf.w <- mba.surf(cbind(coords,w.Mg), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(c(surf.r[["z"]], surf.w[["z"]]), na.rm=TRUE)

image.plot(surf.r, zlim=z.lim, main="Mg lm residuals"); points(coords)
image.plot(surf.w, zlim=z.lim, main="Mg spatial effects"); points(coords)


###################################################
### code chunk number 8: spPredict
###################################################
x.range <- range(coords[,1])
y.range <- range(coords[,2])

pred.coords <- expand.grid(seq(x.range[1], x.range[2], by=4),
                           seq(y.range[1], y.range[2], by=4))

m <- nrow(pred.coords)

pred.X <- mkMvX(list(matrix(1,m,1), matrix(1,m,1), matrix(1,m,1)))

nut.pred <- spPredict(m.1, start=burn.in, thin=10, pred.coords=pred.coords, pred.covars=pred.X)


###################################################
### code chunk number 9: spPredMean
###################################################
y.pred.mu <- apply(nut.pred$p.y.predictive.samples, 1, mean)
y.pred.sd <- apply(nut.pred$p.y.predictive.samples, 1, sd)

Ca.pred.mu <- y.pred.mu[seq(1,length(y.pred.mu),q)]
K.pred.mu <- y.pred.mu[seq(2,length(y.pred.mu),q)]
Mg.pred.mu <- y.pred.mu[seq(3,length(y.pred.mu),q)]

Ca.pred.sd <- y.pred.sd[seq(1,length(y.pred.sd),q)]
K.pred.sd <- y.pred.sd[seq(2,length(y.pred.sd),q)]
Mg.pred.sd <- y.pred.sd[seq(3,length(y.pred.sd),q)]

nut.pred.grid <- as.data.frame(list(x=pred.coords[,1], y=pred.coords[,2],
                                     Ca.mu=Ca.pred.mu, K.mu=K.pred.mu, Mg.mu=Mg.pred.mu,
                                     Ca.sd=Ca.pred.sd, K.sd=K.pred.sd, Mg.sd=Mg.pred.sd))
                                     
coordinates(nut.pred.grid) <- c("x", "y") # promote to SpatialPointsDataFrame
gridded(nut.pred.grid) <- TRUE            # promote to SpatialGridDataFrame

toImage <- function(x){as.image.SpatialGridDataFrame(x)}

res <- 100
par(mfrow=c(3,2))

surf <- mba.surf(cbind(coords,log.nut[,"Ca"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(surf[["z"]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", main="Interpolation of observed Ca")
points(coords)

image.plot(toImage(nut.pred.grid["Ca.mu"]), xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean of Ca prediction")
points(coords)

surf <- mba.surf(cbind(coords,log.nut[,"K"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(surf[["z"]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", main="Interpolation of observed K")
points(coords)

image.plot(toImage(nut.pred.grid["K.mu"]), xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean of K prediction")
points(coords)

surf <- mba.surf(cbind(coords,log.nut[,"Mg"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(surf[["z"]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", main="Interpolation of observed Mg")
points(coords)

image.plot(toImage(nut.pred.grid["Mg.mu"]), xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean of Mg prediction")
points(coords)


###################################################
### code chunk number 10: spPredSD
###################################################
par(mfrow=c(3,2))
surf <- mba.surf(cbind(coords,log.nut[,"Ca"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Interpolation of observed Ca")
points(coords)

surf <- mba.surf(cbind(pred.coords,Ca.pred.sd), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Standard deviation of Ca prediction")
points(coords)

surf <- mba.surf(cbind(coords,log.nut[,"K"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Interpolation of observed K")
points(coords)

surf <- mba.surf(cbind(pred.coords,K.pred.sd), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Standard deviation of K prediction")
points(coords)

surf <- mba.surf(cbind(coords,log.nut[,"Mg"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Interpolation of observed Mg")
points(coords)

surf <- mba.surf(cbind(pred.coords,Mg.pred.sd), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Standard deviation of Mg prediction")
points(coords)


