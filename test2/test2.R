library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)
library(geoR)
library(spBayes)

columbus.poly <- readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep")[1])


coords <- as.matrix(cbind(columbus.poly$X,columbus.poly$Y))
HOVAL <- columbus.poly$HOVAL
bins = 50
max.dist <- max(iDist(coords))
HOVAL.vario <- variog(coords = coords, data = HOVAL, 
                       uvec = (seq(0, max.dist, length = bins)))
plot(HOVAL.vario)

eyefit(HOVAL.vario,silent=TRUE)


```{r 6x_Kriging}
point<-krige.conv(coords = coords, data = lgcatch,loc=c(length(lgcatch),1),
                  krige=krige.control(cov.pars=c(488.5,0.1767),
                                      cov.model="exponential",
                                      nugget=0.3802))
point
pred_low <-point$predict - 2*sqrt(point$krige.var)
pred_high <-point$predict + 2*sqrt(point$krige.var)
print(paste("The PI is between",pred_low,"and",pred_high))
```

