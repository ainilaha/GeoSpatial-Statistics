## R programs for BCG Second Edition, pages 88--94.
##set.seed(1234)
library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)

##Draw the maps using the maps function
library(maps)
library(classInt)



## Constructing neighbors using distances: Columbus example
columbus.poly <- readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep")[1])

##Distance based neighbors in spdep
columbus.coords = coordinates(columbus.poly)
columbus.knn = knearneigh(columbus.coords)
columbus.knn2nb = knn2nb(columbus.knn)
columbus.dist.list = nbdists(columbus.knn2nb, columbus.coords)
columbus.dist.vec = unlist(columbus.dist.list)
columbus.dist.max = max(columbus.dist.vec)
columbus.dnn.nb = dnearneigh(columbus.coords, 0, 0.25*columbus.dist.max)

##Form a listw object using the distance-based nearest neighbors 
columbus.dnn.listw = nb2listw(columbus.dnn.nb, style="B", zero.policy=TRUE)

##SAR model regressing HOUSE_VAL+INCOME using distance-based nearest neighbors
columbus.dnn.sar.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="SAR", listw=columbus.dnn.listw, zero.policy=TRUE)


##CAR model regressing HOUSE_VAL+INCOME using distance-based nearest neighbors
columbus.dnn.car.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="CAR", listw=columbus.dnn.listw, zero.policy=TRUE)
columbus.dnn.car.fitted = fitted(columbus.dnn.car.out)
columbus.poly$fitted.dnn.car = columbus.dnn.car.fitted
