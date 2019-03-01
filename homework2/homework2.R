library(spBayes)
library(classInt)
library(geoR)
url <- "https://www.counterpointstat.com/uploads/1/1/9/3/119383887/myscallops.txt"
myscallops <- read.table(url,header = T)

coords <- as.matrix(myscallops[,c("lat","long")])

strata <- myscallops$strata
lgcatch<- myscallops$lgcatch

plot(coords, pch=1, cex=sqrt(strata)/100, col="darkgreen", xlab="Lat", ylab="Long")
leg.vals <- round(quantile(strata),0)
legend("topleft", pch=1, legend=leg.vals, col="darkgreen",
       pt.cex=sqrt(leg.vals)/1000, bty="n", title="strata")

plot(coords, pch=1, cex=sqrt(lgcatch)/2, col="darkgreen", xlab="Lat", ylab="Long")
leg.vals <- round(quantile(lgcatch),0)
legend("topleft", pch=1, legend=leg.vals, col="darkgreen",
       pt.cex=sqrt(leg.vals)/1000, bty="n", title="lgcatch")



## Create a color palette for subsequent plots.
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

fixed <- classIntervals(strata, n = 4, style = "fixed", fixedBreaks = c(0,6260, 62290, 6310, max(strata) + 1))
fixed.col <- findColours(fixed, col.pal)
plot(coords, col = fixed.col, pch = 19, cex = 0.5, main = "Forestry tree size classes", xlab = "Lat", ylab = "Long")
legend("topleft", fill = attr(fixed.col, "palette"), legend = c("strata1", "strata2", "strata3", "strata4"), bty = "n")


## Load the MBA and fields libraries for creating surface interpolation plots
library(MBA)
library(fields) ## For using the image.plot function
x.res <- 100
y.res <- 100
surf <- mba.surf(cbind(coords, lgcatch), no.X = x.res, no.Y = y.res, h = 5, m = 2, extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "lat", ylab = "long", col = col.br(25))
contour(surf, add=T) 


library(rgl)
col <- rbind(0, cbind(matrix(drape.color(surf[[3]],
                                         col = col.br(25)), x.res - 1, y.res-1), 0))
surface3d(surf[[1]], surf[[2]], surf[[3]], col = col)
axes3d()

title3d(main = "strata", xlab = "Lat", ylab = "Long", zlab = "Strata")
drape.plot(surf[[1]], surf[[2]], surf[[3]], col = col.br(150),theta = 225, phi = 50, border = FALSE, add.legend = FALSE, xlab = "Lat", ylab = "Long", zlab = "Strata")
image.plot(zlim = range(surf[[3]], na.rm = TRUE), legend.only = TRUE, horizontal = FALSE)



library(rgl)
col <- rbind(0, cbind(matrix(drape.color(surf[[3]],
                                         col = col.br(25)), x.res - 1, y.res - 1), 0))
surface3d(surf[[1]], surf[[2]], surf[[3]], col = col)
axes3d()

title3d(main = "DBH", xlab = "Easting (m)", ylab = "Northing (m)", zlab = "DBH (cm)")
drape.plot(surf[[1]], surf[[2]], surf[[3]], col = col.br(150), theta = 225, phi = 50, border = FALSE, add.legend = FALSE, xlab = "Easting (m)", ylab = "Northing (m)", zlab = "DBH (cm)")
image.plot(zlim = range(surf[[3]], na.rm = TRUE), legend.only = TRUE, horizontal = FALSE)