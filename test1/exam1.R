ppm <- c(14.3, 10.05,9.91,
         10.15,12.31,10.39,
         9.61,11.10,9.16,
         9.39,12.05)
X <- c(0.35,0.80,0.41,
       0.05,0.21,0.04,
       0.55,0.11,0.49,
       0.15,0.38)
Y <- c(0.85,0.44,0.79,
       0.08,0.58,0.92,
       0.32,0.85,0.13,
       0.08,0.83)
ppm.df <- data.frame("ppm"=ppm,"X"=X,"Y"=Y)
head(ppm.df)
coords <- as.matrix(ppm.df[,c("X","Y")]) 

library(akima)
int.ppm <-interp.new(ppm.df$X,ppm.df$Y,ppm.df$ppm,extrap=TRUE)
image(int.ppm)
contour(int.ppm, add=T)
persp(int.ppm)
persp(int.ppm,theta = -70, phi = 45, d=1)
persp(int.ppm,theta = -70, phi = 45, d=1,col="Yellow",shade=0.75,border=NA)

library(MBA)
library(fields) ## For using the image.plot function
x.res <- 100
y.res <- 100
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
col.pal <- col.br(5)


surf <- mba.surf(cbind(coords, ppm),
                 no.X = x.res, no.Y = y.res, h = 5,
                 m = 2, extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r",
           xlab = "X", ylab = "Y", col = col.br(25))
contour(surf, add=T)

library(rgl)
col <- rbind(0, cbind(matrix(drape.color(surf[[3]],
                                         col = col.br(25)), x.res - 1, y.res - 1), 0))
surface3d(surf[[1]], surf[[2]], surf[[3]], col = col)
axes3d()

title3d(main = "PPM", xlab = "X", ylab = "Y", zlab = "ppm")
drape.plot(surf[[1]], surf[[2]], surf[[3]], col = col.br(150),
           theta = 225, phi = 50, border = FALSE, add.legend = FALSE, 
           xlab = "X", ylab = "Y", zlab = "ppm")
image.plot(zlim = range(surf[[3]], na.rm = TRUE), legend.only = TRUE, horizontal = FALSE)

