library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)

states <- map_data("state")
states <- data.frame(states)

state.sat.scores <- read.table("data/state-sat.dat", header=F)
colnames(state.sat.scores) <- c("STATE","VERBAL","MATH","ELIGIBLE")

state.sat.scores[]
x = ((state.sat.scores$STATE=="alaska") |
       (state.sat.scores$STATE=="hawaii") | 
       (state.sat.scores$STATE=="us"))
index = c(1:nrow(state.sat.scores))[x]

state.sat.scores <- data.frame(state.sat.scores[-index,])

usa.state <- map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"),
                   function(x) x[1])

usa.poly <- map2SpatialPolygons(usa.state, 
                                IDs=state.ID)
usa.nb <- poly2nb(usa.poly)
usa.listb <- nb2listw(usa.nb, style="B")
usa.listw <- nb2listw(usa.nb, style="W")

# train CAR model

x = ((state.sat.scores$STATE=="alaska") |
       (state.sat.scores$STATE=="hawaii") | 
       (state.sat.scores$STATE=="us"))

index = c(1:nrow(state.sat.scores))[x]
state.sat.scores = state.sat.scores[-index,]

# binnary weights
stat.sat.car.b <- spautolm(ELIGIBLE~ VERBAL,
                           data=state.sat.scores,
                           family="CAR", 
                           listw=usa.listb, 
                           zero.policy=TRUE)
summary(stat.sat.car.b)
# row-normalized weights
stat.sat.car.w <- spautolm(ELIGIBLE~ VERBAL,
                           data=state.sat.scores,
                           family="CAR", 
                           listw=usa.listw, 
                           zero.policy=TRUE)
summary(stat.sat.car.w)


state.sat.scores.fitted = fitted(state.sat.scores)
state.sat.scores$fitted.car = state.sat.scores.fitted




