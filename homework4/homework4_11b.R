library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)

states <- map_data("state")
states <- data.frame(states)

state.sat.scores <- read.table("data/state-sat.dat", header=F)
colnames(state.sat.scores) <- c("STATE","VERBAL","MATH","ELIGIBLE")
x = ((state.sat.scores$STATE=="alaska") |
       (state.sat.scores$STATE=="hawaii") | 
       (state.sat.scores$STATE=="us"))
index = c(1:nrow(state.sat.scores))[x]

state.sat.scores <- data.frame(state.sat.scores[-index,])




states[,"sat_range"] <- 0
names(states)[6] <- "sat"
state.sat.scores$STATE <- unique(states$region)




for(i in 1:nrow(states))
{
  for (j in 1:nrow(state.sat.scores))
  {
    if(grepl(state.sat.scores[j,]$STATE, states[i,]$region))
    {
      sat <-  state.sat.scores[j,]$VERBAL2
      states[i,]$sat <- sat
      if(sat <= 503)
      {
        states[i,]$sat_range <- "<=503" 
      }else if(sat >503 && sat <=525 )
      {
        states[i,]$sat_range <- "504-525"
      }else if(sat >=525 && sat <562 )
      {
        states[i,]$sat_range <- "526-562"
      }else
      {
        states[i,]$sat_range <- ">563"
      }
    }
  }
}

ggplot(data = states) + 
  geom_polygon(aes(x = long, 
                   y = lat, fill = sat_range,group=group),color = "white") + 
  coord_fixed(1.3)

