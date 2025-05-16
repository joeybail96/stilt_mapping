library(ggplot2)
library(maps)
library(patchwork)
library(raster)
library(ncdf4)


# Events
t_start <- c('2021-10-18 05:00:00', '2022-03-04 08:00:00', '2022-03-25 02:00:00', '2022-04-25 10:00:00')
t_end   <- c('2021-10-18 19:00:00', '2022-03-04 22:00:00', '2022-03-26 10:00:00', '2022-04-28 22:00:00')
names<-c("18102021.png", "04032022.png","25032022.png","26042022.png")

#### Basemap
map<-map_data("state")
p<-ggplot()+geom_map(data = map, map = map, 
                     aes(long, lat, map_id = region),  
                     color = "black", fill = NA, size = 0.2)+lims( y = c(30,50), x = c(-125, -100))+ coord_fixed(2)+
                    theme(panel.background = element_blank())

### 

for(i in 1 :length(t_start)){
  print(i)
  # Get foots
  run_times <- seq(from = as.POSIXct(t_start[i], tz = 'UTC'),
                   to   = as.POSIXct(t_end[i], tz = 'UTC'),
                   by   = 'hour')
  
  run_times<-run_times + 7*3600
  run_times<-gsub("-","", run_times)
  run_times<-gsub(" ","", run_times)
  run_times<-gsub(":","", run_times)
  run_times<-substr(run_times, 1, 12)
  foots<-paste0("/uufs/chpc.utah.edu/common/home/lin-group19/KaiW/dust_SPL/out/footprints/",run_times, "_-106.744_40.455_5_foot.nc")
  
  before_foots<-foots[1:4]
  after_foots<-foots[(length(foots) - 3) : length(foots)]
  foots<-foots[!( foots %in% c(before_foots, after_foots))]
  
  
  # Before
  before<-NULL
  for(i in before_foots){
  part<-readRDS(particles)
  save<-part$particle
  part<-part$particle
  
  # Subset particles
  part<-part[part$time %in% seq(-1, -14400, -60), ]
  part<-part[part$indx %in% sample(1:1000, 200, replace = F), ]
  part$zagl[part$zagl > 15000]<-15000
  
  }
  
  
  
  # Collapse to average foots
  before<-stack(before_foots)
  before<-stackApply(before, rep(1, dim(before)[3]), mean)
  before<-as.data.frame(rasterToPoints(before))
  colnames(before)[3]<-"foot"
  before$alp<-0.7
  before$alp[before$foot == 0]<-0
  
  for(j in foots){
    j_foot<-brick(j)
    j_foot<-stackApply(j_foot, rep(1, dim(j_foot)[3]), mean)
    if(j == foots[1]){
      during<-j_foot
      next()
    }
    during<-stack(during, j_foot)
  }
  during<-stackApply(during, rep(1, dim(during)[3]), mean)
  during<-as.data.frame(rasterToPoints(during))
  colnames(during)[3]<-"foot"
  during$alp<-0.7
  during$alp[during$foot == 0]<-0
  
  after<-stack(after_foots)
  after<-stackApply(after, rep(1, dim(after)[3]), mean)
  after<-as.data.frame(rasterToPoints(after))
  colnames(after)[3]<-"foot"
  after$alp<-0.7
  after$alp[after$foot == 0]<-0
  
  max<-max(c(max(before$foot), max(during$foot), max(after$foot)))
  
  b<-p+geom_raster(aes(x = x, y = y, fill = log10(foot), alpha = alp), before)+ coord_cartesian()+
    scale_fill_gradientn(colors = c("transparent","green", "yellow", "orange", "red", "purple"), 
                         values = c(0,0.04, 0.196,0.3, 0.668,1), na.value = "transparent", limits = c(-8, log10(max)))+
    scale_alpha_identity()+
   ggtitle("before")
  
  d<-p+geom_raster(aes(x = x, y = y, fill = log10(foot), alpha = alp), during)+ coord_cartesian()+
    scale_fill_gradientn(colors = c("transparent","green", "yellow", "orange", "red", "purple"), 
                         values = c(0,0.04, 0.196,0.3, 0.668,1), na.value = "transparent", limits = c(-8, log10(max)))+
    scale_alpha_identity()+
    ggtitle("during")
  
  a<-p+geom_raster(aes(x = x, y = y, fill = log10(foot), alpha = alp), after)+ coord_cartesian()+
    scale_fill_gradientn(colors = c("transparent","green", "yellow", "orange", "red", "purple"), 
                         values = c(0,0.04, 0.196,0.3, 0.668,1), na.value = "transparent", limits = c(-8, log10(max)))+
    scale_alpha_identity()+
    ggtitle("after")
  
  # Combine
  plot<-b+d+a
  
  ggsave(names[i], width = 12, height = 4, unit = "in")
  gc()
}







q<-p+geom_point(aes(x = mean_lon, y = mean_lat, color = mean_zagl), trajs, size = 0.5)+
  scale_color_gradientn(colors = c("red", "orange", "yellow", "green", "blue", "purple"), limits = c(0, 7000))+
  geom_map(data = map, map = map, 
           aes(long, lat, map_id = region),  
           color = "black", fill = NA, size = 0.2)+lims( y = c(20,70), x = c(50, 310))+ coord_fixed(2)+labs(title = event)
event<-gsub(" ", "_", event)
filename<-paste0("/uufs/chpc.utah.edu/common/home/lin-group19/KaiW/mercury_SPL/figs/cropped_", event, ".png")
ggsave(filename, height = 8, width = 12, units = "in")