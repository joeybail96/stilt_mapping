library(ggplot2)
library(maps)
library(patchwork)
library(raster)
library(ncdf4)


# Events
t_start <- c('2025-02-04 18:00:00')
t_end   <- c('2025-02-05 18:00:00')
names<-c("04022025.png")

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
  
#  run_times<-run_times + 7*3600
#  run_times<-gsub("-","", run_times)
#  run_times<-gsub(" ","", run_times)
#  run_times<-gsub(":","", run_times)
#  run_times<-substr(run_times, 1, 12)
  run_times <- format(run_times + 7*3600, "%Y%m%d%H%M")
  foots<-paste0("/uufs/chpc.utah.edu/common/home/hallar-group2/climatology/stilt/dust_spl/out/2025_trajectories/footprints/",run_times, "_-106.744_40.455_5_foot.nc")

  before_foots<-foots[1:4]
  after_foots<-foots[(length(foots) - 3) : length(foots)]
  foots<-foots[!( foots %in% c(before_foots, after_foots))]
  
  

  
  
  
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
  
  output_dir <- "/uufs/chpc.utah.edu/common/home/hallar-group2/climatology/stilt/dust_spl/out/scripts/figures"
  ggsave(file.path(output_dir, names[i]), width = 12, height = 4, unit = "in")
  
  gc()
}







#q<-p+geom_point(aes(x = mean_lon, y = mean_lat, color = mean_zagl), trajs, size = 0.5)+
#  scale_color_gradientn(colors = c("red", "orange", "yellow", "green", "blue", "purple"), limits = c(0, 7000))+
#  geom_map(data = map, map = map, 
#           aes(long, lat, map_id = region),  
#           color = "black", fill = NA, size = 0.2)+lims( y = c(20,70), x = c(50, 310))+ coord_fixed(2)+labs(title = event)
#event<-gsub(" ", "_", event)
#filename<-paste0("/uufs/chpc.utah.edu/common/home/lin-group19/KaiW/mercury_SPL/figs/cropped_", event, ".png")
#ggsave(filename, height = 8, width = 12, units = "in")