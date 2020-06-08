library(openair)
library(epitools)
library(lubridate)
library(maps)
library(mapproj)
library(ggplot2)
library(magick)

for (m in 01:31){
folder <- list.dirs("D:/OneDrive - University Of Houston/TG report/alldata/")
filename <- paste0("merged",sprintf("%02.f",m),".RData")
setwd(folder)
load (filename)
}

merged_final01 <- merged_final1
merged_final02 <- merged_final2
merged_final03 <- merged_final3
merged_final04 <- merged_final4
merged_final05 <- merged_final5
merged_final06 <- merged_final6
merged_final07 <- merged_final7
merged_final08 <- merged_final8
merged_final09 <- merged_final9


  
#rm(data11,data12,data13,data14,data15,data16,data17,data18,data19)  
rm(merged_final1,merged_final2,merged_final3,merged_final4,merged_final5,merged_final6,merged_final7,merged_final8,merged_final9)
df_list <- mget(ls(pattern = "merged_final[0-9]"))
all_data <- do.call("rbind",df_list)
lr <- which(grepl('00:00:00', all_data$date))
all_data <- all_data[-c(lr),]
#save("all_data", file = "all_data.RData")
all_data$height <- as.numeric(as.character(all_data$height))/1000
all_data$conc <- as.numeric(as.character(all_data$conc))

dwp <- do.call("rbind", list(merged_final10, merged_final11, merged_final12,
                                   merged_final13,merged_final14, merged_final15,
                                   merged_final16))
lrdwp <- which(grepl('00:00:00', dwp$date))
dwp <- dwp[-c(lrdwp),]

sp <- do.call("rbind", list(merged_final17, merged_final18, merged_final19,
                             merged_final20,merged_final21, merged_final22))
lrsp <- which(grepl('00:00:00', sp$date))
sp <- sp[-c(lrsp),]

epp <- do.call("rbind", list(merged_final25, merged_final26, merged_final27,
                             merged_final28))
lrepp <- which(grepl('00:00:00', epp$date))
epp <- epp[-c(lrepp),]

#all_data2 <- sapply(all_data,mode)
# transform(epp, height = as.numeric(height))
dwp$height <- as.numeric(as.character(dwp$height))/1000
sp$height <- as.numeric(as.character(sp$height))/1000
epp$height <- as.numeric(as.character(epp$height))/1000



# library(mapdata)
# trajPlot(subset(all_data, lon > 124 & lon < 131 & lat >33 & lat <39),
#          pollutant = "conc",lon = "lon",lat = "lat",map.alpha = 0.4,npoints = 5, orientation = c(0,90,0))
# trajPlot(epp,lon = "lon",lat = "lat",map.alpha = 0.4,npoints = NA, orientation = c(0,90,0))
# trajPlot(all_data,pollutant = "height",lon = "lon",lat = "lat",map.alpha = 0.4,npoints = 5, orientation = c(0,90,0))
# trajLevel(all_data, pollutant = "conc",percentile = 50,lon = "lon",lat = "lat",map.alpha = 0.3,npoints = 5, orientation = c(0,90,0))
# trajLevel(all_data, pollutant = "conc",method = "hexbin", statistic = "frequency",lon = "lon",lat = "lat",map.alpha = 0.3,npoints = 5, orientation = c(0,90,0))
# trajLevel(all_data, pollutant = "conc",statistic = "difference",lon = "lon",lat = "lat",map.alpha = 0.3,npoints = 5, orientation = c(0,90,0))
# trajLevel(all_data, pollutant = "conc",statistic = "pscf",lon = "lon",lat = "lat",map.alpha = 0.3,npoints = 5, orientation = c(0,90,0))

## MAIN PLOTS


tiff("epp_large.tiff", width = 5, height = 4, units = 'in', res = 300)
trajPlot(subset(epp, epp$height < 14.5),lon = "lon",lat = "lat",parameters = NULL, pollutant = "height",
         map.alpha = 0.3,npoints = 4, par.settings=list(fontsize=list(text=10)),
         projection = "mercator",xlim = c(105,135),ylim=c(20,55),xlab="Lon",ylab="Lat",
         map = TRUE,cols = c("blue","gold","red"),
         grid.col = "dimgrey",main = "")
dev.off()
folder <- list.dirs("D:/OneDrive - University Of Houston/TG report/alldata/all new results/")
setwd(folder)
frink <- image_read("all_data_large.tiff")
frink<- image_annotate(frink, "Altitude (km)", size = 50, gravity = "northeast", color = "black")
frink<- image_annotate(frink, "(a)", size = 75, gravity = "northwest", color = "black")
  img <- image_draw(frink)
rect(710, 573, 1133, 853, border = "red", lty = "dashed", lwd = 5)
image_write(img, path = "all_data_large_proccesed2.tiff", format = "tiff")
dev.off()

tiff("epp_small.tiff", width = 5, height = 4, units = 'in', res = 300)
trajPlot(subset(epp, epp$height < 14.5),lon = "lon",lat = "lat",parameters = NULL, pollutant = "height",
         map.alpha = 0.3,npoints = 4, par.settings=list(fontsize=list(text=10)),
         projection = "mercator",xlim = c(120,130),ylim=c(30,45),xlab="Lon",ylab="Lat",
         map = TRUE,cols = c("blue","gold","red"),
         grid.col = "grey70",main = "")
dev.off()
frink <- image_read("all_data_small.tiff")
frink<- image_annotate(frink, "Altitude (km)", size = 50, gravity = "northeast", color = "black")
frink<- image_annotate(frink, "(b)", size = 75, gravity = "northwest", color = "black")
image_write(frink, path = "all_data_small_proccesed2.tiff", format = "tiff")
dev.off()



tiff("Plot6.tiff", width = 5, height = 4, units = 'in', res = 300)
clust <- trajCluster(all_data,n.clusters = 5,by.type = TRUE,origin = TRUE,par.settings=list(fontsize=list(text=10)),
                     method = "Euclid",lon = "lon",lat = "lat",parameters = NULL,smooth = TRUE,
                     map.alpha = 0.3, proj = "mercator",npoints = NA,
                     xlim = c(105,135),ylim=c(20,55),xlab="Lon",ylab="Lat",
                     grid.col = "dimgrey",main = "Cluster Euclid")
dev.off()
clust <- trajCluster(subset(all_data, hour.inc == 0),n.clusters = 5, orientation = c(0,90,0), method = "Angle",smooth = TRUE)

trajPlot(all_data, pollutant = "conc",lon = "lon",lat = "lat",projection = "mercator",parameters = NULL,map.alpha = 0.4,npoints = 5)
trajPlot(all_data, group = "height",lon = "lon",lat = "lat",map.alpha = 0.4,npoints = 5, orientation = c(0,90,0))

trajPlot(all_data, pollutant = "conc",type = "height",lon = "lon",lat = "lat",map.alpha = 0.3, orientation = c(0,90,0),grid.col = "transparent")
trajLevel(clust$data,group = "cluster",lon = "lon",lat = "lat",map.alpha = 0.4,npoints = 5, orientation = c(0,90,0))


traj <- importTraj(site = "london", year = 2009)
## calculate clusters
clust2 <- trajCluster(traj, n.clusters = 5)
head(clust2$data) ## note new variable 'cluster'
## use different distance matrix calculation, and calculate by season
traj2 <- trajCluster(traj, method = "Angle", n.clusters = 5)
# }