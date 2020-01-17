library(rjson)
library(digest)
library(glue)
library(devtools)
#data <- read.csv("C:/Users/apouyaei/Desktop/dat.csv")
for (m in 2:32){
  folder <- list.dirs("/tng9/users/apouyaei/8TG/CMAQ-5.2/CCTM/scripts/res/")
  folder[m]
files <- list.files(folder[m],pattern = "*.dat")
setwd(folder[m])
data <- read.table(files[24],header = FALSE, sep = "", dec = ".")
#dat2 <- read.table("C:/Users/Arrman Pouyaei/Desktop/TG report/M01/TG2016123     0.dat",header = FALSE, sep = "", dec = ".")
names(data) <- c("jdate","hour","id","lat","lon","alt","conc")
#data <- data[!(data$id > 100000),]
#subtract the "South Korea" data
#data <- data[( (33.5 < data$lat) & (data$lat < 38.5) & (124.5 < data$lon) & (data$lon < 130.5)),]
#find the "Seoul" point
data_over_seoul <- data[((data$lat - 37.5821) < 0.1) & ((data$lat - 37.5821) > 0) & ((data$lon -  126.9706) < 0.1) & ((data$lon -  126.9706) > 0),]
filename = paste0("",sprintf("%02.f",(m-1)),".RData")
setwd(folder)
save(data_over_seoul,file = filename)
}

