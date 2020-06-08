library(rjson)
library(digest)
library(glue)
library(devtools)

for (m in 01:14) {
folder <- list.dirs("C:/Users/Arrman Pouyaei/Desktop/TG report/")
filename <- paste0("",sprintf("%02.f",m),".RData")
load (filename)
pathloc <- paste0("C:/Users/Arrman Pouyaei/Desktop/TG report/M",sprintf("%02.f",m))
filenames <- list.files(path=pathloc,full.names=TRUE)
All <- lapply(filenames,function(i){
  read.table(i,header = FALSE, sep = "", dec = ".")
})
data <- do.call(rbind.data.frame, All)
data <- data.frame(data)
names(data) <- c("jdate","hour","id","lat","lon","alt","conc")
    

for (n in 1:nrow(data_over_seoul)){
    nam <- paste0("data",n) 
    assign(nam,data[(data$id == data_over_seoul$id[n]),]) 
    }
rm(data,All)
setwd(folder)
save.image(file = paste0("alldata",sprintf("%02.f",m),".RData"))        
}