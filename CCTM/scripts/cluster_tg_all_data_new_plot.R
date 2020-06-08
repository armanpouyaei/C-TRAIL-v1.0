library(openair)
library(epitools)
library(lubridate)
library(maps)
library(mapproj)
library(rjson)
library(digest)
library(glue)
library(devtools)
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/ggmap", ref = "tidyup")
library(ggmap)
library(ggplot2)
library(gganimate)
library(gifski)
theme_set(theme_bw())
library(gapminder)
library(tidyverse)
library(geosphere)

folder <- list.dirs("D:/OneDrive - University Of Houston/TG report/alldata/")
setwd(folder)
load ("hmmm.RData")


# all_data2$date2 <- cut(all_data2$hour.inc, c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
#                    labels = c("24", "23", "22", "21", "20", "19", "18", "17", "16", 
#                               "15", "14", "13", "12", "11", "10", "9", "8" , "7", "6",
#                               "5", "4", "3","2","1"),
#                    ordered_result = T)
all_data$date2[all_data$hour.inc == 0] <- "24:00"
all_data$date2[all_data$hour.inc == -1] <- "23:00"
all_data$date2[all_data$hour.inc == -2] <- "22:00"
all_data$date2[all_data$hour.inc == -3] <- "21:00"
all_data$date2[all_data$hour.inc == -4] <- "20:00"
all_data$date2[all_data$hour.inc == -5] <- "19:00"
all_data$date2[all_data$hour.inc == -6] <- "18:00"
all_data$date2[all_data$hour.inc == -7] <- "17:00"
all_data$date2[all_data$hour.inc == -8] <- "16:00"
all_data$date2[all_data$hour.inc == -9] <- "15:00"
all_data$date2[all_data$hour.inc == -10] <- "14:00"
all_data$date2[all_data$hour.inc == -11] <- "13:00"
all_data$date2[all_data$hour.inc == -12] <- "12:00"
all_data$date2[all_data$hour.inc == -13] <- "11:00"
all_data$date2[all_data$hour.inc == -14] <- "10:00"
all_data$date2[all_data$hour.inc == -15] <- "09:00"
all_data$date2[all_data$hour.inc == -16] <- "08:00"
all_data$date2[all_data$hour.inc == -17] <- "07:00"
all_data$date2[all_data$hour.inc == -18] <- "06:00"
all_data$date2[all_data$hour.inc == -19] <- "05:00"
all_data$date2[all_data$hour.inc == -20] <- "04:00"
all_data$date2[all_data$hour.inc == -21] <- "03:00"
all_data$date2[all_data$hour.inc == -22] <- "02:00"
all_data$date2[all_data$hour.inc == -23] <- "01:00"

all_data2 <- all_data[!duplicated(all_data$conc), ]

dwp$date2[dwp$hour.inc == 0] <- "24:00"
dwp$date2[dwp$hour.inc == -1] <- "23:00"
dwp$date2[dwp$hour.inc == -2] <- "22:00"
dwp$date2[dwp$hour.inc == -3] <- "21:00"
dwp$date2[dwp$hour.inc == -4] <- "20:00"
dwp$date2[dwp$hour.inc == -5] <- "19:00"
dwp$date2[dwp$hour.inc == -6] <- "18:00"
dwp$date2[dwp$hour.inc == -7] <- "17:00"
dwp$date2[dwp$hour.inc == -8] <- "16:00"
dwp$date2[dwp$hour.inc == -9] <- "15:00"
dwp$date2[dwp$hour.inc == -10] <- "14:00"
dwp$date2[dwp$hour.inc == -11] <- "13:00"
dwp$date2[dwp$hour.inc == -12] <- "12:00"
dwp$date2[dwp$hour.inc == -13] <- "11:00"
dwp$date2[dwp$hour.inc == -14] <- "10:00"
dwp$date2[dwp$hour.inc == -15] <- "09:00"
dwp$date2[dwp$hour.inc == -16] <- "08:00"
dwp$date2[dwp$hour.inc == -17] <- "07:00"
dwp$date2[dwp$hour.inc == -18] <- "06:00"
dwp$date2[dwp$hour.inc == -19] <- "05:00"
dwp$date2[dwp$hour.inc == -20] <- "04:00"
dwp$date2[dwp$hour.inc == -21] <- "03:00"
dwp$date2[dwp$hour.inc == -22] <- "02:00"
dwp$date2[dwp$hour.inc == -23] <- "01:00"

dwp2 <- dwp[!duplicated(dwp$conc), ]

sp$date2[sp$hour.inc == 0] <- "24:00"
sp$date2[sp$hour.inc == -1] <- "23:00"
sp$date2[sp$hour.inc == -2] <- "22:00"
sp$date2[sp$hour.inc == -3] <- "21:00"
sp$date2[sp$hour.inc == -4] <- "20:00"
sp$date2[sp$hour.inc == -5] <- "19:00"
sp$date2[sp$hour.inc == -6] <- "18:00"
sp$date2[sp$hour.inc == -7] <- "17:00"
sp$date2[sp$hour.inc == -8] <- "16:00"
sp$date2[sp$hour.inc == -9] <- "15:00"
sp$date2[sp$hour.inc == -10] <- "14:00"
sp$date2[sp$hour.inc == -11] <- "13:00"
sp$date2[sp$hour.inc == -12] <- "12:00"
sp$date2[sp$hour.inc == -13] <- "11:00"
sp$date2[sp$hour.inc == -14] <- "10:00"
sp$date2[sp$hour.inc == -15] <- "09:00"
sp$date2[sp$hour.inc == -16] <- "08:00"
sp$date2[sp$hour.inc == -17] <- "07:00"
sp$date2[sp$hour.inc == -18] <- "06:00"
sp$date2[sp$hour.inc == -19] <- "05:00"
sp$date2[sp$hour.inc == -20] <- "04:00"
sp$date2[sp$hour.inc == -21] <- "03:00"
sp$date2[sp$hour.inc == -22] <- "02:00"
sp$date2[sp$hour.inc == -23] <- "01:00"

sp2 <- sp[!duplicated(sp$conc), ]

epp$date2[epp$hour.inc == 0] <- "24:00"
epp$date2[epp$hour.inc == -1] <- "23:00"
epp$date2[epp$hour.inc == -2] <- "22:00"
epp$date2[epp$hour.inc == -3] <- "21:00"
epp$date2[epp$hour.inc == -4] <- "20:00"
epp$date2[epp$hour.inc == -5] <- "19:00"
epp$date2[epp$hour.inc == -6] <- "18:00"
epp$date2[epp$hour.inc == -7] <- "17:00"
epp$date2[epp$hour.inc == -8] <- "16:00"
epp$date2[epp$hour.inc == -9] <- "15:00"
epp$date2[epp$hour.inc == -10] <- "14:00"
epp$date2[epp$hour.inc == -11] <- "13:00"
epp$date2[epp$hour.inc == -12] <- "12:00"
epp$date2[epp$hour.inc == -13] <- "11:00"
epp$date2[epp$hour.inc == -14] <- "10:00"
epp$date2[epp$hour.inc == -15] <- "09:00"
epp$date2[epp$hour.inc == -16] <- "08:00"
epp$date2[epp$hour.inc == -17] <- "07:00"
epp$date2[epp$hour.inc == -18] <- "06:00"
epp$date2[epp$hour.inc == -19] <- "05:00"
epp$date2[epp$hour.inc == -20] <- "04:00"
epp$date2[epp$hour.inc == -21] <- "03:00"
epp$date2[epp$hour.inc == -22] <- "02:00"
epp$date2[epp$hour.inc == -23] <- "01:00"

epp2 <- epp[!duplicated(epp$conc), ]

#all_data3 <- all_data2[(all_data2$height < 500),]

all_data3 <- all_data2
all_data3$height <- as.numeric(as.character(all_data3$height))/1000
all_data3 <- subset(all_data3,all_data3$height < 14.5)
all_data3$rec_lat <-  37.5821
all_data3$rec_lon <-  126.9706
all_data3$dist<- distVincentyEllipsoid(subset(all_data3, select = c(lon,lat) ), subset(all_data3, select = c(rec_lon,rec_lat)), a=6378137, b=6356752.3142, f=1/298.257223563)

ggplot(data=all_data3, aes(date2, conc)) + geom_boxplot(size = 2,outlier.shape = NA)+
  coord_fixed(ratio = 0.75) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 500))+
  scale_x_discrete(breaks=c("01:00","05:00","09:00","13:00","17:00","21:00","24:00"),
                   labels=c("10:00AM","2:00PM","6:00PM","10:00PM","2:00AM","6:00AM","9:00AM"))+
  theme_bw(base_size = 35) +
  theme(axis.text.x = element_text(face="bold", color="#993333",
                                   size=35, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",
                                   size=35, angle=45),
        legend.key.size = unit(0.5, 'lines'))+
  xlab("Seoul Local Time") + ylab("CO (ppbv)")
ggsave("rev_lower14500_all.png", width = 80, height = 20, units = "cm")


ggplot(data=all_data3, aes(date2, height)) + geom_boxplot(size = 2,outlier.shape = NA)+
  coord_fixed(ratio = 0.75) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 14.5))+
  scale_x_discrete(breaks=c("01:00","05:00","09:00","13:00","17:00","21:00","24:00"),
                   labels=c("10:00AM","2:00PM","6:00PM","10:00PM","2:00AM","6:00AM","9:00AM"))+
  theme_bw(base_size = 35) +
  theme(axis.text.x = element_text(face="bold", color="#993333",
                                   size=35, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",
                                   size=35, angle=45),
    legend.key.size = unit(0.5, 'lines'))+
  xlab("Seoul Local Time") + ylab("Height (km)")
ggsave("height_test_epp.png", width = 80, height = 20, units = "cm")

# ggplot(data=all_data3, aes(date2, dist/1000)) + geom_boxplot(size = 2,outlier.shape = NA)+
#   #coord_fixed(ratio = 100) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 750))+
#   scale_x_discrete(breaks=c("01:00","05:00","09:00","13:00","17:00","21:00","24:00"),
#                    labels=c("10:00AM","2:00PM","6:00PM","10:00PM","2:00AM","6:00AM","9:00AM"))+
#   theme_bw(base_size = 35) +
#   theme(axis.text.x = element_text(face="bold", color="#993333", 
#                                    size=35, angle=0),
#         axis.text.y = element_text(face="bold", color="#993333", 
#                                    size=35, angle=45),
#         legend.key.size = unit(0.5, 'lines'))+
#   xlab("Seoul Local Time") + ylab("Distance (km)")
# 
# ggsave("dist_lower14500_sp.png", width = 80, height = 20, units = "cm")

# ggplot(data=all_data3, aes(date2, height)) + geom_boxplot(size = 2,outlier.shape = NA)+
#   #coord_fixed(ratio = 100) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 15))+
#   scale_x_discrete(breaks=c("01:00","05:00","09:00","13:00","17:00","21:00","24:00"),
#                    labels=c("10:00AM","2:00PM","6:00PM","10:00PM","2:00AM","6:00AM","9:00AM"))+
#   theme_bw(base_size = 35) +
#   theme(axis.text.x = element_text(face="bold", color="#993333", 
#                                    size=35, angle=0),
#         axis.text.y = element_text(face="bold", color="#993333", 
#                                    size=35, angle=45),
#         legend.key.size = unit(0.5, 'lines'))+
#   xlab("Seoul Local Time") + ylab("Height (km)")
# ggsave("height_all.png", width = 80, height = 20, units = "cm")

height4 <- c(25,25,25,26,26,27,27,27,27,30,30,30,28,28,28,28,28,28,25,26,25,25,25,20)*0.001
date4 <- factor(c("24:00","23:00","22:00","21:00","20:00","19:00","18:00","17:00","16:00","15:00"
           ,"14:00","13:00","12:00","11:00","10:00","9:00","8:00","7:00","6:00","5:00"
           ,"4:00","3:00","2:00","1:00"))
pack4 <- factor(c("First","First","First","First","First","First","First","First",
           "First","First","First","First","First","First","First","First",
           "First","First","First","First","First","First","First","First"))
all4 <-data.frame(pack4,date4,height4)
height4 <- c(45,45,45,45,55,55,55,55,55,55,55,55,55,55,45,45,45,45,35,35,35,35,35,35)*0.001
date4 <- factor(c("24:00","23:00","22:00","21:00","20:00","19:00","18:00","17:00","16:00","15:00"
           ,"14:00","13:00","12:00","11:00","10:00","9:00","8:00","7:00","6:00","5:00"
           ,"4:00","3:00","2:00","1:00"))
pack4 <- factor(c("Second","Second","Second","Second","Second","Second","Second","Second",
           "Second","Second","Second","Second","Second","Second","Second","Second",
           "Second","Second","Second","Second","Second","Second","Second","Second"))
all5 <-data.frame(pack4,date4,height4)
height4 <- c(8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,7.9,7.8,7.7,7.6,7.5,7.4)
date4 <- factor(c("24:00","23:00","22:00","21:00","20:00","19:00","18:00","17:00","16:00","15:00"
           ,"14:00","13:00","12:00","11:00","10:00","9:00","8:00","7:00","6:00","5:00"
           ,"4:00","3:00","2:00","1:00"))
pack4 <- factor(c("Third","Third","Third","Third","Third","Third","Third","Third",
           "Third","Third","Third","Third","Third","Third","Third","Third",
           "Third","Third","Third","Third","Third","Third","Third","Third"))
all6 <-data.frame(pack4,date4,height4)
height4 <- c(11,11,11,11,10.8,10.8,10.5,10.5,10,10,10,10,9.8,9.8,9.8,9.8,9.8,9.8,9.7,9.6,9.5,9.4,9.3,9.2)
date4 <- factor(c("24:00","23:00","22:00","21:00","20:00","19:00","18:00","17:00","16:00","15:00"
           ,"14:00","13:00","12:00","11:00","10:00","9:00","8:00","7:00","6:00","5:00"
           ,"4:00","3:00","2:00","1:00"))
pack4 <- factor(c("Fourth","Fourth","Fourth","Fourth","Fourth","Fourth","Fourth","Fourth",
           "Fourth","Fourth","Fourth","Fourth","Fourth","Fourth","Fourth","Fourth",
           "Fourth","Fourth","Fourth","Fourth","Fourth","Fourth","Fourth","Fourth"))
all7 <-data.frame(pack4,date4,height4)

total <- rbind(all4,all5,all6,all7)

ggplot(data=total, aes(x=date4, y=height4, group=pack4, color=pack4)) + 
  geom_line(size=2.0) +
theme_bw(base_size = 35) +
  scale_x_discrete(breaks=c("01:00","05:00","09:00","13:00","17:00","21:00","24:00"),
                   labels=c("10:00AM","2:00PM","6:00PM","10:00PM","2:00AM","6:00AM","9:00AM"))+
  theme(axis.text.x = element_text(face="bold", color="#993333",
                                   size=35, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",
                                   size=35, angle=45))+
  scale_color_manual(values = c("darkblue","blue","grey","green")) +
  theme(legend.position = "none") +
  xlab("Seoul Local Time") + ylab("Height (km)")
