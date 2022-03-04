rm(list = ls())#####################################################
if (!require("devtools")) {#########################################
  install.packages("devtools", dependencies = TRUE)#################
}###################################################################
library(devtools)###################################################
if (!require("QuickEnvironment")) {#################################
  devtools::install_github("DejanDraschkow/QuickEnvironment", ######
                           force=T)################################# 
}###################################################################
library(QuickEnvironment)###########################################
InstallUsefulPackages()#############################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#########

set.seed(22) # for replicability, but inferences don't change if seed is different

#### data import ###############################################################
# list files
files <- list.files(pattern = ".csv", recursive =F)
files
# load the data to three i_dats, one for each experiment
i_dat2 <- NULL
i_dat3 <- NULL
i_dat4 <- NULL
for(file in files){
  print(file)
  jj <- data.table::fread(file, sep = ",", header = T)
  if ((grepl("Loco_MEMexperimenttracking",file))==TRUE) {
    i_dat2 <- rbind(i_dat2, jj)
  }
  if ((grepl("Real_1_MEMexperimenttracking",file))==TRUE) {
    i_dat3 <- rbind(i_dat3, jj)
  }
  if ((grepl("Real_MEMexperimenttracking",file))==TRUE) {
    i_dat4 <- rbind(i_dat4, jj)
  }
  jj<-NULL
}
i_dat2$ExpNr<-2
i_dat3$ExpNr<-3
i_dat3$SubNr<-(i_dat3$SubNr)+24
i_dat4$ExpNr<-4
i_dat4$SubNr<-(i_dat4$SubNr)+48
i_dat <- rbind(i_dat2, i_dat3, i_dat4)
i_dat2 <- NULL
i_dat3 <- NULL
i_dat4 <- NULL

#### Prep time window and exclusion of trials ##################################################

## set values
timepersample = 1000/90 # 1 sample = 11.1111ms
baseline = 19 # samples, i.e. 200ms
baseline_excl = 46 # 500ms
posttrigger = 180 # samples,i.e. 2s
smoothwindow = 4 # samples,i.e. 44ms

##set phase to same since we are not comparing between experiments
i_dat$phase<- "same"

## exclude trials with head or eye samples larger then half the distance to the bars (in meters)
excl <- i_dat[(abs(i_dat$eyeX)>0.5 | abs(i_dat$headBiasX)>0.5 |abs(i_dat$headRoll)>20|abs(i_dat$headYaw)>20) & i_dat$cueFrame>-baseline_excl & i_dat$cueFrame<posttrigger,]
excl <- aggregate(Time~SubNr+item+block+trial+phase, mean, data=excl)
# erase time signiture since it is not useful
excl$Time <- NULL
# add a exclusion tag variable of 1 to all excluded data rows
excl$excl <- 1
# find length (amount) of each subset of excluded data
aggregate(excl~SubNr+item+phase, length, data=excl)
# find number of excluded trial rows per subject
avgexcl<-aggregate(excl~SubNr, sum, data=aggregate(excl~SubNr+item+phase, length, data=excl))
# add extra column for each subject with the proportion of excluded trial rows
avgexcl$prop <- avgexcl$excl/300

avg_excl_per_participant <- (mean(avgexcl$prop))*100
print(paste("Avg excluded trials per participant:",(mean(avgexcl$prop))*100,"%",sep=" "))
summarySE(data=avgexcl,measurevar = "prop")
# combine the excl column variable to i_dat and the column is NA for good data and has a value for excluded data
dat <- full_join(i_dat,excl, by = c("SubNr","block", "trial", "item","phase"))
# make dat only contain non excluded data
dat <- dat[is.na(dat$excl),]
total_excl_trials <- (100-(nrow(dat)/nrow(i_dat))*100)
print(paste("Overall excluded data:",(100-(nrow(dat)/nrow(i_dat))*100),"%",sep=" "))

#### Prep for head direction (called headBiasX in data) position and towardness Figures ######################

# convert data from meters to cm
dat$headBiasX <- dat$headBiasX*100
# smooth data with moving average
dat$headBiasX <- SMA(dat$headBiasX,smoothwindow)

# baseline the headbiasX values
dat_headBiasX_baseline <- aggregate(headBiasX~SubNr+block+trial, data=dat[(dat$cueFrame>-baseline & dat$cueFrame<=0)], mean)
dat_with_baseline <- full_join(dat,dat_headBiasX_baseline, by =c("SubNr","trial","block"))
dat_with_baseline$headBiasX <- dat_with_baseline$headBiasX.x-dat_with_baseline$headBiasX.y
dat <- dat_with_baseline

# aggregate for cue period
jj <- data.table::setDT(dat[dat$cueFrame>-baseline & dat$cueFrame<posttrigger,])[, lapply(.SD, mean), 
                                                                                 by=.(SubNr, block,item,cueFrame,phase), 
                                                                                 .SDcols=c("headBiasX")]
# aggregate across blocks
agg <- data.table::setDF( 
  data.table::setDT(jj)[, lapply(.SD, mean), 
                        by=.(SubNr,item,cueFrame,phase), 
                        .SDcols=c("headBiasX")] 
)

# create towardness values for sensitive comparisons
towardR <- agg[agg$item == "item R",]
towardL <- agg[agg$item == "item L",]
toward <- full_join(towardR,towardL, by =c("SubNr","cueFrame"))
# add left and right  movement and half to find average towardness irrespective of direction
toward$headBiasX <- (toward$headBiasX.x-toward$headBiasX.y)/2
toward$phase <- toward$phase.x
toward$headBiasX.x <- NULL
toward$headBiasX.y <- NULL

#### head direction position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "headBiasX", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(head_direction_position_plot<-ggplot(table, aes(x = time, y = headBiasX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=headBiasX+se, ymin=headBiasX-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-1,1)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))



#### head direction towardness plot #####################################
tableHead<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "headBiasX", 
                 groupvars = c("cueFrame"), na.rm = TRUE)
tableHead$time <- as.numeric(as.character(tableHead$cueFrame))*timepersample

(head_direction_towardness_plot<-ggplot(tableHead, aes(x = time, y = headBiasX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=headBiasX+se, ymin=headBiasX-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-0.2,0.75),breaks=c(0,0.5))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

## stats for head direction towardness plot ##
# change table back to tableHead for stats
table <- tableHead
# change shape of data frame to fit the cluster-based permutation test
headBiasX <- tidyr::spread(toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,]
                           [c("SubNr","cueFrame","headBiasX","phase")], 
                           cueFrame, headBiasX)

# create 0 data for comparison from a normal distribution with mean 0 and sd = paticipants sd
randsample <- NULL
for (i in 1:unique(table$N)){
  r<-rnorm(n = length(table$sd), mean = 0, sd = mean(table$sd))
  randsample <- c(randsample,r)
}
head(randsample)
randsample <- data.frame(matrix(randsample,nrow = unique(table$N), ncol = length(table$sd),byrow = T))
head(randsample[,1:5])
rand <- data.frame(SubNr=1:unique(table$N), phase = "random")
rand <- cbind(rand,randsample)
colnames(rand)<-colnames(headBiasX)
headBiasX <- rbind(headBiasX,rand)

# perform cbpt

(Clust <-permuco::clusterlm(headBiasX[,(3:ncol(headBiasX))] ~ phase + Error(SubNr/phase), 
                            data = headBiasX[,-(3:ncol(headBiasX))],multcomp ="clusterdepth", 
                            np = 5000))
plot(Clust, nbbaselinepts = 65, nbptsperunit = .90)
summary(Clust)

sig_head <- (which(summary(Clust,table_type = "full")$phase$`P(>) clusterdepth`<0.05)-baseline)*timepersample

(head_direction_towardness_plot_stats<-head_direction_towardness_plot + annotate("point", x = sig_head ,  y =-0.1, colour = "black", size=2))


#### head direction density map ##############################

# create density maps with raw, unfiltered, full data
hb_heat <- i_dat[i_dat$cueFrame>45 & i_dat$cueFrame<180]

## convert data from numbers to cm
hb_heat$headBiasX <- hb_heat$headBiasX*100
hb_heat$headBiasY <- hb_heat$headBiasY*100
hb_heat$DisplayHeight <- hb_heat$DisplayHeight*100

# correct Y axis for participant height == display height
hb_heat$headBiasY <- hb_heat$headBiasY-hb_heat$DisplayHeight

hb_heat_table <- data.table::setDF(hb_heat)
itemR <- hb_heat_table[c("headBiasX","headBiasY")][hb_heat_table$item == "item R",]
itemL <- hb_heat_table[c("headBiasX","headBiasY")][hb_heat_table$item == "item L",]

# plot density map
plot(density(itemR$headBiasX ,
             from=-200, to=200))
plot(density(itemL$headBiasX,
             from=-200, to=200))
plot(density(itemR$headBiasX, from=-200, to=200)$y-
       density(itemL$headBiasX, from=-200, to=200)$y)
n <- 200
densR <- MASS::kde2d(itemR$headBiasX, itemR$headBiasY, h = c(40, 40),
                     lims = c(-200,200,-200,200),n=n)
densL <- MASS::kde2d(itemL$headBiasX, itemL$headBiasY, h = c(40, 40),
                     lims = c(-200,200,-200,200),n=n)
dens <- densR$z-densL$z
dens<-matrix(dens,n,n)
sdens <- stack(as.data.frame(dens))
sdens$x <- rep(seq_len(nrow(dens)), ncol(dens))-(n/2)
sdens$y <- as.numeric(sdens$ind)-(n/2)
sdens$density <- sdens$values

ggplot(sdens, aes(x, y)) +
  geom_raster(aes(fill = density)) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",midpoint = 0) +
  theme_set(theme_gray(base_size = 60))

(head_direction_densitymap<-ggplot(sdens, aes(x, y)) +
    geom_raster(aes(fill = density)) + 
    scale_fill_gradient2(low="blue",mid="white", high="red",midpoint = 0) +
    theme_set(theme_gray(base_size = 60))+
    scale_y_continuous(breaks=c(-100,0, 100),limits = NULL)+   
    scale_x_continuous(breaks=c(-100,0, 100),limits = NULL)+ 
    coord_cartesian(xlim=c(-125, 125),ylim=c(-125, 125))  +
    annotate("segment", x = -6, xend = 6, y = 0, yend = 0,
             colour = "black", size=1,alpha=1)+
    annotate("segment", x = 0, xend = 0, y = -6, yend = 6,
             colour = "black", size=1,alpha=1)+
    annotate("point", size=1,x=-100,y=0)+
    annotate("path", size=0.5,
             x=-100+25*cos(seq(0,2*pi,length.out=100)),
             y=0+25*sin(seq(0,2*pi,length.out=100))) +
    annotate("point", size=1,x=100,y=0)+
    annotate("path", size=0.5,
             x=100+25*cos(seq(0,2*pi,length.out=100)),
             y=0+25*sin(seq(0,2*pi,length.out=100))) +
    theme(plot.title = element_text(hjust = 0.5),         
          axis.title.x=element_blank(),          
          axis.title.y=element_blank(),
          axis.text = element_text(size = 18),
          legend.title=element_blank(),
          legend.text=element_blank(),
          legend.position= c(1, 0.9),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )


#### eye direction data prep ####

# convert eye direction data (called eyeX in data) from meters to cm 
dat$eyeX <- dat$eyeX*100
# smooth data with moving average
dat$eyeX <- SMA(dat$eyeX,smoothwindow)

# baseline the eyeX values
dat_eyeX_baseline <- aggregate(eyeX~SubNr+block+trial, data=dat[(dat$cueFrame>-baseline & dat$cueFrame<=0)], mean)
dat_with_baseline <- full_join(dat,dat_eyeX_baseline, by =c("SubNr","trial","block"))
dat_with_baseline$eyeX <- dat_with_baseline$eyeX.x-dat_with_baseline$eyeX.y
dat <- dat_with_baseline

## aggregate for cue period
jj <- data.table::setDT(dat[dat$cueFrame>-baseline & dat$cueFrame<posttrigger,])[, lapply(.SD, mean), 
                                                                                 by=.(SubNr, block,item,cueFrame,phase), 
                                                                                 .SDcols=c("eyeX")]
# aggregate across blocks
agg <- data.table::setDF( 
  data.table::setDT(jj)[, lapply(.SD, mean), 
                        by=.(SubNr,item,cueFrame,phase), 
                        .SDcols=c("eyeX")] 
)

## create towardness values for sensitive comparisons
towardR <- agg[agg$item == "item R",]
towardL <- agg[agg$item == "item L",]
toward <- full_join(towardR,towardL, by =c("SubNr","cueFrame"))
# add left and right  movement and half to find average towardness irrespective of direction
toward$eyeX <- (toward$eyeX.x-toward$eyeX.y)/2
toward$phase <- toward$phase.x
toward$eyeX.x <- NULL
toward$eyeX.y <- NULL

#### eye direction position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "eyeX", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(eye_direction_position_plot<-ggplot(table, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-2,2.2)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-2,-1,0,1,2))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))

#### eye direction overlayed with head direction towardness plot) #####################################
table<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "eyeX", 
                    groupvars = c("cueFrame"), na.rm = TRUE)
table$time <- as.numeric(as.character(table$cueFrame))*timepersample

tableEye<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "eyeX", 
                 groupvars = c("cueFrame"), na.rm = TRUE)
tableEye$time <- as.numeric(as.character(tableEye$cueFrame))*timepersample

tableEye$type <- "eye"
tableHead$type <- "head"
colnames(tableEye) <- c("cueFrame","N", "Bias","sd", "se", "ci", "time", "type" )
colnames(tableHead) <- c("cueFrame","N", "Bias","sd", "se", "ci", "time", "type" )

tableEyeHead <- rbind(tableHead, tableEye)

tableEyeHead$type <- factor(tableEyeHead$type)
levels(tableEyeHead$type)

(toward_eye_N<-ggplot(tableEyeHead, aes(x = time, y = Bias)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(type)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=Bias+se, ymin=Bias-se, fill=factor(type)), alpha = 0.2) +
    coord_cartesian(ylim=c(-0.5, 2))   + 
    theme_set(theme_gray(base_size = 48))+ 
    scale_color_manual(name = "",
                       labels=c("gaze", "horizontal head direction"),
                       values = c("grey68", "black")) + 
    scale_fill_manual(name = "",
                      values = c("grey68", "black")) + guides(fill=FALSE) +
    scale_y_continuous(breaks=c(-1,0,1,2,3))+   
    scale_x_continuous(breaks=c(0, 1000, 2000))+  
    theme(legend.text = element_blank(),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          legend.position= "none", 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

## stats for towardness plot ##
# change column names of tableEye back to be specific to eyeX

# change shape of data frame to fit the cluster-based permutation test
eyeX <- tidyr::spread(toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,]
                       [c("SubNr","cueFrame","eyeX","phase")], 
                       cueFrame, eyeX)

# create 0 data for comparison from a normal distribution with mean 0 and sd = paticipants sd
randsample <- NULL
for (i in 1:unique(table$N)){
  r<-rnorm(n = length(table$sd), mean = 0, sd = mean(table$sd))
  randsample <- c(randsample,r)
}
head(randsample)
randsample <- data.frame(matrix(randsample,nrow = unique(table$N), ncol = length(table$sd),byrow = T))
head(randsample[,1:5])
rand <- data.frame(SubNr=1:unique(table$N), phase = "random")
rand <- cbind(rand,randsample)
colnames(rand)<-colnames(eyeX)
eyeX <- rbind(eyeX,rand)

# perform cbpt
(Clust <-permuco::clusterlm(eyeX[,(3:ncol(eyeX))] ~ phase + Error(SubNr/phase), 
                            data = eyeX[,-(3:ncol(eyeX))],multcomp ="clusterdepth", 
                            np = 5000))
plot(Clust, nbbaselinepts = 65, nbptsperunit = .90)
summary(Clust)

sig_eye <- (which(summary(Clust,table_type = "full")$phase$`P(>) clusterdepth`<0.05)-baseline)*timepersample

(eye_head_towardness_plot_stats<- toward_eye_N + annotate("point", x = sig_eye ,  y =-0.4, colour = "grey68", size=2) + annotate("point", x = sig_head ,  y =-0.2, colour = "black", size=2))


