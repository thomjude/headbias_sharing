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
####################################################################
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
i_dat2$ExpNr<-1
i_dat3$ExpNr<-2
i_dat3$SubNr<-(i_dat3$SubNr)+24
i_dat4$ExpNr<-3
i_dat4$SubNr<-(i_dat4$SubNr)+48
i_dat <- rbind(i_dat2,i_dat3,i_dat4)

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
# aggregate finds summary stats of mean of time for all subsets given
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

dat_all<-dat

#### Prep for exp1 eye direction position and towardness figures ####

## set for data from exp 1
dat<-dat_all[dat_all$ExpNr==1]
# convert data from meters to cm and flip sides so left is positive
dat$eyeX <- dat$eyeX*100
# smooth data with moving average
dat$eyeX <- SMA(dat$eyeX,smoothwindow)

# baseline the eyebisX values
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


#### exp1 eye gaze position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "eyeX", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(eye_position_plot_exp1<-ggplot(table, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-3,3)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-3,-1.5,0,1.5,3))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))


#### exp1 eye gaze towardness plot #####################################
tableEye<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "eyeX", 
                     groupvars = c("cueFrame"), na.rm = TRUE)
tableEye$time <- as.numeric(as.character(tableEye$cueFrame))*timepersample

(eye_towardness_plot_exp1<-ggplot(tableEye, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-1,2.5),breaks=c(-1,0,1,2))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

#### Prep for exp 2 eye direction position and towardness figures ####

## set for data from exp 2
dat<-dat_all[dat_all$ExpNr==2]
# convert data from meters to cm and flip sides so left is positive
dat$eyeX <- dat$eyeX*100
# smooth data with moving average
dat$eyeX <- SMA(dat$eyeX,smoothwindow)

# baseline the EyebisX values
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

#### exp2 eye direction position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "eyeX", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(eye_position_plot_exp2<-ggplot(table, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-3,3)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-3,-1.5,0,1.5,3))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))


#### exp2 eye direction towardness plot #####################################
tableEye<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "eyeX", 
                     groupvars = c("cueFrame"), na.rm = TRUE)
tableEye$time <- as.numeric(as.character(tableEye$cueFrame))*timepersample

(eye_towardness_plot_exp2<-ggplot(tableEye, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-1,2.5),breaks=c(-1,0,1,2))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

#### Prep for exp3 eye direction position and towardness figures ####

## set for data from exp 3
dat<-dat_all[dat_all$ExpNr==3]
# convert data from meters to cm and flip sides so left is positive
dat$eyeX <- dat$eyeX*100
# smooth data with moving average
dat$eyeX <- SMA(dat$eyeX,smoothwindow)

# baseline the EyebisX values
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

#### exp3 eye direction position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "eyeX", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(eye_position_plot_exp3<-ggplot(table, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-3,3)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-3,-1.5,0,1.5,3))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))

#### exp3 eye direction towardness plot #####################################
tableEye<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "eyeX", 
                     groupvars = c("cueFrame"), na.rm = TRUE)
tableEye$time <- as.numeric(as.character(tableEye$cueFrame))*timepersample

(eye_towardness_plot_exp3<-ggplot(tableEye, aes(x = time, y = eyeX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=eyeX+se, ymin=eyeX-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-1,2.5),breaks=c(-1,0,1,2))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

