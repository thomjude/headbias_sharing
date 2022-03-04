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

#### head translation data prep ####

# convert head translation data (called headX in data) from meters to cm
dat$headX <- dat$headX*100
# smooth data with moving average
dat$headX <- SMA(dat$headX,smoothwindow)

# baseline the headX values
dat_headX_baseline <- aggregate(headX~SubNr+block+trial, data=dat[(dat$cueFrame>-baseline & dat$cueFrame<=0)], mean)
dat_with_baseline <- full_join(dat,dat_headX_baseline, by =c("SubNr","trial","block"))
dat_with_baseline$headX <- dat_with_baseline$headX.x-dat_with_baseline$headX.y
dat <- dat_with_baseline

## aggregate for cue period
jj <- data.table::setDT(dat[dat$cueFrame>-baseline & dat$cueFrame<posttrigger,])[, lapply(.SD, mean), 
                                                                                 by=.(SubNr, block,item,cueFrame,phase), 
                                                                                 .SDcols=c("headX")]
# aggregate across blocks
agg <- data.table::setDF( 
  data.table::setDT(jj)[, lapply(.SD, mean), 
                        by=.(SubNr,item,cueFrame,phase), 
                        .SDcols=c("headX")] 
)

## create towardness values for sensitive comparisons
towardR <- agg[agg$item == "item R",]
towardL <- agg[agg$item == "item L",]
toward <- full_join(towardR,towardL, by =c("SubNr","cueFrame"))
# add left and right  movement and half to find average towardness irrespective of direction
toward$headX <- (toward$headX.x-toward$headX.y)/2
toward$phase <- toward$phase.x
toward$headX.x <- NULL
toward$headX.y <- NULL

#### translation sideways position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "headX", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(translation_position_plot<-ggplot(table, aes(x = time, y = headX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=headX+se, ymin=headX-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-0.7,0.7)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-0.5,0,0.5))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))

#### translation towardness plot #####################################
table<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "headX", 
                 groupvars = c("cueFrame"), na.rm = TRUE)
table$time <- as.numeric(as.character(table$cueFrame))*timepersample

(translation_towardness_plot<-ggplot(table, aes(x = time, y = headX)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=headX+se, ymin=headX-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-0.1,0.5),breaks=c(0,0.25,0.5))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

## stats for towardness plot ##
# change shape of data frame to fit the cluster-based permutation test
headX <- tidyr::spread(toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,]
                       [c("SubNr","cueFrame","headX","phase")], 
                       cueFrame, headX)

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
colnames(rand)<-colnames(headX)
headX <- rbind(headX,rand)

# perform cbpt
(Clust <-permuco::clusterlm(headX[,(3:ncol(headX))] ~ phase + Error(SubNr/phase), 
                            data = headX[,-(3:ncol(headX))],multcomp ="clusterdepth", 
                            np = 5000))
plot(Clust, nbbaselinepts = 65, nbptsperunit = .90)
summary(Clust)

sig_translation <- (which(summary(Clust,table_type = "full")$phase$`P(>) clusterdepth`<0.05)-baseline)*timepersample

(translation_towardness_plot_stats <- translation_towardness_plot + annotate("point", x = sig_translation ,  y =-0.05, colour = "black", size=2))

#### roll data prep ####

# reverse roll data (called headRoll in data) to make positive values be rightwards (chin toward left shoulder)
dat$headRoll <- dat$headRoll * -1
# smooth data with moving average
dat$headRoll <- SMA(dat$headRoll,smoothwindow)

# baseline the headRoll values
dat_headRoll_baseline <- aggregate(headRoll~SubNr+block+trial, data=dat[(dat$cueFrame>-baseline & dat$cueFrame<=0)], mean)
dat_with_baseline <- full_join(dat,dat_headRoll_baseline, by =c("SubNr","trial","block"))
dat_with_baseline$headRoll <- dat_with_baseline$headRoll.x-dat_with_baseline$headRoll.y
dat <- dat_with_baseline


## aggregate for cue period
jj <- data.table::setDT(dat[dat$cueFrame>-baseline & dat$cueFrame<posttrigger,])[, lapply(.SD, mean), 
                                                                                 by=.(SubNr, block,item,cueFrame,phase), 
                                                                                 .SDcols=c("headRoll")]
# aggregate across blocks
agg <- data.table::setDF( 
  data.table::setDT(jj)[, lapply(.SD, mean), 
                        by=.(SubNr,item,cueFrame,phase), 
                        .SDcols=c("headRoll")] 
)

## create towardness values for sensitive comparisons
towardR <- agg[agg$item == "item R",]
towardL <- agg[agg$item == "item L",]
toward <- full_join(towardR,towardL, by =c("SubNr","cueFrame"))
# add left and right  movement and half to find average towardness irrespective of direction
toward$headRoll <- (toward$headRoll.x-toward$headRoll.y)/2
toward$phase <- toward$phase.x
toward$headRoll.x <- NULL
toward$headRoll.y <- NULL

#### roll position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger,], measurevar = "headRoll", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(roll_position_plot<-ggplot(table, aes(x = time, y = headRoll)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=headRoll+se, ymin=headRoll-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-0.3,0.3)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-0.25,0,0.25))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))

#### roll towardness plot) #####################################
table<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "headRoll", 
                 groupvars = c("cueFrame"), na.rm = TRUE)
table$time <- as.numeric(as.character(table$cueFrame))*timepersample

(roll_towardness_plot<-ggplot(table, aes(x = time, y = headRoll)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=headRoll+se, ymin=headRoll-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-0.1,0.2),breaks=c(-0.1,0,0.1,0.2))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

## stats for towardness plot ##
# change shape of data frame to fit the cluster-based permutation test
headRoll <- tidyr::spread(toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,]
                          [c("SubNr","cueFrame","headRoll","phase")], 
                          cueFrame, headRoll)

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
colnames(rand)<-colnames(headRoll)
headRoll <- rbind(headRoll,rand)

# perform cbpt
(Clust <-permuco::clusterlm(headRoll[,(3:ncol(headRoll))] ~ phase + Error(SubNr/phase), 
                            data = headRoll[,-(3:ncol(headRoll))],multcomp ="clusterdepth", 
                            np = 5000))
plot(Clust, nbbaselinepts = 65, nbptsperunit = .90)
summary(Clust)

sig_roll <- (which(summary(Clust,table_type = "full")$phase$`P(>) clusterdepth`<0.05)-baseline)*timepersample

(roll_towardness_plot_stats <- roll_towardness_plot + annotate("point", x = sig_roll ,  y =-0.09, colour = "black", size=2))

#### yaw data prep ####

# smooth yaw data (called headYaw in data) with moving average
dat$headYaw <- SMA(dat$headYaw,smoothwindow)

# baseline the headYaw values
dat_headYaw_baseline <- aggregate(headYaw~SubNr+block+trial, data=dat[(dat$cueFrame>-baseline & dat$cueFrame<=0)], mean)
dat_with_baseline <- full_join(dat,dat_headYaw_baseline, by =c("SubNr","trial","block"))
dat_with_baseline$headYaw <- dat_with_baseline$headYaw.x-dat_with_baseline$headYaw.y
dat <- dat_with_baseline

## aggregate for cue period
jj <- data.table::setDT(dat[dat$cueFrame>-baseline & dat$cueFrame<posttrigger,])[, lapply(.SD, mean), 
                                                                                 by=.(SubNr, block,item,cueFrame,phase), 
                                                                                 .SDcols=c("headYaw")]
# aggregate across blocks
agg <- data.table::setDF( 
  data.table::setDT(jj)[, lapply(.SD, mean), 
                        by=.(SubNr,item,cueFrame,phase), 
                        .SDcols=c("headYaw")] 
)

## create towardness values for sensitive comparisons
towardR <- agg[agg$item == "item R",]
towardL <- agg[agg$item == "item L",]
toward <- full_join(towardR,towardL, by =c("SubNr","cueFrame"))
# add left and right  movement and half to find average towardness irrespective of direction
toward$headYaw <- (toward$headYaw.x-toward$headYaw.y)/2
toward$phase <- toward$phase.x
toward$headYaw.x <- NULL
toward$headYaw.y <- NULL

#### yaw position plot ########################################
table<-summarySE(data=agg[agg$cueFrame>-baseline & agg$cueFrame<posttrigger ,], measurevar = "headYaw", 
                 groupvars = c("item","cueFrame"))
table$time <- as.numeric(as.character(table$cueFrame))*timepersample
(yaw_position_plot<-ggplot(table, aes(x = time, y = headYaw)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(aes(color = factor(item)),alpha = 1,size=2) +
    geom_ribbon(aes(ymax=headYaw+se, ymin=headYaw-se, fill=factor(item)), alpha = 0.2)+
    coord_cartesian(ylim=c(-0.3,0.3)) + 
    theme_set(theme_gray(base_size = 60))+
    scale_color_manual(name = "",
                       labels=c("L item", "R item"),
                       values = c("blue", "red")) + 
    scale_fill_manual(name = "",
                      labels=c("L item", "R item"),
                      values = c("blue", "red")) + guides(fill=FALSE) + guides(colour=FALSE) + 
    scale_y_continuous(breaks=c(-0.25,0,0.25))+   
    scale_x_continuous(breaks=c(0,1000,2000))+   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          legend.position= c(.2, .8), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")))

#### yaw towardness plots #####################################
table<-summarySE(data=toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,],measurevar = "headYaw", 
                 groupvars = c("cueFrame"), na.rm = TRUE)
table$time <- as.numeric(as.character(table$cueFrame))*timepersample

(yaw_towardness_plot<-ggplot(table, aes(x = time, y = headYaw)) + 
    geom_hline(yintercept=0, linetype="dashed",size=1, alpha=1) +
    geom_vline(xintercept=0, linetype="dashed",size=0.2, alpha=0.5) +
    geom_line(size=2) +
    geom_ribbon(aes(ymax=headYaw+se, ymin=headYaw-se), alpha = 0.2)+
    theme_set(theme_gray(base_size = 60))+ 
    scale_y_continuous(limits=c(-0.1,0.2),breaks=c(-0.1,0,0.1,0.2))+ 
    scale_x_continuous(limits=c(-200,2000),breaks=c(0,1000,2000))+ 
    theme(axis.text = element_text(size = 18),
          axis.title.x=element_blank(),           
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) )

## stats for towardness plot ##
# change shape of data frame to fit the cluster-based permutation test
headYaw <- tidyr::spread(toward[toward$cueFrame>-baseline & toward$cueFrame<posttrigger,]
                         [c("SubNr","cueFrame","headYaw","phase")], 
                         cueFrame, headYaw)

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
colnames(rand)<-colnames(headYaw)
headYaw <- rbind(headYaw,rand)

# perform cbpt
(Clust <-permuco::clusterlm(headYaw[,(3:ncol(headYaw))] ~ phase + Error(SubNr/phase), 
                            data = headYaw[,-(3:ncol(headYaw))],multcomp ="clusterdepth", 
                            np = 5000))
plot(Clust, nbbaselinepts = 65, nbptsperunit = .90)
summary(Clust)

sig <- (which(summary(Clust,table_type = "full")$phase$`P(>) clusterdepth`<0.05)-baseline)*timepersample

(yaw_towardness_plot_stats<-yaw_towardness_plot + annotate("point", x = sig ,  y =-0.08, colour = "black", size=2))


