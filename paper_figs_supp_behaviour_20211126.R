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
}
i_dat2$ExpNr<-1
i_dat3$ExpNr<-2
i_dat3$SubNr<-(i_dat3$SubNr)+24
i_dat4$ExpNr<-3
i_dat4$SubNr<-(i_dat4$SubNr)+48
i_dat <- rbind(i_dat2,i_dat3,i_dat4)

#### Prep for tie window and exclusion of trials ##################################################

## set values
timepersample = 1000/90 # 1 sample = 11.1111ms
baseline_excl = 46 # 500ms
baseline = 19 # samples, i.e. 200ms
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

#save other experiment data for later
dat_all<- dat


#### Exp1 error and RT plots and stats ####

#set to only experiment 1
dat <- dat_all[dat_all$ExpNr==1,]

behav_dat <- dat[dat$cueFrame==0,]
jj <- aggregate(cbind(ResponseTime,error)~SubNr+block+item+phase, mean, data=behav_dat)
aggBehav <- aggregate(cbind(ResponseTime,error)~SubNr+item+phase, mean, data=jj)

## errors
(table<-summarySEwithin(data=aggBehav,
                        measurevar = "error", withinvars = c("item"), idvar = "SubNr"))

(error_plot_exp1<-ggplot(table, aes(x = item, y = error)) + ggtitle("") +
    geom_hline(yintercept=45, linetype="dashed",size=1, alpha=0.5) +
    xlab("")+ ylab("") + coord_cartesian(ylim=c(0, 50)) +
    geom_bar(stat="identity",alpha = 0.5,size=1, aes(fill = item), colour="black") +
    geom_errorbar(aes(ymax=error+ci, ymin=error-ci), width=.1,size=1)+
    geom_point(data = aggBehav, aes(x = item, y = error), alpha=0.2 )+
    geom_line(data = aggBehav, aes(x = item, y = error,group=SubNr), alpha=0.2 ) +
    scale_fill_manual(labels=c("L item", "R item"), values=c("blue", "red"))+
    annotate("text", label = "chance error", x = 1.5, y = 50, size = 12, colour = "black", alpha=0.5) +
    theme_set(theme_gray(base_size = 44))+
    scale_y_continuous(breaks=c(0, 20, 45))+   
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          legend.position= "none",
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid")) )
## stats
ANOVA_error_exp1 <- ezANOVA(aggBehav,dv="error",within=c(item),wid=SubNr)

aggDiff <- aggBehav$error[aggBehav$item=="item R"]-aggBehav$error[aggBehav$item=="item L"]
BF_error_exp1 <- BayesFactor::ttestBF(x = aggDiff)


## RTs
(table<-summarySEwithin(data=aggBehav,
                        measurevar = "ResponseTime", withinvars = c("item"), idvar = "SubNr"))

(rt_plot_exp1<-ggplot(table, aes(x = item, y = ResponseTime)) + ggtitle("") +
    xlab("")+ ylab("") + coord_cartesian(ylim=c(0, 2)) +
    geom_bar(stat="identity",alpha = 0.5,size=1, aes(fill = item), colour="black") +
    geom_errorbar(aes(ymax=ResponseTime+ci, ymin=ResponseTime-ci), width=.1,size=1)+
    scale_fill_manual(labels=c("L item", "R item"), values=c("blue", "red"))+
    geom_point(data = aggBehav, aes(x = item, y = ResponseTime), alpha=0.2 )+
    geom_line(data = aggBehav, aes(x = item, y = ResponseTime,group=SubNr), alpha=0.2 ) +
    guides(fill=F)+
    theme_set(theme_gray(base_size = 44))+
    scale_y_continuous(breaks=c(0, 1, 2))+   
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          legend.position= c(.8, .8),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid")) )
## stats
ANOVA_RT_exp1 <- ezANOVA(aggBehav,dv="ResponseTime",within=c(item),wid=SubNr)

aggDiff <- aggBehav$ResponseTime[aggBehav$item=="item R"]-aggBehav$ResponseTime[aggBehav$item=="item L"]
BF_RT_exp1 <- BayesFactor::ttestBF(x = aggDiff)




#### Exp2 error and RT plots and stats ####

#set to only experiment 2
dat <- dat_all[dat_all$ExpNr==2,]

behav_dat <- dat[dat$cueFrame==0,]
jj <- aggregate(cbind(ResponseTime,error)~SubNr+block+item+phase, mean, data=behav_dat)
aggBehav <- aggregate(cbind(ResponseTime,error)~SubNr+item+phase, mean, data=jj)

## errors
(table<-summarySEwithin(data=aggBehav,
                        measurevar = "error", withinvars = c("item"), idvar = "SubNr"))

(error_plot_exp2<-ggplot(table, aes(x = item, y = error)) + ggtitle("") +
    geom_hline(yintercept=45, linetype="dashed",size=1, alpha=0.5) +
    xlab("")+ ylab("") + coord_cartesian(ylim=c(0, 50)) +
    geom_bar(stat="identity",alpha = 0.5,size=1, aes(fill = item), colour="black") +
    geom_errorbar(aes(ymax=error+ci, ymin=error-ci), width=.1,size=1)+
    geom_point(data = aggBehav, aes(x = item, y = error), alpha=0.2 )+
    geom_line(data = aggBehav, aes(x = item, y = error,group=SubNr), alpha=0.2 ) +
    scale_fill_manual(labels=c("L item", "R item"), values=c("blue", "red"))+
    annotate("text", label = "chance error", x = 1.5, y = 50, size = 12, colour = "black", alpha=0.5) +
    theme_set(theme_gray(base_size = 44))+
    scale_y_continuous(breaks=c(0, 20, 45))+   
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          legend.position= "none",
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid")) )
## stats
ANOVA_error_exp2 <- ezANOVA(aggBehav,dv="error",within=c(item),wid=SubNr)

aggDiff <- aggBehav$error[aggBehav$item=="item R"]-aggBehav$error[aggBehav$item=="item L"]
BF_error_exp2 <- BayesFactor::ttestBF(x = aggDiff)


## RTs
(table<-summarySEwithin(data=aggBehav,
                        measurevar = "ResponseTime", withinvars = c("item"), idvar = "SubNr"))

(rt_plot_exp2<-ggplot(table, aes(x = item, y = ResponseTime)) + ggtitle("") +
    xlab("")+ ylab("") + coord_cartesian(ylim=c(0, 2)) +
    geom_bar(stat="identity",alpha = 0.5,size=1, aes(fill = item), colour="black") +
    geom_errorbar(aes(ymax=ResponseTime+ci, ymin=ResponseTime-ci), width=.1,size=1)+
    scale_fill_manual(labels=c("L item", "R item"), values=c("blue", "red"))+
    geom_point(data = aggBehav, aes(x = item, y = ResponseTime), alpha=0.2 )+
    geom_line(data = aggBehav, aes(x = item, y = ResponseTime,group=SubNr), alpha=0.2 ) +
    guides(fill=F)+
    theme_set(theme_gray(base_size = 44))+
    scale_y_continuous(breaks=c(0, 1, 2))+   
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          legend.position= c(.8, .8),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid")) )
## stats
ANOVA_RT_exp2 <- ezANOVA(aggBehav,dv="ResponseTime",within=c(item),wid=SubNr)

aggDiff <- aggBehav$ResponseTime[aggBehav$item=="item R"]-aggBehav$ResponseTime[aggBehav$item=="item L"]
BF_RT_exp2 <- BayesFactor::ttestBF(x = aggDiff)


#### Exp3 error and RT plots and stats ####


#set to only experiment 3
dat <- dat_all[dat_all$ExpNr==3,]

behav_dat <- dat[dat$cueFrame==0,]
jj <- aggregate(cbind(ResponseTime,error)~SubNr+block+item+phase, mean, data=behav_dat)
aggBehav <- aggregate(cbind(ResponseTime,error)~SubNr+item+phase, mean, data=jj)

## errors
(table<-summarySEwithin(data=aggBehav,
                        measurevar = "error", withinvars = c("item"), idvar = "SubNr"))

(error_plot_exp3<-ggplot(table, aes(x = item, y = error)) + ggtitle("") +
    geom_hline(yintercept=45, linetype="dashed",size=1, alpha=0.5) +
    xlab("")+ ylab("") + coord_cartesian(ylim=c(0, 50)) +
    geom_bar(stat="identity",alpha = 0.5,size=1, aes(fill = item), colour="black") +
    geom_errorbar(aes(ymax=error+ci, ymin=error-ci), width=.1,size=1)+
    geom_point(data = aggBehav, aes(x = item, y = error), alpha=0.2 )+
    geom_line(data = aggBehav, aes(x = item, y = error,group=SubNr), alpha=0.2 ) +
    scale_fill_manual(labels=c("L item", "R item"), values=c("blue", "red"))+
    annotate("text", label = "chance error", x = 1.5, y = 50, size = 12, colour = "black", alpha=0.5) +
    theme_set(theme_gray(base_size = 44))+
    scale_y_continuous(breaks=c(0, 20, 45))+   
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          legend.position= "none",
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid")) )
## stats
ANOVA_error_exp3 <- ezANOVA(aggBehav,dv="error",within=c(item),wid=SubNr)

aggDiff <- aggBehav$error[aggBehav$item=="item R"]-aggBehav$error[aggBehav$item=="item L"]
BF_error_exp3 <- BayesFactor::ttestBF(x = aggDiff)


## RTs
(table<-summarySEwithin(data=aggBehav,
                        measurevar = "ResponseTime", withinvars = c("item"), idvar = "SubNr"))

(rt_plot_exp3<-ggplot(table, aes(x = item, y = ResponseTime)) + ggtitle("") +
    xlab("")+ ylab("") + coord_cartesian(ylim=c(0, 2)) +
    geom_bar(stat="identity",alpha = 0.5,size=1, aes(fill = item), colour="black") +
    geom_errorbar(aes(ymax=ResponseTime+ci, ymin=ResponseTime-ci), width=.1,size=1)+
    scale_fill_manual(labels=c("L item", "R item"), values=c("blue", "red"))+
    geom_point(data = aggBehav, aes(x = item, y = ResponseTime), alpha=0.2 )+
    geom_line(data = aggBehav, aes(x = item, y = ResponseTime,group=SubNr), alpha=0.2 ) +
    guides(fill=F)+
    theme_set(theme_gray(base_size = 44))+
    scale_y_continuous(breaks=c(0, 1, 2))+   
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          legend.position= c(.8, .8),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid")) )
## stats
ANOVA_RT_exp3 <- ezANOVA(aggBehav,dv="ResponseTime",within=c(item),wid=SubNr)

aggDiff <- aggBehav$ResponseTime[aggBehav$item=="item R"]-aggBehav$ResponseTime[aggBehav$item=="item L"]
BF_RT_exp3 <- BayesFactor::ttestBF(x = aggDiff)



#### stats for all experiments for print ####
ANOVA_RT_exp1
BF_RT_exp1
ANOVA_error_exp1
BF_error_exp1
ANOVA_RT_exp2
BF_RT_exp2
ANOVA_error_exp2
BF_error_exp2
ANOVA_RT_exp3
BF_RT_exp3
ANOVA_error_exp3
BF_error_exp3



