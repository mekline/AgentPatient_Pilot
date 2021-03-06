#This file reads in ALL the %-signal-change values, per-participant, per-parcel, per-contrast,
# Those %-signal-change calculations are produced by the awesome toolbox analyses, and represent a single overall calculation
#derived for the whole parcel region (not individual voxels, as mk sometimes forgets)

rm(list = ls())
library(bootstrap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

##SET YOUR DIRECTORY!!!

########
#READ IN DATA
########
#Here, we read in all those files, calculate a whole passle of mean and standard error bars, and then make graphs

myResults = read.csv('loc_langloc_crit_agpat_20171230_no_430.csv', 
                     colClasses = c("factor","factor", "factor","numeric","numeric")) %>%
  mutate(Localizer = 'langlocSN') %>%
  mutate(ROISystem = 'LHLang') %>%
  mutate(TaskCrit = 'AgPat')


allSigChange = read.csv('loc_langloc_crit_langloc_20171230_no_430.csv', 
                        colClasses = c("factor","factor", "factor","numeric","numeric")) %>%
  mutate(Localizer = 'langlocSN') %>%
  mutate(ROISystem = 'LHLang') %>%
  mutate(TaskCrit = 'langlocSN')

allSigChange = rbind(allSigChange, myResults)



#########
# TRANSFORMATIONS
#########

avgSigChange = allSigChange %>%
  group_by(Subject,Contrast,ROISystem,Localizer,TaskCrit) %>%
  summarize('SignalChange' = mean(SignalChange)) %>%
  mutate(ROIName = 'Localizer Average') %>%
  mutate(nVoxels = 0)

allSigChange <- union(allSigChange, avgSigChange)
#(Gets mad that there's a new ROIName level, re-factorize :p)
allSigChange$ROIName <- as.factor(allSigChange$ROIName)

#Next, get the table that we'll be making the graphs from: for each region (including the average region), take all 
#the individual signal changes and calculate a mean, a standard error (incase we want it) 
#and bootstrapped CIs (which we'll actually use)

sterr <- function(mylist){
  my_se = sd(mylist)/sqrt(length(mylist)) 
  
  return(my_se)
}

#bootstrapped 95% confidence intervals! calculate them from allSigChange
#then merge into mystats
bootup <- function(mylist){
  foo <- bootstrap(mylist, 1000, mean)
  return(quantile(foo$thetastar, 0.975)[1])
}
bootdown <- function(mylist){
  foo <- bootstrap(mylist, 1000, mean)
  return(quantile(foo$thetastar, 0.025)[1])
}

toGraph <- allSigChange %>%
  group_by(ROIName, Contrast, Localizer, ROISystem, TaskCrit) %>%
  summarize(meanSig = mean(SignalChange), sterr = sterr(SignalChange),
            bootup = bootup(SignalChange), bootdown = bootdown(SignalChange))%>%
  mutate(ROIGroup = ifelse(ROIName == 'Localizer Average','across ROIs','individual ROIs'))


#Force some factor orderings here, they are finicky!  
toGraph$ROIName <- factor(toGraph$ROIName, levels = c("LIFGorb",
                                                      "LIFG",
                                                      "LMFG",
                                                      "LAntTemp",
                                                      "LPostTemp",
                                                      "LAngG",
                                                      "Localizer Average"))

toGraph$ROIGroup <- factor(toGraph$ROIGroup, levels = c("across ROIs",
                                                        "individual ROIs"))

toGraph <- toGraph %>%
  arrange(ROIGroup, ROIName)

#########
# Effect size reports
#########
#report  a simple measure of effect size: the
#mean signal change in each system.  (DO LATER)
            
#########
# Graphs!
#########
#Now we can use the information stored in mystats to make pretty graphs! This could be done in excel too by printing mystats
#Change to figs output folder

figdir = paste(getwd(),'figs')
setwd(figdir)

#Subset and rename for language localiser
LangLoc_LangCrit <- filter(toGraph, Localizer =='langlocSN',TaskCrit =='langlocSN') %>%
  filter(Contrast %in% c('S','N'))
LangLoc_APCrit <- filter(toGraph, Localizer =='langlocSN',TaskCrit =='AgPat') %>%
  filter(Contrast %in% c("agt","pat"))

#More factor reordering
LangLoc_LangCrit$Contrast <- factor(LangLoc_LangCrit$Contrast, levels = c("S","N"))
LangLoc_LangCrit <- LangLoc_LangCrit %>%
  arrange(ROIGroup, ROIName, Contrast)

LangLoc_APCrit$Contrast <- factor(LangLoc_APCrit$Contrast, levels = c("agt","pat"))
LangLoc_APCrit <- LangLoc_APCrit %>%
  arrange(ROIGroup, ROIName, Contrast)

#Graphing function!

makeBar = function(plotData, fileName = 'TEST NAME', ylow=-0.5,yhigh=2.5, mycolors = c("gray35", "gray60")) {
  
  #freeze factor orders, AGAIN
  plotData$ROIName <- factor(plotData$ROIName, levels = unique(plotData$ROIName))
  plotData$ROIGroup <- factor(plotData$ROIGroup, levels = unique(plotData$ROIGroup))
  myfi = paste(fileName, '.jpg', sep="")#filename
  print(myfi)
  
  ggplot(data=plotData, aes(x=ROIName, y=meanSig, fill=Contrast)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=bootdown, ymax=bootup), colour="black", width=.1, position=position_dodge(.9)) +
    coord_cartesian(ylim=c(ylow-0.5,yhigh+0.5)) +
    scale_y_continuous(breaks = seq(ylow-0.5, yhigh+0.5, 0.5))+
    xlab('') +
    ylab(str_wrap('% signal change over fixation', width=18)) +
    theme_bw() +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size = 25)) +
    facet_grid(~ROIGroup, scale='free_x', space='free_x') +
    theme(strip.background = element_blank()) +
    theme(strip.text = element_blank()) 
  
  
  ggsave(filename=myfi, width=length(unique(plotData$ROIName))*2.2, height=6.1)
  
}

makeBar(LangLoc_LangCrit, 'LangLoc_LangCrit_no_430')
makeBar(LangLoc_APCrit, 'LangLoc_APCrit_no_430')

#Old contrast lists

# Add in the contrast and ROI names so it's not just numbers!!!!! (This ordering comes from the 
# standard ordering produced by the 2nd level analyses; we'll arrange differently in the plots)
# 
# RHLangROI.Names = c('RPost Temp', 'RAnt Temp', 'RAngG', 'RIFG',      'RMFG',     'RIFG orb');
# LangROI.Names = c('LPost Temp', 'LAnt Temp', 'LAngG', 'LIFG',      'LMFG',     'LIFG orb');
# 
# MDROI.Names = c('LIFG op',  'RIFG op', 'LMFG',    'RMFG',    'LMFG orb',
#                 'RMFG orb', 'LPrecG', 'RPrecG',  'LInsula', 'RInsula',
#                 'LSMA',    'RSMA',   'LPar Inf', 'RPar Inf', 'LPar Sup',
#                 'RPar Sup', 'LACC',   'RACC');
# 
# ToMROI.Names = c('DM PFC', 'LTPJ',  'MM PFC', 'PC',
#                  'RTPJ',  'VM PFC', 'RSTS');
# 
# normal.contrasts = c('agt', 'pat')
