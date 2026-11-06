#=== === === === === === === ===
# Script started by Rebekah Stiling October 2026
# This script uses B-IBI scores from PSSB to calculate trends over time. Much of the script is recycled from the 04_BIBI-trends_for-SFS-2026.R script I made in my, based on scripts Beth Sosik left behind.
# rstiling@kingcounty.gov
#=== === === === === === === ===

#load relevant libraries
library(tidyverse) #for wrangling
library(janitor) #for cleaning names.


#Read in the data #This is all scores from PSSB from all of time from Boise, Vashon, and Ambient
scores<-list.files("Inputs/Scores_Oct2026", full.names = TRUE) |> #List all files names in the folder
  discard(~ grepl("metadata", .x, ignore.case = TRUE)) |> #drop the metadata files which are just info
    map(read_delim)  |> #this "read_delim()" to each item in the list of file names.
    bind_rows() |>  #this puts all the files together
  clean_names() #clean column names

head(scores)

#create empty folder for graphs and csvs
output_folder <- paste0("output",Sys.Date()) 
dir.create(output_folder)

## prepare unique id for actions taken today
today<-Sys.Date()

# # Specifications:
# # only  2022 and after
# # Projects: Ambient, Vashon, Miller-Walker, Boise 
# # BIBI scores
# 
# #The code runs through trends analysis at multiple spatial scales: REgional, subbasin, WRIA and site level trends.
# 
# ###Raw Data Processing####
# #library(Microsoft365R) #for connecting to SharePoint
# library(readxl) #for reading excel file
# library(tidyverse) #for data wrangling
# 
# 
# library(janitor) #for cleaning names
# library(plyr)
# library(openxlsx)
# library(lubridate)
# library(stringr)
# library(tidyverse)
# library(dplyr)
# getwd()
# 
# atts<-read.xlsx("data/inputs/2012_taxa_attributes.xlsx")
# lookup<-read.csv("data/inputs/ORWA_TaxaTranslator_20250812.csv")##update this to latest translator table.
# BCG_atts<-read.csv("data/inputs/ORWA_Attributes_20250812.csv")#load the BCG taxonomic hierarchy
# 
# #create empty folder for graphs and csvs
# output_folder <- paste0("data/output/BIBItrends_",Sys.Date()) 
# dir.create(output_folder)
# 
# 
# ###This function takes the raw taxa data .txt output from PSSB and binds it all into one data object
# taxaBind <- function(file.path) {
#   
#   path.files <- list.files(file.path)
#   # read in files
#   list.with.each.file <- lapply(paste(file.path, list.files(file.path), sep = ''), function(y) read.delim(y, header=TRUE))
#   taxa<-do.call("rbind.data.frame", list.with.each.file)
#   return(taxa)
#   
#   
# }
# 
# #bind text files into one large dataframe
# file.path="data/pssb_2002-2025_BIBIscores/" #where are the text files stored?
# BIBI_pssb_messy <-taxaBind(file.path)
# write.csv(BIBI_pssb_messy, paste0(output_folder, "/BIBI-scores_from-PSSB.csv"))
BIBI_pssb_messy<-scores

#Preliminary exploration of the data set
BIBI_pssb_messy |> select(project) |> unique() 

BIBI_pssb <- BIBI_pssb_messy |> filter(project != "Macroinvertebrates")

BIBI_pssb |> select(project) |> unique() 

BIBI_pssb |> select(agency) |> unique()
BIBI_pssb |> select(site_code) |> unique()
BIBI_pssb |> select(wria_number) |> unique()
BIBI_pssb |> select(stream) |> unique()
BIBI_pssb |> select(subbasin) |> unique()
BIBI_pssb |> select(visit_id) |> n_distinct()


####Trends####

library(stringr)
library(plyr)
library(reshape2)
library(stats)
library(ggplot2)
library(ggpmisc)
library(EnvStats)
library(rkt)
library(plyr)
library(tidyverse)
####note: we've changed subsampling number to 450 and removed samples with total abundance <450 for greater statistical consistency in richness between samples over time

bibi<-BIBI_pssb

#########  Format, remove data from wrong months, average scores from samples with same site+year, and subset B-IBI data to data since 2001 with more than 9 years of data
bibi$month <- format(as.Date(bibi$event_date, '%m/%d/%Y'), "%m") # add month column
bibi$year <- format(as.Date(bibi$event_date, '%m/%d/%Y'), "%Y") # add year column

bibi$sampnum <- paste(bibi$site_code, format(as.Date(bibi$event_date, '%m/%d/%Y'), "%y"), sep = "_") # add sample identifier

bibi <- bibi |> filter(month %in% c('07', '08', '09', '10'))



dups<-bibi[duplicated(bibi[,c("project" , "wria_number"  , "basin" , "subbasin" , "site_code", "year")]),]

names(bibi)

bibi_means <- aggregate(cbind(overall_score, taxa_richness_score, ephemeroptera_richness_score, plecoptera_richness_score, trichoptera_richness_score, clinger_richness_score, long_lived_richness_score, intolerant_richness_score, percent_dominant_score,  predator_percent_score, tolerant_percent_score) ~   site_code + subbasin + wria_number + year, data = bibi, FUN = "mean", na.action = na.exclude)

write.csv(bibi_means, paste0(output_folder,"/BIBI_means.csv"))

bibi.lr<-bibi_means
bibi.lr1<-bibi.lr

ggplot(bibi.lr, aes(x=site_code, y=year))+geom_tile()

ggplot(bibi.lr, aes(x=year, y=overall_score))+geom_smooth(se=F)+theme(legend.position = "none")+
  geom_point()



library(ggpmisc)

meanscores_working<-bibi.lr1 |> group_by(year) |> 
  dplyr::summarise(
    meanScore=mean(overall_score), 
    nyears=length(unique(overall_score)), 
    sd=sd(overall_score), 
    min=min(overall_score), 
    max=max(overall_score), 
    median=median(overall_score)) |> 
  ungroup()


#this has code that has been deprecated, needs fixing.
meanscore<-ggplot(meanscores_working, aes(x=year, y=meanScore)) + 
  geom_smooth(method="lm",se=F)+
  theme(legend.position = "none")+
  geom_point()+
  ylim(c(0,100))+labs(y="BIBI Score")+
  stat_poly_eq(formula = y ~ x, rr.digits = 2, coef.digits = 3, size = 4, color = "black", aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), parse = TRUE) 

median(bibi.lr1$overall_score)

bibi.lr1$year_num <- (as.numeric(bibi.lr1$year))
median(bibi.lr1$year_num)

##calculating intercept for RKT derived Sen-Theil line for Overall Score (0.38 is the Sen Theil slope for the RKT of overall score, calculated later in the code), from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf pdf page 15, and https://rdrr.io/github/USGS-R/HASP/src/R/statistics.R
median(bibi.lr1$overall_score)-(0.38*median(bibi.lr1$year_num))


##These lines are manually drawn and need updating! 11/5/2026
scores<-ggplot(bibi.lr1 )+
  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=overall_score, x=as.numeric(paste0(year_num)), group=as.factor(year_num)), fill="white")+
  geom_point(data=meanscores_working, mapping=aes(x=as.numeric(paste0(year)), y=meanScore))+
  # geom_smooth(data=meanscores, mapping=aes(x=as.numeric(paste0(Year)), y=meanScore), method="lm",se=F)+
  labs(y="Overall Score", x="Year")+ theme(text = element_text(size = 28))+
  # geom_quantile(data=bibi.lr1, mapping=aes(x=as.numeric(paste0(Year)), y=Overall.Score), quantiles=.5, color="red", size=1)#+
  geom_abline(slope=0.38, intercept = -711.2844, color="red", lwd=1)##adding Sen-Theil line and intercept calculated above


names(bibi.lr1)

l=unique(c(as.character(bibi.lr1$site_code)))
bibi.lr1 <- bibi.lr1[complete.cases(bibi.lr1[, 1]), ]
str(bibi.lr1)
bibi.lr1$site_num<-as.numeric(factor(bibi.lr1$site_code, levels=l))

metrics<-colnames(bibi.lr1)[5:(ncol(bibi.lr1)-1)]

library(rkt)

for (each in metrics) {
  bibi.lr1$rkt.pval <- with(bibi.lr1, rkt(date = year_num, y = get(each), block = site_num, correct = FALSE)[1])
  bibi.lr1$rkt.Sen <- with(bibi.lr1, rkt(date = year_num, y = get(each), block = site_num, correct = FALSE)[3])
  bibi.lr1$rkt.tau <- with(bibi.lr1, rkt(date = year_num, y = get(each), block = site_num, correct = FALSE)[12])
  bibi.lr1$RKTtrend <- with(bibi.lr1, ifelse(rkt.pval < 0.05 & rkt.tau < 0, "negative", ifelse(rkt.pval < 0.05 & rkt.tau > 0, "positive", "none")))
  colnames(bibi.lr1)[ncol(bibi.lr1)-3]<-paste0("rkt.pval_", each)
  colnames(bibi.lr1)[ncol(bibi.lr1)-2]<-paste0("rkt.Sen_", each)
  colnames(bibi.lr1)[ncol(bibi.lr1)-1]<-paste0("rkt.tau_", each)
  colnames(bibi.lr1)[ncol(bibi.lr1)]<-paste0("RKTtrend_", each)
}

overall_trends<-unique(bibi.lr1[(ncol(bibi.lr)+2):(ncol(bibi.lr1))])
overall_trends = as.matrix(overall_trends)
write.csv(overall_trends, paste0(output_folder,"/RKT_overall_trends_PSSB.csv"))
trend<-subset(overall_trends, select=c(str_subset(colnames(overall_trends), "RKTtrend_")))

##plot the trends for each metric
for (i in 5:15) {
  ylabel<-names(bibi.lr1[i])
  # p<-ggplot(bibi.lr1, aes(as.factor(year), bibi.lr1[,i]))+geom_violin(draw_quantiles=c(.25, .5, .75))+ stat_summary(fun.y=median, geom="point", size=2, color="red")+labs(y=ylabel)
  o<-ggplot(bibi.lr1, aes(as.factor(year_num), bibi.lr1[,i]))+geom_boxplot()+labs(y=ylabel)#+geom_abline(aes(intercept=get(paste0("MKinter_",names(bibi.lr1[i]))),slope=get(paste0("MKSlope_",names(bibi.lr1[i])))))
  # print(p)
  print(o)
  ggsave(paste0(output_folder,"/metrics_boxplots_", ylabel, "PSSB.png"), plot = o, width=20, height=10)
}


# ########## basin specific RKT regressions ####
# 
# bibi.lr2 <- plyr::ddply(bibi.lr, 'Subbasin', mutate, timerange = (as.numeric(max(Year))-as.numeric(min(Year))))
# bibi.lr2 <- plyr::ddply(bibi.lr2, 'Subbasin', mutate, por = length(unique(Year))) # add on period of record col
# bibi.lr2<-subset(bibi.lr2, Subbasin!="")
# 
# bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$por > 9,])
# bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$timerange > 9,])# look only at streams with >9 years of data
# plyr::ddply(bibi.lr2, 'Subbasin', summarize, count=length(unique(Year)))
# 
# unique(bibi.lr2$Site.Code)
# 
# bibi.lr2<-plyr::ddply(bibi.lr2, 'Subbasin', mutate, max = max(Year))
# 
# scal<-"Subbasin" # "Site.Code" or "Stream.or.River"
# 
# form<-names(bibi.lr2)[5:(ncol(bibi.lr2)-3)]
# bibi.lr2$Year<-as.numeric(paste0(bibi.lr2$Year))
# 
# 
# l=unique(c(as.character(bibi.lr2$Site.Code)))
# bibi.lr2$Site.Number<-as.numeric(factor(bibi.lr2$Site.Code, levels=l))
# 
# library(plyr)
# for (var in form) {
#   bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.tau=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$tau)
#   bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.Sen=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$B)
#   bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.pval=rkt( date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$sl)
#   bibi.lr2$MKtrend <- ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau < 0, "negative", ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau > 0, "positive", "none"))
#   colnames(bibi.lr2)[ncol(bibi.lr2)-3]<-paste0("RKTtau_", var)
#   colnames(bibi.lr2)[ncol(bibi.lr2)-2]<-paste0("RKTSlope_", var)
#   colnames(bibi.lr2)[ncol(bibi.lr2)-1]<-paste0("RKTpval_", var)
#   colnames(bibi.lr2)[ncol(bibi.lr2)]<-paste0("RKTtrend_", var)
# }
# 
# names(bibi.lr2)
# overall_trends<-unique(bibi.lr2[c(3,21:(ncol(bibi.lr2)))])
# overall_trends<-unique(subset(overall_trends, select=c("Subbasin" ,str_subset(colnames(overall_trends), "RKTtau"), str_subset(colnames(overall_trends), "RKTpval"),str_subset(colnames(overall_trends), "RKTtrend"))))
# 
# 
# 
# write.csv(overall_trends, paste0(output_folder,"/Subbasin_RKTtrends_PSSB.csv"))
# 
# 
# bibi_sum<-ddply(bibi.lr2, "Subbasin", summarise, mean_Overall.Score=mean(Overall.Score), SD_overallscore=sd(Overall.Score), SE_overallscore=sd(Overall.Score)/sqrt(length(Overall.Score)), len=length(Overall.Score), RKT_trend=unique(RKTtrend_Overall.Score), Avg_RKTSlope_Overall.Score=mean(RKTSlope_Overall.Score), Avg_RKTtau_Overall.Score=mean(RKTtau_Overall.Score), Avg_RKTpval_Overall.Score=mean(RKTpval_Overall.Score), sites=length(unique(Site.Code)), min=min(Overall.Score), max=max(Overall.Score), med=median(Overall.Score))
# write.csv(bibi_sum, paste0(output_folder,"/score_RKTtrend_summarized_by_subbasin_PSSB.csv"))
# 
# pp<-ggplot(bibi.lr2, aes(Year, Overall.Score, group=Year))+geom_boxplot(aes(fill=RKTtrend_Overall.Score))+facet_wrap(~Subbasin, ncol=4)+labs(y="B-IBI Score", x="Year")+ theme(axis.text.x = element_text(angle = 90))
# ggsave(paste0(output_folder,"/SubbasinRKTtrend_PSSB.png"), plot = pp, width=20, height=10)
# 
# 
# 
# ########################
# ########## watershed specific RKT regressions ####
# 
# bibi.lr2 <- ddply(bibi.lr, 'WRIA.Number', mutate, timerange = (as.numeric(max(Year))-as.numeric(min(Year))))
# bibi.lr2 <- ddply(bibi.lr2, 'WRIA.Number', mutate, por = length(unique(Year))) # add on period of record col
# bibi.lr2<-subset(bibi.lr2, WRIA.Number!="")
# 
# bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$por > 9,])
# bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$timerange > 9,])# look only at streams with >9 years of data
# ddply(bibi.lr2, 'WRIA.Number', summarize, count=length(unique(Year)))
# 
# unique(bibi.lr2$Site.Code)
# 
# bibi.lr2<-ddply(bibi.lr2, 'WRIA.Number', mutate, max = max(Year))
# 
# scal<-"WRIA.Number" # "Site.Code" or "Stream.or.River"
# 
# form<-names(bibi.lr2)[5:(ncol(bibi.lr2)-3)]
# bibi.lr2$Year<-as.numeric(paste0(bibi.lr2$Year))
# 
# 
# l=unique(c(as.character(bibi.lr2$Site.Code)))
# bibi.lr2$Site.Number<-as.numeric(factor(bibi.lr2$Site.Code, levels=l))
# 
# for (var in form) {
#   bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.tau=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$tau)
#   bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.Sen=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$B)
#   bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.pval=rkt( date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$sl)
#   bibi.lr2$MKtrend <- ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau < 0, "negative", ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau > 0, "positive", "none"))
#   colnames(bibi.lr2)[ncol(bibi.lr2)-3]<-paste0("RKTtau_", var)
#   colnames(bibi.lr2)[ncol(bibi.lr2)-2]<-paste0("RKTSlope_", var)
#   colnames(bibi.lr2)[ncol(bibi.lr2)-1]<-paste0("RKTpval_", var)
#   colnames(bibi.lr2)[ncol(bibi.lr2)]<-paste0("RKTtrend_", var)
# }
# 
# names(bibi.lr2)
# overall_trends<-unique(bibi.lr2[c(4,21:(ncol(bibi.lr2)))])
# overall_trends<-unique(subset(overall_trends, select=c("WRIA.Number" ,str_subset(colnames(overall_trends), "RKTtau"), str_subset(colnames(overall_trends), "RKTpval"),str_subset(colnames(overall_trends), "RKTtrend"))))
# 
# 
# 
# write.csv(overall_trends, paste0(output_folder,"/Trend_watershed_RKTtrends_PSSB.csv"))
# 
# 
# bibi_sum<-ddply(bibi.lr2, "WRIA.Number", summarise, mean_Overall.Score=mean(Overall.Score), SD_overallscore=sd(Overall.Score), SE_overallscore=sd(Overall.Score)/sqrt(length(Overall.Score)), len=length(Overall.Score), RKT_trend=unique(RKTtrend_Overall.Score), Avg_RKTSlope_Overall.Score=mean(RKTSlope_Overall.Score), Avg_RKTtau_Overall.Score=mean(RKTtau_Overall.Score), Avg_RKTpval_Overall.Score=mean(RKTpval_Overall.Score), sites=length(unique(Site.Code)), min=min(Overall.Score), max=max(Overall.Score), med=median(Overall.Score))
# write.csv(bibi_sum, paste0(output_folder,"/score_RKTtrend_summarized_by_watershed_PSSB.csv"))
# 
# pp<-ggplot(bibi.lr2, aes(Year, Overall.Score, group=Year))+geom_boxplot(aes(fill=RKTtrend_Overall.Score))+facet_wrap(~WRIA.Number, ncol=4)+labs(y="B-IBI Score", x="Year")+ theme(axis.text.x = element_text(angle = 90))
# ggsave(paste0(output_folder,"/WatershedRKTtrend_PSSB.png"), plot = pp, width=20, height=10)
# 



#####Site trends ####
library(EnvStats)

bibi<-bibi_means

bibi<-bibi %>%
  dplyr::group_by(site_code) %>%
  dplyr::mutate(por = n_distinct(year)) %>% # add on period of record col
  ungroup()

unique(bibi$site_code)
bibi <- bibi |> filter(por > 9) # look only at sites with >9 years of data

bibi$year_num <- as.numeric(bibi$year)
bibi <- 
  bibi %>%
  dplyr::group_by(site_code) %>%
  dplyr::mutate(timerange = max(year_num) - min(year_num)) %>% # add on period of record col
  ungroup()

bibi<-bibi |> filter(timerange > 9)# look only at streams with >9 years of data

unique(bibi$site_code)

ggplot(bibi, aes(x=site_code, y=year))+geom_tile()

ggplot(bibi, aes(x=year, y=overall_score))+geom_smooth(se=F)+theme(legend.position = "none")+
  geom_point()



library(ggpmisc)

meanscores_sites<- bibi |> select(site_code, overall_score, year) |> 
  dplyr::group_by(year) |> 
  dplyr::summarise(meanScore = mean(overall_score),
                   nyears = n_distinct(site_code),
                   sd = sd(overall_score),
                   min = min(overall_score),
                   max = max(overall_score),
                   median = median(overall_score))

meanscore<-ggplot(meanscores_sites, aes(x=year, y=meanScore))+geom_smooth(method="lm",se=F)+theme(legend.position = "none")+
  geom_point()+ylim(c(0,100))+labs(y="BIBI Score")+
  stat_poly_eq(formula = y ~ x, rr.digits = 2, coef.digits = 3, size = 4, color = "black", aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), parse = TRUE) 

scores<-ggplot(bibi )+
  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=overall_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  geom_point(data=meanscores_working, mapping=aes(x=as.numeric(paste0(year)), y=meanScore))+
  geom_smooth(data=meanscores_sites, mapping=aes(x=as.numeric(paste0(year)), y=meanScore), method="lm",se=F)+
  labs(y="Overall Score", x="Year")+ theme(text = element_text(size = 28))+
  geom_quantile(data=bibi, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.5, color="red", linewidth = 1)#+
# geom_abline(slope=0.48, intercept = -914.2032, color="red", lwd=1)##adding Sen-Theil line and intercept calculated above



scal<-"site_code" # "Site.Code" or "Subbasin"

bibi<-subset(bibi, site_code!="")
form<-names(bibi)[5:(ncol(bibi)-3)]
# bibi.lr[bibi.lr$Year<2012,"Overall.Score"]<-(bibi.lr[bibi.lr$Year<2012,"Overall.Score"]*.972)+5.01
# for (i in 9:(ncol(bibi.lr))) {
for (var in form) {
  bibi<- ddply(bibi, .variables =scal, mutate, slope =lm(get(paste0(var)) ~ as.numeric(year))[[1]][2])
  bibi<-ddply(bibi, .variables =scal, mutate, pval = summary(lm(get(paste0(var)) ~ as.numeric(year)))[[4]][2,4])
  bibi$trend <- ifelse(bibi$pval < 0.05 & bibi$slope < 0, "negative", ifelse(bibi$pval < 0.05 & bibi$slope > 0, "positive", "none")) # add trend col
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.tau=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[1])
  # bibi<-ddply(bibi,.(get(scal)),mutate, mk.tau=kendallTrendTest(as.formula(paste0(colnames(bibi[i]), " ~ ", "as.numeric(Year)")))$estimate[1])
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.Sen=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[2])
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.pval=kendallTrendTest( get(paste0(var))~ as.numeric(year))$p.value)
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.inter=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[3])
  bibi$MKtrend <- ifelse(bibi$mk.pval < 0.05 & bibi$mk.tau < 0, "negative", ifelse(bibi$mk.pval < 0.05 & bibi$mk.tau > 0, "positive", "none"))
  colnames(bibi)[ncol(bibi)-7]<-paste0("slope_",var)
  colnames(bibi)[ncol(bibi)-6]<-paste0("LRpval_", var)
  colnames(bibi)[ncol(bibi)-5]<-paste0("LRtrend_", var)
  colnames(bibi)[ncol(bibi)-4]<-paste0("MKtau_", var)
  colnames(bibi)[ncol(bibi)-3]<-paste0("MKSlope_", var)
  colnames(bibi)[ncol(bibi)-2]<-paste0("MKpval_", var)
  colnames(bibi)[ncol(bibi)-1]<-paste0("MKinter_", var)
  colnames(bibi)[ncol(bibi)]<-paste0("MKtrend_", var)
}

sampletimes<- ggplot(bibi, aes(x=site_code, y=year, fill= MKSlope_overall_score), fill="black")+
  geom_tile(color="black")+
  scale_fill_gradient2(name="Sen-Theil Slope")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# ggsave("sampletimes.png", plot=sampletimes, width=20, height=10)

names(bibi)
overall_trends<-unique(bibi[c(2:3,19:(ncol(bibi)))])

# overall_trends<-unique(subset(overall_trends, select=c("Site.Code" , "Subbasin", "Stream.or.River", str_subset(colnames(overall_trends), "MKtau"), str_subset(colnames(overall_trends), "MKpval"),str_subset(colnames(overall_trends), "MKtrend"))))
write.csv(overall_trends, paste0(output_folder,"/Site_Trends_coarse_PSSB.csv"))

bibi_sum_site<-ddply(bibi, "Site.Code", summarise, mean_Overall.Score=mean(overall_score), SD_overallscore=sd(Overall.Score), SE_overallscore=sd(Overall.Score)/sqrt(length(Overall.Score)), len=length(Overall.Score), MK_trend=unique(MKtrend_Overall.Score), Avg_MKSlope_Overall.Score=mean(MKSlope_Overall.Score), Avg_MKtau_Overall.Score=mean(MKtau_Overall.Score), Avg_MKpval_Overall.Score=mean(MKpval_Overall.Score), sites=length(unique(Site.Code)), min=min(Overall.Score), max=max(Overall.Score), med=median(Overall.Score))
write.csv(bibi_sum_site,paste0(output_folder,"/score_RKTtrend_summarized_by_site_PSSB.csv"))

sum(bibi_sum_site$len)

############  PLOTTING

############  Plot all positive and negative trends for overall and individual B-IBI scores
# for (var in colnames(bibi[str_subset(colnames(bibi), "MKtrend")])) {
#   if(nrow(bibi[bibi[var]=='positive',])>=1) {
#     pp<-ggplot(bibi[bibi[var]=='positive',], aes(as.numeric(Year), get(str_subset(colnames(bibi[,c(1:18)]), str_split(var, "_", n=2)[[1]][2])), group = get(scal), fill = get(scal), colour = get(scal))) + ylim(c(0,100)) + geom_point(colour = 'black', size = 5, pch = 21) + theme(legend.position="none") + 
#       facet_wrap(~ WRIA.Number +Subbasin+get(scal), ncol = 4) +  #factor(stream2) + 
#       theme(axis.title=element_text(size=20,face="bold"), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90), strip.text = element_text(size=12)) + labs(y = var, x = "Year")
#     pp<-pp+geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])),  colour=get(scal)))
#     pp<-pp+geom_text(mapping= aes(x=2012, y=6, label=paste("Slope = ", round(get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])), 3))), colour="black")+theme(text=element_text(size=20), plot.title=element_text(hjust=.5, size=24, face="bold"))
#     # print(pp)
#     ggsave(paste0("BIBI~timeUp_", scal, "_", var,   "Bellevue.png"), plot = pp, width=20, height=10)
#   }
#   if(nrow(bibi[bibi[var]=='negative',])>=1) {
#     pn<-ggplot(bibi[bibi[var]=='negative',], aes(as.numeric(Year), get(str_subset(colnames(bibi[,c(1:18)]), str_split(var, "_", n=2)[[1]][2])), group = get(scal), fill = get(scal), colour = get(scal))) + 
#       ylim(c(0,100)) + geom_point(colour = 'black', size = 5, pch = 21) + theme(legend.position="none") + 
#       #stat_smooth(method = "lm", formula = y ~ x) + 
#       facet_wrap(~ WRIA.Number +Subbasin + get(scal), ncol = 4) + #factor(stream2) +
#       #stat_poly_eq(formula = y ~ x, label.y = 92, rr.digits = 2, coef.digits = 3, size = 4, color = "black", aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~")), parse = TRUE) + 
#       theme(axis.title=element_text(size=20,face="bold"), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90), strip.text = element_text(size=12)) + 
#       labs(y = var, x = "Year")
#     pn<-pn+geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])),  colour=get(scal)))
#     pn<-pn+geom_text(mapping= aes(x=2012, y=6, label=paste("Slope = ", round(get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])), 3))), colour="black")+theme(text=element_text(size=20), plot.title=element_text(hjust=.5, size=24, face="bold"))
#     # print(pp)
#     ggsave(paste0("BIBI~timeDown_", scal, "_",var, "Bellevue.png"), plot = pn, width=20, height=10)
#   }
# }
library(RColorBrewer)

# for (var in colnames(bibi[str_subset(colnames(bibi), "trend")])) {
# if(nrow(bibi[bibi[var]=='positive',])>=1) {
var="MKtrend_Overall.Score"

bibi$pchOS<-as.factor(bibi$MKtrend_Overall.Score)
shapes = c(6, 4, 2) 
shapes <- shapes[as.numeric(bibi$pchOS)]


pp<-ggplot(bibi, aes(as.numeric(Year), get(str_subset(colnames(bibi[,c(1:18)]), str_split(var, "_", n=2)[[1]][2])), group = get(scal), fill = get(scal), colour = get(scal))) + 
  ylim(c(0,100)) + 
  geom_point( aes(shape = c(shapes))) + scale_shape_identity()+ #theme(legend.position="none") +
  facet_wrap(~ Subbasin, ncol = 4) +  #factor(stream2) + 
  theme(legend.position="none",axis.title=element_text(size=20,face="bold"), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90), strip.text = element_text(size=12)) + labs(y = var, x = "Year")
pp<-pp+geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])),  colour=get(scal)))+
  scale_colour_manual(values=rep(brewer.pal(8,"Dark2"), times=30))
# pp<-pp+geom_text(mapping= aes(x=2012, y=6, label=paste("Slope = ", round(get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])), 3))), colour="black")+theme(text=element_text(size=20), plot.title=element_text(hjust=.5, size=24, face="bold"))
# ggsave(paste0("WRIA9_BIBI~timeUp_", scal, "_", var,   "-2.png"), plot = pp, width=20, height=10)
# }
print(pp)
# }
ggsave(paste0(output_folder,"/OverallScores_Sites_PSSB.png"), plot = pp, width=20, height=10)


ggplot(bibi, aes(x=Year, y=Overall.Score))+  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  # geom_smooth(se=F, method="lm")+
  geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2]))), color="blue", size=1)+
  theme(legend.position = "none")+
  facet_wrap(~Site.Code)+geom_point( aes(shape = c(shapes))) + scale_shape_identity()+labs(y="Overall B-IBI Score")


ggsave(paste0(output_folder,"/SiteScores_PSSB.png"), width=20, height=10)
