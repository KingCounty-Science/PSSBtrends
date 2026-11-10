#=== === === === === === === ===
# Script started by Rebekah Stiling October 2026
# This script uses B-IBI scores from PSSB to calculate trends over time. Much of the script is recycled from the 04_BIBI-trends_for-SFS-2026.R script I made in my, based on scripts Beth Sosik left behind.
# rstiling@kingcounty.gov
#=== === === === === === === ===

#load relevant libraries
library(tidyverse) #for wrangling
library(janitor) #for cleaning names.

# Specifications:
# only  2002 and after
# Projects: Ambient, Vashon, Miller-Walker, Boise 
# BIBI scores
# SHOULD I REDOWNLOAD AND ONLY USE CORSE??? ####

#Read in the data #This is all scores from PSSB from all of time from Boise, Vashon, Miller Walker, and Ambient

scores<-list.files("Inputs/Scores_Oct2026", full.names = TRUE) |> #List all files names in the folder
  discard(~ grepl("metadata", .x, ignore.case = TRUE)) |> #drop the metadata files which are just info
    map(read_delim)  |> #this "read_delim()" to each item in the list of file names.
    bind_rows() |>  #this puts all the files together
  clean_names() #clean column names

head(scores)

#create empty folder for graphs and csvs
output_folder <- paste0("output",Sys.Date()) 
dir.create(output_folder)

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



bibi_cleaning <-BIBI_pssb

## cleaning ####
bibi_cleaning <-BIBI_pssb

# remove data from wrong months, average scores from samples with same site+year, and subset B-IBI data to data since 2001 with more than 9 years of data 
bibi_cleaning$month <- format(as.Date(bibi_cleaning$event_date, '%m/%d/%Y'), "%m") # add month column
bibi_cleaning$year <- format(as.Date(bibi_cleaning$event_date, '%m/%d/%Y'), "%Y") # add year column

bibi_cleaning$sampnum <- paste(bibi_cleaning$site_code, format(as.Date(bibi_cleaning$event_date, '%m/%d/%Y'), "%y"), sep = "_") # add sample identifier

bibi_cleaning <- bibi_cleaning |> filter(month %in% c('07', '08', '09', '10'))



dups<-bibi_cleaning[duplicated(bibi_cleaning[,c("project" , "wria_number"  , "basin" , "subbasin" , "site_code", "year")]),]

names(bibi_cleaning)

#bibi_means has site year metrics
bibi_means <- aggregate(cbind(overall_score, taxa_richness_score, ephemeroptera_richness_score, plecoptera_richness_score, trichoptera_richness_score, clinger_richness_score, long_lived_richness_score, intolerant_richness_score, percent_dominant_score,  predator_percent_score, tolerant_percent_score) ~   site_code + subbasin + wria_number + year, data = bibi_cleaning, FUN = "mean", na.action = na.exclude)

#Get year informaion associated with bibi_means

bibi_years<-bibi_means |> #the new table will be 
  dplyr::group_by(site_code) |> # group the data by site code
  dplyr::summarize(n_years = n(), #for each site code, count the number of years in the group, then
                   min_years = as.numeric(min(year)), #record the minimum year (turn to numeric), then
                   max_years = as.numeric(max(year)), #record the maximum year (turn to numeric), then
                   duration = max_years - min_years) #determine the duration of years

bibi_means_years_10 <- left_join(bibi_means, bibi_years, by = "site_code") |> #keep all rows of bibi_means, and add the info from bibi_years)
  filter(n_years > 9)


write.csv(bibi_means, paste0(output_folder,"/BIBI_means.csv"))
write.csv(bibi_means_years_10, paste0(output_folder,"/bibi_means_years.csv"))


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
  geom_point(data=meanscores_sites, mapping=aes(x=as.numeric(paste0(year)), y=meanScore))+
  geom_smooth(data=meanscores_sites, mapping=aes(x=as.numeric(paste0(year)), y=meanScore), method="lm",se=F)+
  labs(y="Overall Score", x="Year")+ theme(text = element_text(size = 28))+
  geom_quantile(data=bibi, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.5, color="red", linewidth = 1)#+
# geom_abline(slope=0.48, intercept = -914.2032, color="red", lwd=1)##adding Sen-Theil line and intercept calculated above



scal<-"site_code" # "Site.Code" or "Subbasin"

bibi<-subset(bibi, site_code!="")
form<-names(bibi)[5:(ncol(bibi)-3)]
# bibi.lr[bibi.lr$Year<2012,"Overall.Score"]<-(bibi.lr[bibi.lr$Year<2012,"Overall.Score"]*.972)+5.01
# for (i in 9:(ncol(bibi.lr))) {

#this loops through each of the 11 BiBi metrics and calculates the trend according to site, then records the output repeatedly for that site along all years. Why do it this way?

library(plyr) ## Figure out how to re-write this without the plyr package

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

#This is a creative way of showing both the duration of the data and the overall slope.
sampletimes<- ggplot(bibi, 
                     aes(x=site_code, y=year, fill= MKSlope_overall_score))+
  geom_tile()+
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

##This step gets one trend per site.
overall_trends <- bibi |> 
  select(site_code, subbasin, timerange:(ncol(bibi))) |> 
  unique()

# overall_trends<-unique(subset(overall_trends, select=c("Site.Code" , "Subbasin", "Stream.or.River", str_subset(colnames(overall_trends), "MKtau"), str_subset(colnames(overall_trends), "MKpval"),str_subset(colnames(overall_trends), "MKtrend"))))
write.csv(overall_trends, paste0(output_folder,"/Site_Trends_coarse_PSSB.csv"))

bibi_sum_site<-ddply(bibi, "site_code", summarise, mean_Overall.Score=mean(overall_score), SD_overallscore=sd(overall_score), SE_overallscore=sd(overall_score)/sqrt(length(overall_score)), len=length(overall_score), MK_trend=unique(MKtrend_overall_score), Avg_MKSlope_Overall.Score=mean(MKSlope_overall_score), Avg_MKtau_Overall.Score=mean(MKtau_overall_score), Avg_MKpval_Overall.Score=mean(MKpval_overall_score), sites=length(unique(site_code)), min=min(overall_score), max=max(overall_score), med=median(overall_score))
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
var="MKtrend_overall_score"

bibi$pchOS<-as.factor(bibi$MKtrend_overall_score)
shapes = c(6, 4, 2) 
shapes <- shapes[as.numeric(bibi$pchOS)]


pp<-ggplot(bibi, aes(as.numeric(year), get(str_subset(colnames(bibi[,c(1:18)]), str_split(var, "_", n=2)[[1]][2])), group = get(scal), fill = get(scal), colour = get(scal))) + 
  ylim(c(0,100)) + 
  geom_point( aes(shape = c(shapes))) + scale_shape_identity()+ #theme(legend.position="none") +
  facet_wrap(~ subbasin, ncol = 4) +  #factor(stream2) + 
  theme(legend.position="none",axis.title=element_text(size=20,face="bold"), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90), strip.text = element_text(size=12)) + labs(y = var, x = "Year")

pp<-pp+geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])),  colour=get(scal)))+
  scale_colour_manual(values=rep(brewer.pal(8,"Dark2"), times=30))
# pp<-pp+geom_text(mapping= aes(x=2012, y=6, label=paste("Slope = ", round(get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])), 3))), colour="black")+theme(text=element_text(size=20), plot.title=element_text(hjust=.5, size=24, face="bold"))
# ggsave(paste0("WRIA9_BIBI~timeUp_", scal, "_", var,   "-2.png"), plot = pp, width=20, height=10)
# }
print(pp)
# }
ggsave(paste0(output_folder,"/OverallScores_Sites_PSSB.png"), plot = pp, width=20, height=10)


ggplot(bibi, aes(x=year, y=overall_score))+  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  # geom_smooth(se=F, method="lm")+
  geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2]))), color="blue", size=1)+
  theme(legend.position = "none")+
  facet_wrap(~site_code)+geom_point( aes(shape = c(shapes))) + scale_shape_identity()+labs(y="Overall B-IBI Score")


ggsave(paste0(output_folder,"/SiteScores_PSSB.png"), width=20, height=10)
