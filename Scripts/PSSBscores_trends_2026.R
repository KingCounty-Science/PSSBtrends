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
# COARSE BIBI scores

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
BIBI_pssb |> select(number_of_organisms) |> unique()

## cleaning ####
bibi_cleaning <-BIBI_pssb

# remove data from wrong months, average scores from samples with same site+year, and subset B-IBI data to data since 2001 with more than 9 years of data 
bibi_cleaning$month <- format(as.Date(bibi_cleaning$event_date, '%m/%d/%Y'), "%m") # add month column
bibi_cleaning$year <- format(as.Date(bibi_cleaning$event_date, '%m/%d/%Y'), "%Y") # add year column

bibi_cleaning$sampnum <- paste(bibi_cleaning$site_code, format(as.Date(bibi_cleaning$event_date, '%m/%d/%Y'), "%y"), sep = "_") # add sample identifier

bibi_cleaning <- bibi_cleaning |> 
  filter(month %in% c('07', '08', '09', '10')) |> #only keep samples collected in sample window
  filter_out(number_of_organisms <400) # removes samples with fewer than 400 organisms - review this. It was 450 orgs, but PSSB is doing some weird things in giving 449 orgs in sample when we ask for 450 - see 09DUW0091 2002 sample

#This only gives us the 2nd and above observation, we don't subset the initial instance of a sitecode & year.
dups<-bibi_cleaning[duplicated(bibi_cleaning[,c("project" , "wria_number"  , "basin" , "subbasin" , "site_code", "year")]),]

names(bibi_cleaning)

#bibi_means has site year metrics; this takes the average of metrics (including scores) if there is more than one sample per year per site
bibi_means <- aggregate(cbind(overall_score, taxa_richness_score, ephemeroptera_richness_score, plecoptera_richness_score, trichoptera_richness_score, clinger_richness_score, long_lived_richness_score, intolerant_richness_score, percent_dominant_score,  predator_percent_score, tolerant_percent_score, ephemeroptera_richness_quantity, plecoptera_richness_quantity, trichoptera_richness_quantity, ept_richness_quantity) ~   site_code + subbasin + wria_number + year + number_of_organisms, data = bibi_cleaning, FUN = "mean", na.action = na.exclude)

#Get year information associated with bibi_means

bibi_years<-bibi_means |> #the new table will be 
  dplyr::group_by(site_code) |> # group the data by site code
  dplyr::summarize(n_years = n(), #for each site code, count the number of years in the group, then
                   min_years = as.numeric(min(year)), #record the minimum year (turn to numeric), then
                   max_years = as.numeric(max(year)), #record the maximum year (turn to numeric), then
                   duration = max_years - min_years) #determine the duration of years

#keep all rows of bibi_means, and add the info from bibi_years, and filters out sites with less than 10 years of data
bibi_means_years_10 <- left_join(bibi_means, bibi_years, by = "site_code") |> 
  filter(n_years > 9)

write.csv(bibi_means, paste0(output_folder,"/BIBI_means.csv")) #this is whole data set
write.csv(bibi_means_years_10, paste0(output_folder,"/bibi_means_w_at_least_10yrs.csv")) #this is data set with sites that have at least 10 years of data

#plot of "min_years" to see the starting year of record for each site; this allows us to see that some sites are started late in record (WRIA 7 and BOise largely)
ggplot(bibi_years, aes(x = min_years))+
  geom_histogram(binwidth = 1)

#create a vector with the sites that had collection begin in 2014 or after
later_starts<-bibi_years |> filter(min_years >= 2014) |> select(site_code) |> unique() |> pull()


#####Site trends ####
library(EnvStats)

####Use the data set "bibi_means" ONLY for site level trends and ONLY if you want all sites; don't use this dataset for overall trends as this still includes sites in WRIA 7 and Boise that were started later
bibi <- bibi_means |> 
  filter_out(year == 2026) #drop incomplete 2026 data

bibi<-bibi %>%
  dplyr::group_by(site_code) %>%
  dplyr::mutate(por = n_distinct(year)) %>% # add on period of record col
  ungroup()

unique(bibi$site_code)
bibi <- bibi |> filter(por > 9) # look only at sites with >9 years of data

bibi$year_num <- as.numeric(bibi$year) #add a year column as numeric
bibi <- 
  bibi %>%
  dplyr::group_by(site_code) %>%
  dplyr::mutate(timerange = max(year_num) - min(year_num)) %>% # add on range of record col
  ungroup()

#####Use this code to generate the "no late starts" data set for looking at 125 sites that were started at are near 2002
bibi_no_late_starts <- bibi_means |> 
  filter_out((year==2026)) |> 
  filter_out(site_code %in% later_starts) |> 
  dplyr::group_by(site_code) %>%
  dplyr::mutate(por = n_distinct(year)) %>% # add on period of record col
  ungroup() |> 
  filter(por > 9)
bibi_no_late_starts<-bibi_no_late_starts %>%  ### are these next 4 lines redundant with the previous few lines? are these needed?
  dplyr::group_by(site_code) %>%
  dplyr::mutate(por = n_distinct(year)) %>% # add on period of record col
  ungroup()

unique(bibi_no_late_starts$site_code)

bibi_no_late_starts$year_num <- as.numeric(bibi_no_late_starts$year) #add a year column as numeric
bibi_no_late_starts <- 
  bibi_no_late_starts %>%
  dplyr::group_by(site_code) %>%
  dplyr::mutate(timerange = max(year_num) - min(year_num)) %>% # add on range of record col
  ungroup()

unique(bibi_no_late_starts$site_code)

ggplot(bibi_no_late_starts, aes(x=site_code, y=year))+geom_tile()

ggplot(bibi_no_late_starts, aes(x=year, y=overall_score))+geom_smooth(se=F)+theme(legend.position = "none")+
  geom_point()

write.csv(bibi_no_late_starts, paste0(output_folder,"/BIBI_means_no_late_starts.csv"))

library(ggpmisc)

meanscores_sites<- bibi |> select(site_code, overall_score, year) |> 
  dplyr::group_by(year) |> 
  dplyr::summarise(meanScore = mean(overall_score),
                   nyears = n_distinct(site_code),
                   sd = sd(overall_score),
                   min = min(overall_score),
                   max = max(overall_score),
                   median = median(overall_score))

meanscores_sites_fewer<- bibi_no_late_starts |> select(site_code, overall_score, year) |>
  dplyr::group_by(year) |>
  dplyr::summarise(meanScore = mean(overall_score),
                   nyears = n_distinct(site_code),
                   sd = sd(overall_score),
                   min = min(overall_score),
                   max = max(overall_score),
                   median = median(overall_score))

#plot for trends across years of all sites, based on the appropriate number of sites (e.g., those that have a long duration of data, started in 2002-2013, only including samples with at least 400 orgs[see note above])
#just the mean scores
ggplot(meanscores_sites_fewer, aes(x=year, y=meanScore))+geom_smooth(method="lm",se=F)+theme(legend.position = "none")+
  geom_point()+ylim(c(0,100))+labs(y="BIBI Score")+
  stat_poly_eq(formula = y ~ x, rr.digits = 2, coef.digits = 3, size = 4, color = "black", aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), parse = TRUE) 
#scores and trend with background and with 10th, 50th, and 90th quantiles
plot_all_sites_trends <-ggplot(bibi_no_late_starts)+
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=overall_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  geom_point(data=meanscores_sites_fewer, mapping=aes(x=as.numeric(paste0(year)), y=meanScore))+
  #geom_smooth(data=meanscores_sites_fewer, mapping=aes(x=as.numeric(paste0(year)), y=meanScore), method="lm",se=F)+
  labs(y="BIBI Score", x="Year")+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  theme_classic()+
  theme(text = element_text(size = 50))+
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.5, color="blue", linewidth = 1.5) +
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.9, color="darkgreen", linewidth = 1.5) +
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.1, color="white", linewidth = 1.5) #+
  #theme_minimal() #takes away grey and shrinks font a little
plot_all_sites_trends
ggsave(paste0(output_folder,"/all_sites_trends_w_10th_50th_90th_quantiles.png"), plot = plot_all_sites_trends, width = 20, height = 15, units="in")

#scores and trend with background and with 50th  quantiles
plot_all_sites_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  #geom_boxplot(aes(y=overall_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  geom_violin(aes(y=overall_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  #geom_point(data=meanscores_sites_fewer, mapping=aes(x=as.numeric(paste0(year)), y=meanScore))+
  #geom_smooth(data=meanscores_sites_fewer, mapping=aes(x=as.numeric(paste0(year)), y=meanScore), method="lm",se=F)+
  labs(y="BIBI Score", x="Year")+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  theme_classic()+
  theme(text = element_text(size = 50))+
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.5, color="blue", linewidth = 1.5)
plot_all_sites_trends_50
ggsave(paste0(output_folder,"/all_sites_trends_w_50th_quantile_line_violin.png"), plot = plot_all_sites_trends_50, width = 20, height = 15, units="in")

#scores and trend with background and with 10th, 50th and 90th quantiles; violin plot
plot_all_sites_trends_10_50_90 <-ggplot(bibi_no_late_starts)+
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_violin(aes(y=overall_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  labs(y="BIBI Score", x="Year")+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  theme_classic()+
  theme(text = element_text(size = 50))+
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.5, color="blue", linewidth = 1.5)+ 
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.9, color="darkgreen", linewidth = 1.5) +
  geom_quantile(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), quantiles=.1, color="white", linewidth = 1.5) #+
plot_all_sites_trends_10_50_90
ggsave(paste0(output_folder,"/all_sites_trends_w_10_50_90th_quantile_lines_violin.png"), plot = plot_all_sites_trends_10_50_90, width = 20, height = 15, units="in")

#stonefly score box plot
plot_plecop_score_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=2, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=2, ymax=4, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=4, ymax=6, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=6, ymax=8, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=8, ymax=10, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=plecoptera_richness_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=plecoptera_richness_score), method="lm",se=F)+
  labs(y="Plecoptera richness score", x="Year")+ 
  theme_classic()+
  theme(text = element_text(size = 50))
plot_plecop_score_trends_50
ggsave(paste0(output_folder,"/all_sites_plecop_score_trend_box.png"), plot = plot_plecop_score_trends_50, width = 20, height = 15, units="in")

# stonelfy score violin
plot_plecop_score_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_violin(aes(y=plecoptera_richness_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="blue", alpha = 0.4)+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=plecoptera_richness_score), method="lm",se=T)+
  labs(y="Plecoptera richness score", x="Year")+ 
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14))
  theme_classic()+
  theme(text = element_text(size = 50))
plot_plecop_score_trends_50
ggsave(paste0(output_folder,"/all_sites_plecop_score_trend_violin.png"), plot = plot_plecop_score_trends_50, width = 20, height = 15, units="in")

#stonefly richness violin plot
plot_plecop_rich_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_violin(aes(y=plecoptera_richness_quantity, x=as.numeric(paste0(year)), group=as.factor(year)), fill="blue", alpha = 0.4)+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=plecoptera_richness_quantity), method="lm",se=T)+
  labs(y="Plecoptera richness", x="Year")+ 
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14))
  theme_classic()+
  theme(text = element_text(size = 50))
plot_plecop_rich_trends_50
ggsave(paste0(output_folder,"/all_sites_plecop_richness_trend_violin.png"), plot = plot_plecop_rich_trends_50, width = 20, height = 15, units="in")

#mayfly score box plot
plot_mayfly_score_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=2, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=2, ymax=4, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=4, ymax=6, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=6, ymax=8, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=8, ymax=10, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=ephemeroptera_richness_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="white")+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=ephemeroptera_richness_score), method="lm",se=F)+
  labs(y="Ephemeroptera richness score", x="Year")+ 
  theme_classic()+
  theme(text = element_text(size = 50))
plot_mayfly_score_trends_50
ggsave(paste0(output_folder,"/all_sites_ephemeroptera_score_trend.png"), plot = plot_mayfly_score_trends_50, width = 20, height = 15, units="in")

# mayfly score violin plots
plot_ephem_score_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_violin(aes(y=ephemeroptera_richness_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="grey", alpha = 0.5)+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=ephemeroptera_richness_score), method="lm",se=T)+
  labs(y="Ephemeroptera richness score", x="Year")+ 
  theme_classic()+
  theme(text = element_text(size = 50))
plot_ephem_score_trends_50
ggsave(paste0(output_folder,"/all_sites_ephemeroptera_score_trend_violin.png"), plot = plot_ephem_score_trends_50, width = 20, height = 15, units="in")

# mayfly richness violin plots
plot_ephem_rich_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_violin(aes(y=ephemeroptera_richness_quantity, x=as.numeric(paste0(year)), group=as.factor(year)), fill="grey", alpha = 0.5)+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=ephemeroptera_richness_quantity), method="lm",se=T)+
  labs(y="Ephemeroptera richness", x="Year")+ 
  theme_classic()+
  theme(text = element_text(size = 50))
plot_ephem_rich_trends_50
ggsave(paste0(output_folder,"/all_sites_ephemeroptera_richness_trend_violin.png"), plot = plot_ephem_rich_trends_50, width = 20, height = 15, units="in")

#EPT richness violin plot
plot_ept_rich_trends_50 <-ggplot(bibi_no_late_starts)+
  geom_violin(aes(y=ept_richness_quantity, x=as.numeric(paste0(year)), group=as.factor(year)), fill="pink", alpha = 0.8)+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=ept_richness_quantity), method="lm",se=T)+
  labs(y="EPT richness", x="Year")+ 
  theme_classic()+
  theme(text = element_text(size = 50))
plot_ept_rich_trends_50
ggsave(paste0(output_folder,"/all_sites_ept_richness_trend.png"), plot = plot_ept_rich_trends_50, width = 20, height = 15, units="in")

# overall scores, violin plots without score category colors
plot_overall_50 <-ggplot(bibi_no_late_starts)+
  geom_violin(aes(y=overall_score, x=as.numeric(paste0(year)), group=as.factor(year)), fill="black", alpha = 0.8)+
  geom_smooth(data=bibi_no_late_starts, mapping=aes(x=as.numeric(paste0(year)), y=overall_score), method="lm",se=T)+
  labs(y="BIBI Score", x="Year")+ 
  theme_classic()+
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  theme(text = element_text(size = 50))
plot_overall_50
ggsave(paste0(output_folder,"/all_sites_overall_trend_violin.png"), plot = plot_overall_50, width = 20, height = 15, units="in")


####
scal<-"site_code" # "Site.Code" or "Subbasin"

### code for calculating trends using all data in "bibi"; see below for code calculating trends using only data in "bibi_no_late_starts", which is the most filtered (most conservative) set of data
bibi<-subset(bibi, site_code!="")
form<-names(bibi)[5:(ncol(bibi)-3)]
# bibi.lr[bibi.lr$Year<2012,"Overall.Score"]<-(bibi.lr[bibi.lr$Year<2012,"Overall.Score"]*.972)+5.01
# for (i in 9:(ncol(bibi.lr))) {

#this loops through each of the 11 BiBi metrics and calculates the trend according to site, then records the output repeatedly for that site along all years. Why do it this way?

library(plyr) ## Figure out how to re-write this without the plyr package

for (var in form) {
  #First linear regression (slope, p-value, and interpretation)
  bibi<- ddply(bibi, .variables =scal, mutate, slope =lm(get(paste0(var)) ~ as.numeric(year))[[1]][2]) #linear regression slope
  bibi<-ddply(bibi, .variables =scal, mutate, pval = summary(lm(get(paste0(var)) ~ as.numeric(year)))[[4]][2,4]) #linear regression pvalue
  bibi$trend <- ifelse(bibi$pval < 0.05 & bibi$slope < 0, "negative", ifelse(bibi$pval < 0.05 & bibi$slope > 0, "positive", "none")) # interpret linear regression trend direction as positive, negative or none
  
  #Then Mann-kendall
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.tau=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[1]) #What is the strength of the trend?
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.Sen=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[2]) #what is the slope (points per year)
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.pval=kendallTrendTest( get(paste0(var))~ as.numeric(year))$p.value) #how significant is the pvalue for mannkendall?
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.inter=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[3]) #This is B0 or the intercept
  bibi$MKtrend <- ifelse(bibi$mk.pval < 0.05 & bibi$mk.tau < 0, "negative", ifelse(bibi$mk.pval < 0.05 & bibi$mk.tau > 0, "positive", "none")) #this is an interpretation of the pvalue and tau (if it is less that .05, what direction is the trend.)
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

#######

### code for calculating trends using only data in "bibi_no_late_starts", which is the most filtered (most conservative) set of data; see code above for all code calculating trends using all data in "bibi"

bibi_no_late_starts<-subset(bibi_no_late_starts, site_code!="")
form<-names(bibi_no_late_starts)[5:(ncol(bibi_no_late_starts)-3)]
# bibi.lr[bibi.lr$Year<2012,"Overall.Score"]<-(bibi.lr[bibi.lr$Year<2012,"Overall.Score"]*.972)+5.01
# for (i in 9:(ncol(bibi.lr))) {

#this loops through each of the 11 BiBi metrics and calculates the trend according to site, then records the output repeatedly for that site along all years. Why do it this way?

library(plyr) ## Figure out how to re-write this without the plyr package

for (var in form) {
  #First linear regression (slope, p-value, and interpretation)
  bibi_no_late_starts<- ddply(bibi_no_late_starts, .variables =scal, mutate, slope =lm(get(paste0(var)) ~ as.numeric(year))[[1]][2]) #linear regression slope
  bibi_no_late_starts<-ddply(bibi_no_late_starts, .variables =scal, mutate, pval = summary(lm(get(paste0(var)) ~ as.numeric(year)))[[4]][2,4]) #linear regression pvalue
  bibi_no_late_starts$trend <- ifelse(bibi_no_late_starts$pval < 0.05 & bibi_no_late_starts$slope < 0, "negative", ifelse(bibi_no_late_starts$pval < 0.05 & bibi_no_late_starts$slope > 0, "positive", "none")) # interpret linear regression trend direction as positive, negative or none
  
  #Then Mann-kendall
  bibi_no_late_starts<-ddply(bibi_no_late_starts,.(get(scal)),mutate, mk.tau=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[1]) #What is the strength of the trend?
  bibi_no_late_starts<-ddply(bibi_no_late_starts,.(get(scal)),mutate, mk.Sen=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[2]) #what is the slope (points per year)
  bibi_no_late_starts<-ddply(bibi_no_late_starts,.(get(scal)),mutate, mk.pval=kendallTrendTest( get(paste0(var))~ as.numeric(year))$p.value) #how significant is the pvalue for mannkendall?
  bibi_no_late_starts<-ddply(bibi_no_late_starts,.(get(scal)),mutate, mk.inter=kendallTrendTest( get(paste0(var))~ as.numeric(year))$estimate[3]) #This is B0 or the intercept
  bibi_no_late_starts$MKtrend <- ifelse(bibi_no_late_starts$mk.pval < 0.05 & bibi_no_late_starts$mk.tau < 0, "negative", ifelse(bibi_no_late_starts$mk.pval < 0.05 & bibi_no_late_starts$mk.tau > 0, "positive", "none")) #this is an inpterpretation of the pvalue and tau (if it is less that .05, what direction is the trend.)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-7]<-paste0("slope_",var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-6]<-paste0("LRpval_", var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-5]<-paste0("LRtrend_", var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-4]<-paste0("MKtau_", var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-3]<-paste0("MKSlope_", var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-2]<-paste0("MKpval_", var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)-1]<-paste0("MKinter_", var)
  colnames(bibi_no_late_starts)[ncol(bibi_no_late_starts)]<-paste0("MKtrend_", var)
}

#This is a creative way of showing both the duration of the data and the overall slope.
sampletimes<- ggplot(bibi_no_late_starts, 
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
sampletimes
# ggsave("sampletimes.png", plot=sampletimes, width=20, height=10)

names(bibi_no_late_starts)

##This step gets one trend per site.
overall_trends <- bibi_no_late_starts |> 
  select(site_code, subbasin, timerange:(ncol(bibi_no_late_starts))) |> 
  unique()

# overall_trends<-unique(subset(overall_trends, select=c("Site.Code" , "Subbasin", "Stream.or.River", str_subset(colnames(overall_trends), "MKtau"), str_subset(colnames(overall_trends), "MKpval"),str_subset(colnames(overall_trends), "MKtrend"))))
write.csv(overall_trends, paste0(output_folder,"/Site_Trends_coarse_PSSB_no_late_starts.csv"))

bibi_sum_site<-ddply(bibi_no_late_starts, "site_code", summarise, mean_Overall.Score=mean(overall_score), SD_overallscore=sd(overall_score), SE_overallscore=sd(overall_score)/sqrt(length(overall_score)), len=length(overall_score), MK_trend=unique(MKtrend_overall_score), Avg_MKSlope_Overall.Score=mean(MKSlope_overall_score), Avg_MKtau_Overall.Score=mean(MKtau_overall_score), Avg_MKpval_Overall.Score=mean(MKpval_overall_score), sites=length(unique(site_code)), min=min(overall_score), max=max(overall_score), med=median(overall_score))
write.csv(bibi_sum_site,paste0(output_folder,"/score_MKtrend_summarized_by_site_PSSB_no_late_starts.csv"))

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


## Here it appears beth is trying to make a plot for each site, where the shape indicates the significance of the trend, the line is the trend, and the colors indicate general stream health. 
library(RColorBrewer)

# for (var in colnames(bibi[str_subset(colnames(bibi), "trend")])) {
# if(nrow(bibi[bibi[var]=='positive',])>=1) {
var="MKtrend_overall_score"

bibi_small <- bibi_no_late_starts|> filter_out(site_code %in% later_starts)

bibi_small$pchOS<-as.factor(bibi_small$MKtrend_overall_score)
shapes = c(6, 4, 2) 
shapes <- shapes[as.numeric(bibi_small$pchOS)]


pp<-ggplot(bibi_small, aes(year_num, get(str_subset(colnames(bibi_small[,c(1:18)]), str_split(var, "_", n=2)[[1]][2])), group = get(scal), fill = get(scal), colour = get(scal))) + 
  ylim(c(0,100)) + 
  geom_point( aes(shape = c(shapes))) + scale_shape_identity()+ #theme(legend.position="none") +
  facet_wrap(~ subbasin, ncol = 4) +  #factor(stream2) + 
  theme(legend.position="none",axis.title=element_text(size=20,face="bold"), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90), strip.text = element_text(size=12)) + labs(y = var, x = "Year")
pp
pp<-pp+geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])),  colour=get(scal)))+
  scale_colour_manual(values=rep(brewer.pal(8,"Dark2"), times=30))
# pp<-pp+geom_text(mapping= aes(x=2012, y=6, label=paste("Slope = ", round(get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])), 3))), colour="black")+theme(text=element_text(size=20), plot.title=element_text(hjust=.5, size=24, face="bold"))
# ggsave(paste0("WRIA9_BIBI~timeUp_", scal, "_", var,   "-2.png"), plot = pp, width=20, height=10)
# }
pp
# }
ggsave(paste0(output_folder,"/OverallScores_Sites_PSSB_plots_by_basin.png"), plot = pp, width=20, height=10)


all_sites_facet<- ggplot(bibi_small, aes(x=year_num, y=overall_score))+  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  # geom_smooth(se=F, method="lm")+
  geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2]))), color="blue", size=1)+
  theme(legend.position = "none")+
  facet_wrap(~site_code)+geom_point( aes(shape = c(shapes))) + scale_shape_identity()+labs(y="Overall B-IBI Score")
all_sites_facet
ggsave(paste0(output_folder,"/OverallScores_Sites_PSSB_plots_by_site.png"), plot = all_sites_facet, width=20, height=10)

#for graph of Duwamish trib - declining site
site09DUW0091 <- bibi_small |> 
  filter(site_code == "09DUW0091")

plot_09DUW0091<- ggplot(site09DUW0091, aes(x=year_num, y=overall_score))+  
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_point(aes(shape = MKtrend_overall_score, fill = MKtrend_overall_score), size  = 5) +
  scale_shape_manual(values=c(21,24,25), breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  scale_fill_manual(values=c("white", "gray", "black"),breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  geom_smooth(method="lm", formula = y ~ x, se=F, color = "black" )+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  labs(
    x ="Year",
    y = "BIBI Score") +
  theme_bw()+
  theme(
    text=element_text(size=25), 
    legend.position = "none",
  )
plot_09DUW0091
ggsave(paste0(output_folder,"/08WES1340_trend.png"), plot = plot_08WES1340, width=20, height=10)


#for graph of Taylor Creek 08WES1340 [need to figure out how to loop this for all sites]

site08WES1340 <- bibi_small |> 
  filter(site_code == "08WES1340")

plot_08WES1340<- ggplot(site08WES1340, aes(x=year_num, y=overall_score))+  
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_point(aes(shape = MKtrend_overall_score, fill = MKtrend_overall_score), size  = 5) +
  scale_shape_manual(values=c(21,24,25), breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  scale_fill_manual(values=c("white", "gray", "black"),breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  geom_smooth(method="lm", formula = y ~ x, se=F, color = "black" )+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  labs(
    x ="Year",
    y = "BIBI Score") +
  theme_bw()+
  theme(
    text=element_text(size=25), 
    legend.position = "none",
    )
plot_08WES1340
ggsave(paste0(output_folder,"/08WES1340_trend.png"), plot = plot_08WES1340, width=20, height=10)


site08WES1490 <- bibi_small |> 
  filter(site_code == "08WES1490")
plot_08WES1490<- ggplot(site08WES1490, aes(x=year_num, y=overall_score))+  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_point(aes(shape = MKtrend_overall_score, fill = MKtrend_overall_score), size  = 6) +
  scale_shape_manual(values=c(21,24,25), breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  scale_fill_manual(values=c("white", "gray", "black"),breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  geom_smooth(method="lm", formula = y ~ x, se=F, color = "black" )+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  labs(
    x ="Year",
    y = "BIBI Score") +
  theme_bw()+
  theme(
    text=element_text(size=25), 
    legend.position = "none")
plot_08WES1490
ggsave(paste0(output_folder,"/08WES1490_trend.png"), plot = plot_08WES1490, width=20, height=10)

# Can't figure out how to plot two lines (by site) on one graph; trying to plot 08WES1340 and 08WES1490 together
sites08WES_1340_1490 <- bibi_small |> 
  filter(site_code %in% c("08WES1490", "08WES1340"))
plot_08WES_1340_1490<-ggplot(sites08WES_1340_1490, aes(x=year_num, y=overall_score, color = site_code))+  
  geom_rect(xmin=2001, xmax=2026, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2026, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_point(aes(shape = MKtrend_overall_score, fill = MKtrend_overall_score), size  = 5) +
  scale_shape_manual(values=c(21,24,25), breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  scale_fill_manual(values=c("white", "darkgray", "black"),breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  geom_smooth(method="lm", formula = y ~ x, se=F, color = "black")+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  labs(
    x ="Year",
    y = "BIBI Score") +
  theme(text=element_text(size=25)) #removed "legend.position = "none"
plot_08WES_1340_1490


## code for from FWS paper - pasted here just as an example for testing/borrowing code
imp_change_02_25 <-ggplot(imp_change, aes(x=change_2002_2025, y=avg_mk_slope_overall_score))+
  geom_point(aes(shape=mk_trend,  fill=mk_trend)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black")+
  stat_poly_eq(formula = y ~ x, rr.digits = 2, coef.digits = 3,label.x=5, label.y=-1, position="identity", color = "black", aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "*\", \"*")), parse = TRUE)+
  scale_shape_manual(values=c(21,24, 25), breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  scale_fill_manual(values=c("white", "darkgray", "black"),breaks=c("none","positive", "negative"), labels=c("non-significant", "increasing", "decreasing"))+
  geom_line(y=0, lty="dashed") +
  theme_bw()+
  theme(text=element_text(family="sans"), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))+
  xlab("Change in % impervious cover") +
  ylab("Slope estimate (BIBI trend)") +
  labs(shape="Trend Direction", fill="Trend Direction")+
  guides(fill = "none", shape="none") 
imp_change_02_25



# summary table of sites that were found to be trending using MK tests (negative, none, or positive); includes overall score and metrics 
bibi_small %>%
  select(site_code, 
        MKtrend_overall_score, 
        MKtrend_taxa_richness_score,
        
        MKtrend_ephemeroptera_richness_score,
        MKtrend_plecoptera_richness_score,
        MKtrend_trichoptera_richness_score,

        MKtrend_clinger_richness_score,
        MKtrend_long_lived_richness_score,
        MKtrend_intolerant_richness_score,
        
        MKtrend_percent_dominant_score,
        MKtrend_predator_percent_score,
        MKtrend_tolerant_percent_score) |> 
  unique () |> 
  pivot_longer(!site_code, names_to = "metric", values_to = "trend") |> 
  dplyr::count(metric,trend) |>
  pivot_wider(names_from = trend, values_from = n )



#CART analysis ####
# Load the library
library(rpart)

year_overall <-bibi_small |> 
  select(site_code, year_num, overall_score) |> 
  unique()

# CART model using only one predictor variable: Year
# method = "anova" for continuous variable (?)
cart_model <- rpart(
  overall_score ~ year_num,
  data = year_overall,
  method = "anova",
  cp = 0.001) # control parameters

# Print model summary
print(cart_model)

# Plot the tree
plot(cart_model, uniform = TRUE, main = "CART: overall_score ~ Year")
text(cart_model, use.n = TRUE, all = TRUE, cex = 0.8)

