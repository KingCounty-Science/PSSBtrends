
library(plyr)
library(openxlsx)
library(lubridate)
library(stringr)
library(tidyverse)
library(dplyr)
getwd()

atts<-read.xlsx("G:/GreenWQA/Biota/Contracts/Bellevue/2023 Report/for trends analysis/2012_taxa_attributes.xlsx")
lookup<-read.csv("./Inputs/ORWA_TaxaTranslator_20240417.csv")##update this to latest translator table.
BCG_atts<-read.csv("./Inputs/ORWA_Attributes_20240417.csv")#load the BCG taxonomic hierarchy

##This function takes the raw taxa data .txt output from PSSB and binds it all into one data object
taxaBind <- function(file.path) {

  path.files <- list.files(file.path)
  # read in files
  list.with.each.file <- lapply(paste(file.path, list.files(file.path), sep = ''), function(y) read.delim(y, header=TRUE))
  taxa<-do.call("rbind.data.frame", list.with.each.file)
  return(taxa)


}

file.path="./Inputs/taxonomy data/"
raw<-taxaBind(file.path)

length(unique(raw$Project))
length(unique(raw$Agency))
length(unique(raw$Site.Code))
length(unique(raw$WRIA.Number))
length(unique(raw$Stream.or.River))
length(unique(raw$Subbasin))

###Kate has identified what sites she wants to run trends on. This subsets the taxa data to just those sites
# site<-read.xlsx("./Inputs/ScoresByYear_all streams and rivers default selection.xlsx", detectDates = T, sheet="NEW - sites for trends")
# raw2<-subset(raw, Site.Code %in% site$Site.Code)

####Fix some names to match the translator better
raw[raw$Taxon=="Lepidotoma-panel case larvae","Taxon"]<-"Lepidostoma-panel case larvae" ###fix this in PSSB

PSSB_taxa<-unique(raw[,c(28, 29, 48:69)])
##there are some repeat entries that somewhere in the hierarchy have an NA instead of "". This yields multiples of the same taxa. Fix this.
PSSB_taxa[is.na(PSSB_taxa)]<-""
PSSB_taxa<-unique(PSSB_taxa) #we're generating a list of all taxa in PSSB samples
any(duplicated(PSSB_taxa$Taxon.Serial.Number))
any(duplicated(PSSB_taxa$Taxon))
PSSB_taxa[which(duplicated(PSSB_taxa$Taxon)),] ### two entries in PSSB-- fix this in PSSB
PSSB_taxa<-subset(PSSB_taxa, Taxon.Serial.Number !="-40")

####merge the translator lookup with the raw data, rename OTU_MetricCalc, and look for any taxa missing a translation
OTU<-merge(raw, subset(lookup, select=c(Taxon, OTU_MetricCalc, NonTarget)), by.x="Taxon", by.y="Taxon", all.x=T)
OTU[which(is.na(OTU$OTU_MetricCalc)),]
colnames(OTU)[ncol(OTU)-1]<-"OTU"
missing<-unique(OTU[which(is.na(OTU$OTU)), "Taxon"])## screening step to see if any taxa aren't mapped

###prepare the data
OTU$Visit.Date<-format(as.Date(OTU$Visit.Date, "%m/%d/%Y"))
OTU$Year<-year(OTU$Visit.Date)

OTU<-subset(OTU, is.na(QC.Replicate.Of.Sample.Code)|QC.Replicate.Of.Sample.Code=="")##remove QC replicates
OTU<-subset(OTU, Non.B.IBI=="False")##remove non-target organisms
OTU[which(OTU$OTU=="DNI"),"OTU"]<-OTU[which(OTU$OTU=="DNI"),"Taxon"]###These are marked as "DNI" in BCG translation table, but they aren't on B-IBI exclusion list. Adding back in for now.

OTU$Unique<-as.logical(OTU$Unique)

##collapse to Visit.ID, because 1998-2015 samples were often three reps of 3 sq ft with different sample names for each rep
OTU_collapsed<-ddply(OTU, .(Visit.ID, OTU, WRIA.Number, Agency, Basin, Subbasin, Stream.or.River, Project, Visit.Date, Year, Latitude, Longitude, Lab.Name, Site.Code), summarize, Quantity_OTU = sum(Quantity), Unique_OTU=any(Unique))

OTU_collapsed$Visit.Date<-as.Date(OTU_collapsed$Visit.Date, "%Y-%m-%d")

########create lookup table for taxa hierarchy. ####
#we need to  get the BCG hierarchy and the PSSB taxa hierarchy in the same format and combine them

names(BCG_atts)
names(BCG_atts[,c(1, 22,24:40)])
hierarchy<-BCG_atts[,c(1, 22,24:40)]
names(PSSB_taxa)
names(hierarchy)
##the BCG hierarchy isn't the same as PSSB. Need to consolidate some levels into Species Group, then rename columns
hierarchy[which(hierarchy$SpeciesComplex!=""&hierarchy$SpeciesGroup==""),"SpeciesGroup"]<-hierarchy[which(hierarchy$SpeciesComplex!=""&hierarchy$SpeciesGroup==""),"SpeciesComplex"]
hierarchy[which(hierarchy$SpeciesSubGroup!=""&hierarchy$SpeciesGroup==""),"SpeciesGroup"]<-hierarchy[which(hierarchy$SpeciesSubGroup!=""&hierarchy$SpeciesGroup==""),"SpeciesSubGroup"]
hierarchy$Superclass<-NA
hierarchy$Infraclass<-NA
hierarchy$Superorder<-NA
hierarchy$Infraorder<-NA
hierarchy$Custom.Subfamily<-NA
hierarchy$Subtribe<-NA
hierarchy$Subspecies<-NA

##need to get the hierarchy in the right order
names(PSSB_taxa)
names(hierarchy)
names(hierarchy[,c(2,1, 3, 4, 20, 5:6,21,22,  7:8, 23,9:10, 11,24, 12,13,25,14,16, 15,19, 26)])
hierarchy<-hierarchy[,c(2,1, 3, 4, 20, 5:6,21,22,  7:8, 23,9:10, 11,24, 12,13,25,14,16, 15,19, 26)]

names(hierarchy)<-names(PSSB_taxa)
hierarchy[is.na(hierarchy)]<-""
PSSB_taxa$Taxon.Serial.Number<-as.character(PSSB_taxa$Taxon.Serial.Number)

append<-setdiff(hierarchy$Taxon, PSSB_taxa$Taxon)### look for taxa that are in the BCG hierarchy that aren't in the PSSB hierarchy
append<-hierarchy[hierarchy$Taxon %in% append,] ###restrict the BCG hierarchy to just those not in the PSSB hierarchy
any(duplicated(PSSB_taxa$Taxon))
PSSB_taxa[which(duplicated(PSSB_taxa$Taxon)),] 

hier_combined<-rbind(PSSB_taxa, append) ###combine
any(duplicated(hier_combined$Taxon))

# ##Append correct hierarchy to taxa
OTU_collapsed<-merge(OTU_collapsed, hier_combined, by.x="OTU", by.y="Taxon", all.x=T)
any(is.na(OTU_collapsed$Phylum))
OTU_collapsed[which(is.na(OTU_collapsed$Phylum)),]
any(OTU_collapsed$Phylum=="")
###################################### Roll-up by broad rules ##########

KC_taxa_coarse<-OTU_collapsed

KC_taxa_coarse$OTU_COARSE<-""
names(KC_taxa_coarse)

coarse_rules<-data.frame(taxa=c("Oligochaeta", "Acari", "Gastropoda","Dytiscidae", "Simuliidae", "Chironomidae", "Trichoptera"), ranktouse=c("Subclass", "Subclass", "Family", "Family", "Family", "Family", "Genus"), rank=c("Subclass", "Subclass", "Class", "Family", "Family", "Family", "Order"))

for (j in 1:nrow(coarse_rules)){
  
  STE_rank<-coarse_rules[j, "ranktouse"]
  rank<-coarse_rules[j, "rank"]
  taxa<-coarse_rules[j, "taxa"]
  index<-which(names(KC_taxa_coarse)== STE_rank)
  rankindex<-which(names(KC_taxa_coarse)== rank)
  halt1<-which(names(KC_taxa_coarse)== "Phylum")
  halt2<-which(names(KC_taxa_coarse)== "Species")
  
  for (i in halt2:halt1){
    if(i > index) {
      KC_taxa_coarse[which(KC_taxa_coarse[,rankindex]==taxa),i]<-""
    } else if (i<= index) {
      
      KC_taxa_coarse[which(KC_taxa_coarse[,rankindex]==taxa&KC_taxa_coarse$OTU_COARSE==""),"OTU_COARSE"]<-KC_taxa_coarse[which(KC_taxa_coarse[,rankindex]==taxa&KC_taxa_coarse$OTU_COARSE==""),i]
    }
    
  }
}
KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==""),"OTU_COARSE"]<-KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==""),"OTU"]

KC_taxa_coarse[is.na(KC_taxa_coarse)]<-""

################read in PSSB attribute table, do rolling lookup between coarse taxa hierarchy and attribute table ############

names(KC_taxa_coarse)
missing_atts<-unique(subset(KC_taxa_coarse, select=c(18:40)))
attribs2<-data.frame(Taxon.Name=character(), Predator=character(), Long.Lived=character(), Tolerant=character(), Intolerant=character(), Clinger=character(), OTU_COARSE=character(),  iter=numeric())

for (i in 1:ncol(missing_atts)){
  k<-(ncol(missing_atts)+1)-i
  attribs<-merge(subset(atts, select=c("Taxon.Name", "2012.Clinger", "2012.Intolerant", "2012.Long.Lived", "2012.Predator", "2012.Tolerant")), missing_atts[, c(k, 23)], by.x="Taxon.Name", by.y=names(missing_atts[k]))
  names(attribs)<-str_replace(names(attribs), "2012.", "")
  names(attribs)<-str_replace(names(attribs), ".1", "")
  attribs2<-rbind(attribs, attribs2)
  names(attribs2)<-str_replace(names(attribs2), ".1", "")
  missing_atts<-subset(missing_atts, !OTU_COARSE %in% attribs2$OTU_COARSE)
}

attribs2<-attribs2[, -1]
names(attribs2)[6]<-"Taxon"
atts<-unique(attribs2)
any(duplicated(atts$Taxon))
atts[(which(duplicated(atts$Taxon))),]
missing<-unique(KC_taxa_coarse$OTU_COARSE)[!unique(KC_taxa_coarse$OTU_COARSE) %in% atts$Taxon] ##These taxa do not have a match in the PSSB attribute table

library(plyr)
OTU_collapsed2<-ddply(KC_taxa_coarse, .(Visit.ID, OTU_COARSE,Agency, WRIA.Number, Basin, 
                                        Subbasin, Stream.or.River, Project, Visit.Date, 
                                        Year, Latitude, Longitude, Lab.Name, Site.Code
), summarize, Quantity_OTU = sum(Quantity_OTU), Unique_OTU=any(Unique_OTU))

##append attributes from the lookup table we made
OTU_collapsed3<-left_join(OTU_collapsed2, hier, by=c("OTU_COARSE"="Taxon"))
any(is.na(OTU_collapsed3$Clinger))
##Fill in attributes as "FALSE" for the taxa with no attribute matches
OTU_collapsed3[which(is.na(OTU_collapsed3$Clinger)),"Clinger"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Intolerant)),"Intolerant"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Long.Lived)),"Long.Lived"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Predator)),"Predator"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Tolerant)),"Tolerant"]<-FALSE

##append on hierarchy
###some taxa have multiple unique hierarchies because different entries got rolled up to the same level, and taxonomists have been uneven in entering in 'infraorder', 'suborder' and 'infraclass'. Need to run some loops to clean this up by selecting the most complete hierarchy available
new_hierarchy<-unique(KC_taxa_coarse[,c(18:40)])
countlength<-ddply(new_hierarchy, .(OTU_COARSE), summarize, count=length(OTU_COARSE))
fixthese<-countlength[countlength$count>1,]
check<-new_hierarchy[new_hierarchy$OTU_COARSE %in% fixthese$OTU_COARSE,]
for (u in 1:nrow(check)){
  check$sum[u]<-sum(check[u,]!="")
}
for (j in unique(check$OTU_COARSE)) {
  
  test<-check[check$OTU_COARSE==j,]
    summ<-colSums(test == "")
    ilen<-max(summ)
    fixlevels<-summ[which(summ<ilen& summ>0)]
    
     for (z in names(fixlevels)){
       update<-test[which(test[,z]!=""), z]
       if (length(update)>1) print("Error in taxonomy agreement")
       check[check$OTU_COARSE==j,z]<-update
     }
}
check<-unique(check)
check<-check[,-which(names(check) %in% "sum")]

new_hierarchy<-subset(new_hierarchy, ! OTU_COARSE %in% check$OTU_COARSE)
new_hierarchy<-rbind(new_hierarchy, check)
countlength<-ddply(new_hierarchy, .(OTU_COARSE), summarize, count=length(OTU_COARSE))
fixthese<-countlength[countlength$count>1,]

OTU_collapsed3<-merge(OTU_collapsed3, new_hierarchy, by.x="OTU_COARSE", by.y="OTU_COARSE", all.x=T)
any(is.na(OTU_collapsed3$Phylum))
unique(OTU_collapsed3[which(is.na(OTU_collapsed3$Phylum)),]$OTU_COARSE)
any(OTU_collapsed3$Phylum=="")

write.csv(OTU_collapsed3, "Collapsed_Coarse_Taxa.csv")
OTU_collapsed3<-read.csv( "Collapsed_Coarse_Taxa.csv")

##The rarify function below is similar to how subsampling was done for the trends report.
rarify<- function (inbug, sample.ID, abund, subsize, taxa, mySeed = NA)
{ #set.seed(mySeed)
  KC_rarified<-inbug[0,]
  for(visit in unique(inbug[,sample.ID])){
    set.seed(mySeed, "Mersenne-Twister",normal.kind = "Inversion",  sample.kind="Rounding") ##5/15/2024, BAS moved setseed here so that scores are repeatable when dataset changes. Discovered quirk of set.seed that made results inconsistent when dataset was updated to include additional years of data. 7/9/2024 BAS added additional parameters so that results are consistent between versions of R.
    print(visit)
    test<-subset(inbug, get(sample.ID)==visit)
    testsample<-rep(test[,taxa], test[,abund])
    if(sum(test[,abund])>=subsize){
      subsamp<-sample(x = testsample, size = subsize,replace=F)
    }
    else {subsamp<-sample(x = testsample, size = sum(test[,abund]),replace=F)}
    subsamp<-as.data.frame(table(subsamp))
    names(subsamp)<-c(taxa,abund)
    subsamp_meta<-test[c(match(subsamp[,taxa], test[,taxa])),]
    subsamp_meta<-subset(subsamp_meta, select=-c(get(abund)))
    subsamp_comp<-merge(subsamp_meta, subsamp, by=taxa)
    KC_rarified<-rbind(subsamp_comp, KC_rarified)
    
  }
  return(KC_rarified)
}

KC_rarified<-rarify(inbug=OTU_collapsed3,
                    sample.ID="Visit.ID",
                    abund="Quantity_OTU",
                    subsize<-450,
                    mySeed=17760704,
                    taxa="OTU_COARSE")



detach("package:plyr", unload = TRUE)

str(KC_rarified)
KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
KC_rarified$Predator<-as.logical(KC_rarified$Predator)


Tot_Richness<-KC_rarified %>% dplyr::filter(Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Total_Richness=length(OTU_COARSE))
E_Richness<-KC_rarified %>% dplyr::filter(Order=="Ephemeroptera", Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(E_Richness=length(OTU_COARSE))
P_Richness<-KC_rarified %>% dplyr::filter(Order=="Plecoptera", Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(P_Richness=length(OTU_COARSE))
T_Richness<-KC_rarified %>% dplyr::filter(Order=="Trichoptera", Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(T_Richness=length(OTU_COARSE))
Cling_Richness<-KC_rarified %>% dplyr::filter(`Clinger`==TRUE, Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Clin_Richness=length(OTU_COARSE))
LL_Richness<-KC_rarified %>% dplyr::filter(`Long.Lived`==TRUE, Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(LL_Richness=length(OTU_COARSE))
Intol_Richness<-KC_rarified %>% dplyr::filter(`Intolerant`==TRUE, Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Intol_Richness=length(OTU_COARSE))
Tot_Abund<-KC_rarified   %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Tot_Abund=sum(Quantity_OTU)) 
Tol_Abund<-KC_rarified   %>% dplyr::group_by(Visit.ID)%>% dplyr::filter(`Tolerant`==TRUE) %>% dplyr::summarise(Tol_Abund=sum(Quantity_OTU))
Pred_Abund<-KC_rarified   %>% dplyr::group_by(Visit.ID)%>% dplyr::filter(`Predator`==TRUE) %>% dplyr::summarise(Pred_Abund=sum(Quantity_OTU))
Dom_abund<-KC_rarified   %>% dplyr::group_by(Visit.ID) %>% dplyr::arrange(desc(Quantity_OTU)) %>% dplyr::slice(1:3) %>% dplyr::summarise(Dom_Abund=sum(Quantity_OTU))

KC_results<-Reduce(function(x, y, ...) merge(x,y, all=TRUE, ...), list(Tot_Richness, E_Richness, P_Richness, T_Richness, Cling_Richness, LL_Richness, Intol_Richness, Dom_abund,
                                                                       Tol_Abund, Pred_Abund,Tot_Abund))
KC_results<-KC_results %>% mutate(Tol_Percent=Tol_Abund/Tot_Abund*100, Pred_Percent= Pred_Abund/Tot_Abund*100, Dom_Percent= Dom_Abund/Tot_Abund*100)
KC_results[is.na(KC_results)]<-0

KC_results<-KC_results %>% mutate(Tot_Richness_Score= 10 * (Total_Richness-16)/(37-16))
KC_results<-KC_results %>% mutate(E_Richness_Score= 10 * (E_Richness-1)/(8-1))
KC_results<-KC_results %>% mutate(P_Richness_Score= 10 * (P_Richness-1)/(8-1))
KC_results<-KC_results %>% mutate(T_Richness_Score= 10 * (T_Richness-1)/(9-1))
KC_results<-KC_results %>% mutate(Cling_Richness_Score= 10 * (Clin_Richness-5)/(22-5))
KC_results<-KC_results %>% mutate(LL_Richness_Score= 10 * (LL_Richness-2)/(10-2))
KC_results<-KC_results %>% mutate(Intol_Richness_Score= 10 * (Intol_Richness-0)/(7-0))
KC_results<-KC_results %>% mutate(Dom_Percent_Score= 10-(10 * (Dom_Percent-44)/(82-44)))
KC_results<-KC_results %>% mutate(Pred_Percent_Score= 10 * (Pred_Percent-1)/(21-1))
KC_results<-KC_results %>% mutate(Tol_Percent_Score= 10-(10 * (Tol_Percent-0)/(43-0)))

KC_results[which(KC_results$Tot_Richness_Score>10),"Tot_Richness_Score"]<-10
KC_results[which(KC_results$Tot_Richness_Score<0),"Tot_Richness_Score"]<-0
KC_results[which(KC_results$E_Richness_Score>10),"E_Richness_Score"]<-10
KC_results[which(KC_results$E_Richness_Score<0),"E_Richness_Score"]<-0
KC_results[which(KC_results$P_Richness_Score>10),"P_Richness_Score"]<-10
KC_results[which(KC_results$P_Richness_Score<0),"P_Richness_Score"]<-0
KC_results[which(KC_results$T_Richness_Score>10),"T_Richness_Score"]<-10
KC_results[which(KC_results$T_Richness_Score<0),"T_Richness_Score"]<-0
KC_results[which(KC_results$Cling_Richness_Score>10),"Cling_Richness_Score"]<-10
KC_results[which(KC_results$Cling_Richness_Score<0),"Cling_Richness_Score"]<-0
KC_results[which(KC_results$LL_Richness_Score>10),"LL_Richness_Score"]<-10
KC_results[which(KC_results$LL_Richness_Score<0),"LL_Richness_Score"]<-0
KC_results[which(KC_results$Intol_Richness_Score>10),"Intol_Richness_Score"]<-10
KC_results[which(KC_results$Intol_Richness_Score<0),"Intol_Richness_Score"]<-0
KC_results[which(KC_results$Dom_Percent_Score>10),"Dom_Percent_Score"]<-10
KC_results[which(KC_results$Dom_Percent_Score<0),"Dom_Percent_Score"]<-0
KC_results[which(KC_results$Pred_Percent_Score>10),"Pred_Percent_Score"]<-10
KC_results[which(KC_results$Pred_Percent_Score<0),"Pred_Percent_Score"]<-0
KC_results[which(KC_results$Tol_Percent_Score>10),"Tol_Percent_Score"]<-10
KC_results[which(KC_results$Tol_Percent_Score<0),"Tol_Percent_Score"]<-0

KC_results<-KC_results %>% mutate(Overall.Score=Tot_Richness_Score+ E_Richness_Score+P_Richness_Score+T_Richness_Score+ Cling_Richness_Score+ LL_Richness_Score+ Intol_Richness_Score+ Dom_Percent_Score+ Pred_Percent_Score+ Tol_Percent_Score)

KC_results<-left_join(KC_results, unique(subset(KC_rarified, select=c("Visit.ID","Agency", "WRIA.Number", "Basin", 
                                                                      "Subbasin", "Stream.or.River", "Project", "Visit.Date", 
                                                                      "Year", "Lab.Name", "Site.Code"))), by="Visit.ID")

KC_results<-left_join(KC_results, unique(subset(raw, select=c(Visit.ID, Latitude, Longitude))), by="Visit.ID")

write.csv(KC_results, "B-IBI_results_PSSB_Scores.csv")


###############


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
# bibi<-KC_results
unique(bibi$Site.Code)
bibi <-read.csv("B-IBI_results_PSSB_Scores.csv")
bibi<-subset(bibi, Tot_Abund>=450)
unique(bibi$Site.Code)

#########  Format, remove data from wrong months, average scores from samples with same site+year, and subset B-IBI data to data since 2001 with more than 9 years of data
bibi$month <- format(as.Date(bibi$Visit.Date, '%Y-%m-%d'), "%m") # add month column

bibi$sampnum <- paste(bibi$Trend_station, format(as.Date(bibi$Visit.Date, '%Y-%m-%d'), "%y"), sep = "_") # add sample identifier
bibi <- droplevels(bibi[bibi$month %in% c('07', '08', '09', '10'),]) # exclude samples collected outside of July - Oct
dups<-bibi[duplicated(bibi[,c("Project" , "WRIA.Number"  , "Basin" , "Subbasin" , "Stream.or.River" , "Site.Code", "Year")]),]

bibi <-aggregate(cbind( Overall.Score, Tot_Richness_Score, E_Richness_Score, P_Richness_Score, T_Richness_Score, Cling_Richness_Score, LL_Richness_Score, Intol_Richness_Score, Dom_Percent_Score,  Pred_Percent_Score, Tol_Percent_Score) ~   Site.Code + Subbasin +WRIA.Number+ Year, data = bibi, FUN = "mean", na.action = na.exclude)

bibi.lr<-bibi
bibi.lr1<-bibi.lr

ggplot(bibi.lr, aes(x=Site.Code, y=Year))+geom_tile()

ggplot(bibi.lr, aes(x=Year, y=Overall.Score))+geom_smooth(se=F)+theme(legend.position = "none")+
  geom_point()



library(ggpmisc)

meanscores<-ddply(bibi.lr1, .(Year), summarize, meanScore=mean(Overall.Score), nyears=length(unique(Site.Code)), sd=sd(Overall.Score), min=min(Overall.Score), max=max(Overall.Score), median=median(Overall.Score))
meanscore<-ggplot(meanscores, aes(x=Year, y=meanScore))+geom_smooth(method="lm",se=F)+theme(legend.position = "none")+
  geom_point()+ylim(c(0,100))+labs(y="BIBI Score")+
  stat_poly_eq(formula = y ~ x, rr.digits = 2, coef.digits = 3, size = 4, color = "black", aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), parse = TRUE) 

median(bibi.lr1$Overall.Score)
median(bibi.lr1$Year)

##calculating intercept for RKT derived Sen-Theil line for Overall Score (0.38; calculated later in the code), from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf pdf page 15, and https://rdrr.io/github/USGS-R/HASP/src/R/statistics.R
median(bibi.lr1$Overall.Score)-(0.38*median(bibi.lr1$Year))


scores<-ggplot(bibi.lr1 )+
  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=Overall.Score, x=as.numeric(paste0(Year)), group=as.factor(Year)), fill="white")+
  geom_point(data=meanscores, mapping=aes(x=as.numeric(paste0(Year)), y=meanScore))+
  # geom_smooth(data=meanscores, mapping=aes(x=as.numeric(paste0(Year)), y=meanScore), method="lm",se=F)+
  labs(y="Overall Score", x="Year")+ theme(text = element_text(size = 28))+
  # geom_quantile(data=bibi.lr1, mapping=aes(x=as.numeric(paste0(Year)), y=Overall.Score), quantiles=.5, color="red", size=1)#+
geom_abline(slope=0.38, intercept = -713.7596, color="red", lwd=1)##adding Sen-Theil line and intercept calculated above


names(bibi.lr1)

l=unique(c(as.character(bibi.lr1$Site.Code)))
bibi.lr1 <- bibi.lr1[complete.cases(bibi.lr1[, 1]), ]
str(bibi.lr1)
bibi.lr1$Site.Number<-as.numeric(factor(bibi.lr1$Site.Code, levels=l))

metrics<-colnames(bibi.lr1)[5:(ncol(bibi.lr1)-1)]

library(rkt)

for (each in metrics) {
  bibi.lr1$rkt.pval <- with(bibi.lr1, rkt(date = Year, y = get(each), block = Site.Number, correct = FALSE)[1])
  bibi.lr1$rkt.Sen <- with(bibi.lr1, rkt(date = Year, y = get(each), block = Site.Number, correct = FALSE)[3])
  bibi.lr1$rkt.tau <- with(bibi.lr1, rkt(date = Year, y = get(each), block = Site.Number, correct = FALSE)[12])
  bibi.lr1$RKTtrend <- with(bibi.lr1, ifelse(rkt.pval < 0.05 & rkt.tau < 0, "negative", ifelse(rkt.pval < 0.05 & rkt.tau > 0, "positive", "none")))
  colnames(bibi.lr1)[ncol(bibi.lr1)-3]<-paste0("rkt.pval_", each)
  colnames(bibi.lr1)[ncol(bibi.lr1)-2]<-paste0("rkt.Sen_", each)
  colnames(bibi.lr1)[ncol(bibi.lr1)-1]<-paste0("rkt.tau_", each)
  colnames(bibi.lr1)[ncol(bibi.lr1)]<-paste0("RKTtrend_", each)
}

overall_trends<-unique(bibi.lr1[(ncol(bibi.lr)+2):(ncol(bibi.lr1))])
overall_trends = as.matrix(overall_trends)
write.csv(overall_trends,"RKT_overall_trends_PSSB.csv")
trend<-subset(overall_trends, select=c(str_subset(colnames(overall_trends), "RKTtrend_")))

##plot the trends for each metric
for (i in 5:15) {
  ylabel<-names(bibi.lr1[i])
  # p<-ggplot(bibi.lr1, aes(as.factor(year), bibi.lr1[,i]))+geom_violin(draw_quantiles=c(.25, .5, .75))+ stat_summary(fun.y=median, geom="point", size=2, color="red")+labs(y=ylabel)
  o<-ggplot(bibi.lr1, aes(as.factor(Year), bibi.lr1[,i]))+geom_boxplot()+labs(y=ylabel)#+geom_abline(aes(intercept=get(paste0("MKinter_",names(bibi.lr1[i]))),slope=get(paste0("MKSlope_",names(bibi.lr1[i])))))
  # print(p)
  print(o)
  ggsave(paste0("metrics_boxplots_", ylabel, "PSSB.png"), plot = o, width=20, height=10)
}


########## basin specific RKT regressions ####

bibi.lr2 <- ddply(bibi.lr, 'Subbasin', mutate, timerange = (as.numeric(max(Year))-as.numeric(min(Year))))
bibi.lr2 <- ddply(bibi.lr2, 'Subbasin', mutate, por = length(unique(Year))) # add on period of record col
bibi.lr2<-subset(bibi.lr2, Subbasin!="")

bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$por > 9,])
bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$timerange > 9,])# look only at streams with >9 years of data
ddply(bibi.lr2, 'Subbasin', summarize, count=length(unique(Year)))

unique(bibi.lr2$Site.Code)

bibi.lr2<-ddply(bibi.lr2, 'Subbasin', mutate, max = max(Year))

scal<-"Subbasin" # "Site.Code" or "Stream.or.River"

form<-names(bibi.lr2)[5:(ncol(bibi.lr2)-3)]
bibi.lr2$Year<-as.numeric(paste0(bibi.lr2$Year))


l=unique(c(as.character(bibi.lr2$Site.Code)))
bibi.lr2$Site.Number<-as.numeric(factor(bibi.lr2$Site.Code, levels=l))

for (var in form) {
  bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.tau=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$tau)
  bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.Sen=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$B)
  bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.pval=rkt( date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$sl)
  bibi.lr2$MKtrend <- ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau < 0, "negative", ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau > 0, "positive", "none"))
  colnames(bibi.lr2)[ncol(bibi.lr2)-3]<-paste0("RKTtau_", var)
  colnames(bibi.lr2)[ncol(bibi.lr2)-2]<-paste0("RKTSlope_", var)
  colnames(bibi.lr2)[ncol(bibi.lr2)-1]<-paste0("RKTpval_", var)
  colnames(bibi.lr2)[ncol(bibi.lr2)]<-paste0("RKTtrend_", var)
}

names(bibi.lr2)
overall_trends<-unique(bibi.lr2[c(3,21:(ncol(bibi.lr2)))])
overall_trends<-unique(subset(overall_trends, select=c("Subbasin" ,str_subset(colnames(overall_trends), "RKTtau"), str_subset(colnames(overall_trends), "RKTpval"),str_subset(colnames(overall_trends), "RKTtrend"))))



write.csv(overall_trends, "Subbasin_RKTtrends_PSSB.csv")


bibi_sum<-ddply(bibi.lr2, "Subbasin", summarise, mean_Overall.Score=mean(Overall.Score), SD_overallscore=sd(Overall.Score), SE_overallscore=sd(Overall.Score)/sqrt(length(Overall.Score)), len=length(Overall.Score), RKT_trend=unique(RKTtrend_Overall.Score), Avg_RKTSlope_Overall.Score=mean(RKTSlope_Overall.Score), Avg_RKTtau_Overall.Score=mean(RKTtau_Overall.Score), Avg_RKTpval_Overall.Score=mean(RKTpval_Overall.Score), sites=length(unique(Site.Code)), min=min(Overall.Score), max=max(Overall.Score), med=median(Overall.Score))
write.csv(bibi_sum,"score_RKTtrend_summarized_by_subbasin_PSSB.csv")

pp<-ggplot(bibi.lr2, aes(Year, Overall.Score, group=Year))+geom_boxplot(aes(fill=RKTtrend_Overall.Score))+facet_wrap(~Subbasin, ncol=4)+labs(y="B-IBI Score", x="Year")+ theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("SubbasinRKTtrend_PSSB.png"), plot = pp, width=20, height=10)



########################
########## watershed specific RKT regressions ####

bibi.lr2 <- ddply(bibi.lr, 'WRIA.Number', mutate, timerange = (as.numeric(max(Year))-as.numeric(min(Year))))
bibi.lr2 <- ddply(bibi.lr2, 'WRIA.Number', mutate, por = length(unique(Year))) # add on period of record col
bibi.lr2<-subset(bibi.lr2, WRIA.Number!="")

bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$por > 9,])
bibi.lr2<-droplevels(bibi.lr2[bibi.lr2$timerange > 9,])# look only at streams with >9 years of data
ddply(bibi.lr2, 'WRIA.Number', summarize, count=length(unique(Year)))

unique(bibi.lr2$Site.Code)

bibi.lr2<-ddply(bibi.lr2, 'WRIA.Number', mutate, max = max(Year))

scal<-"WRIA.Number" # "Site.Code" or "Stream.or.River"

form<-names(bibi.lr2)[5:(ncol(bibi.lr2)-3)]
bibi.lr2$Year<-as.numeric(paste0(bibi.lr2$Year))


l=unique(c(as.character(bibi.lr2$Site.Code)))
bibi.lr2$Site.Number<-as.numeric(factor(bibi.lr2$Site.Code, levels=l))

for (var in form) {
  bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.tau=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$tau)
  bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.Sen=rkt(date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$B)
  bibi.lr2<-ddply(bibi.lr2,.(get(scal)),mutate, rkt.pval=rkt( date = Year, y= get(paste0(var)), block = Site.Number, correct = FALSE )$sl)
  bibi.lr2$MKtrend <- ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau < 0, "negative", ifelse(bibi.lr2$rkt.pval < 0.05 & bibi.lr2$rkt.tau > 0, "positive", "none"))
  colnames(bibi.lr2)[ncol(bibi.lr2)-3]<-paste0("RKTtau_", var)
  colnames(bibi.lr2)[ncol(bibi.lr2)-2]<-paste0("RKTSlope_", var)
  colnames(bibi.lr2)[ncol(bibi.lr2)-1]<-paste0("RKTpval_", var)
  colnames(bibi.lr2)[ncol(bibi.lr2)]<-paste0("RKTtrend_", var)
}

names(bibi.lr2)
overall_trends<-unique(bibi.lr2[c(4,21:(ncol(bibi.lr2)))])
overall_trends<-unique(subset(overall_trends, select=c("WRIA.Number" ,str_subset(colnames(overall_trends), "RKTtau"), str_subset(colnames(overall_trends), "RKTpval"),str_subset(colnames(overall_trends), "RKTtrend"))))



write.csv(overall_trends, "Trend_watershed_RKTtrends_PSSB.csv")


bibi_sum<-ddply(bibi.lr2, "WRIA.Number", summarise, mean_Overall.Score=mean(Overall.Score), SD_overallscore=sd(Overall.Score), SE_overallscore=sd(Overall.Score)/sqrt(length(Overall.Score)), len=length(Overall.Score), RKT_trend=unique(RKTtrend_Overall.Score), Avg_RKTSlope_Overall.Score=mean(RKTSlope_Overall.Score), Avg_RKTtau_Overall.Score=mean(RKTtau_Overall.Score), Avg_RKTpval_Overall.Score=mean(RKTpval_Overall.Score), sites=length(unique(Site.Code)), min=min(Overall.Score), max=max(Overall.Score), med=median(Overall.Score))
write.csv(bibi_sum,"score_RKTtrend_summarized_by_watershed_PSSB.csv")

pp<-ggplot(bibi.lr2, aes(Year, Overall.Score, group=Year))+geom_boxplot(aes(fill=RKTtrend_Overall.Score))+facet_wrap(~WRIA.Number, ncol=4)+labs(y="B-IBI Score", x="Year")+ theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("WatershedRKTtrend_PSSB.png"), plot = pp, width=20, height=10)




#####Site trends ####
library(EnvStats)

bibi <- ddply(bibi, 'Site.Code', mutate, por = length(unique(Year))) # add on period of record col
unique(bibi$Site.Code)
bibi<-droplevels(bibi[bibi$por > 9,])# look only at sites with >9 years of data
bibi <- ddply(bibi, 'Site.Code', mutate, timerange = (max(Year)-min(Year))) # add on period of record col
bibi<-droplevels(bibi[bibi$timerange > 9,])# look only at streams with >9 years of data
unique(bibi$Site.Code)

ggplot(bibi, aes(x=Site.Code, y=Year))+geom_tile()

ggplot(bibi, aes(x=Year, y=Overall.Score))+geom_smooth(se=F)+theme(legend.position = "none")+
  geom_point()



library(ggpmisc)

meanscores<-ddply(bibi, .(Year), summarize, meanScore=mean(Overall.Score), nyears=length(unique(Site.Code)), sd=sd(Overall.Score), min=min(Overall.Score), max=max(Overall.Score), median=median(Overall.Score))
meanscore<-ggplot(meanscores, aes(x=Year, y=meanScore))+geom_smooth(method="lm",se=F)+theme(legend.position = "none")+
  geom_point()+ylim(c(0,100))+labs(y="BIBI Score")+
  stat_poly_eq(formula = y ~ x, rr.digits = 2, coef.digits = 3, size = 4, color = "black", aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), parse = TRUE) 



scores<-ggplot(bibi )+
  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  geom_boxplot(aes(y=Overall.Score, x=as.numeric(paste0(Year)), group=as.factor(Year)), fill="white")+
  geom_point(data=meanscores, mapping=aes(x=as.numeric(paste0(Year)), y=meanScore))+
  geom_smooth(data=meanscores, mapping=aes(x=as.numeric(paste0(Year)), y=meanScore), method="lm",se=F)+
  labs(y="Overall Score", x="Year")+ theme(text = element_text(size = 28))+
  geom_quantile(data=bibi, mapping=aes(x=as.numeric(paste0(Year)), y=Overall.Score), quantiles=.5, color="red", size=1)#+
# geom_abline(slope=0.48, intercept = -914.2032, color="red", lwd=1)##adding Sen-Theil line and intercept calculated above



scal<-"Site.Code" # "Site.Code" or "Subbasin"

bibi<-subset(bibi, Site.Code!="")
form<-names(bibi)[5:(ncol(bibi)-2)]
# bibi.lr[bibi.lr$Year<2012,"Overall.Score"]<-(bibi.lr[bibi.lr$Year<2012,"Overall.Score"]*.972)+5.01
# for (i in 9:(ncol(bibi.lr))) {
for (var in form) {
  bibi<- ddply(bibi, .variables =scal, mutate, slope =lm(get(paste0(var)) ~ as.numeric(Year))[[1]][2])
  bibi<-ddply(bibi, .variables =scal, mutate, pval = summary(lm(get(paste0(var)) ~ as.numeric(Year)))[[4]][2,4])
  bibi$trend <- ifelse(bibi$pval < 0.05 & bibi$slope < 0, "negative", ifelse(bibi$pval < 0.05 & bibi$slope > 0, "positive", "none")) # add trend col
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.tau=kendallTrendTest( get(paste0(var))~ as.numeric(Year))$estimate[1])
  # bibi<-ddply(bibi,.(get(scal)),mutate, mk.tau=kendallTrendTest(as.formula(paste0(colnames(bibi[i]), " ~ ", "as.numeric(Year)")))$estimate[1])
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.Sen=kendallTrendTest( get(paste0(var))~ as.numeric(Year))$estimate[2])
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.pval=kendallTrendTest( get(paste0(var))~ as.numeric(Year))$p.value)
  bibi<-ddply(bibi,.(get(scal)),mutate, mk.inter=kendallTrendTest( get(paste0(var))~ as.numeric(Year))$estimate[3])
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

sampletimes<-ggplot(bibi, aes(x=Site.Code, y=Year, fill=MKSlope_Overall.Score),color="black")+geom_tile(color="black")+scale_fill_gradient2(name="Sen-Theil Slope")+
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
write.csv(overall_trends,"Site_Trends_coarse_PSSB.csv")

bibi_sum_site<-ddply(bibi, "Site.Code", summarise, mean_Overall.Score=mean(Overall.Score), SD_overallscore=sd(Overall.Score), SE_overallscore=sd(Overall.Score)/sqrt(length(Overall.Score)), len=length(Overall.Score), MK_trend=unique(MKtrend_Overall.Score), Avg_MKSlope_Overall.Score=mean(MKSlope_Overall.Score), Avg_MKtau_Overall.Score=mean(MKtau_Overall.Score), Avg_MKpval_Overall.Score=mean(MKpval_Overall.Score), sites=length(unique(Site.Code)), min=min(Overall.Score), max=max(Overall.Score), med=median(Overall.Score))
write.csv(bibi_sum_site,"score_RKTtrend_summarized_by_site_PSSB.csv")

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
  scale_colour_manual(values=rep(brewer.pal(8,"Dark2"), times=20))
# pp<-pp+geom_text(mapping= aes(x=2012, y=6, label=paste("Slope = ", round(get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2])), 3))), colour="black")+theme(text=element_text(size=20), plot.title=element_text(hjust=.5, size=24, face="bold"))
# ggsave(paste0("WRIA9_BIBI~timeUp_", scal, "_", var,   "-2.png"), plot = pp, width=20, height=10)
# }
print(pp)
# }
ggsave(paste0("OverallScores_Sites_PSSB.png"), plot = pp, width=20, height=10)


ggplot(bibi, aes(x=Year, y=Overall.Score))+  geom_rect(xmin=2001, xmax=2022, ymin=0, ymax=20, fill="firebrick1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=20, ymax=40, fill="tan1", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=40, ymax=60, fill="lightgoldenrod", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=60, ymax=80, fill="lightgreen", color="black")+
  geom_rect(xmin=2001, xmax=2022, ymin=80, ymax=100, fill="steelblue1", color="black")+
  # geom_smooth(se=F, method="lm")+
  geom_abline(aes(intercept=get(paste0("MKinter_",str_split(var, "_", n=2)[[1]][2])),slope=get(paste0("MKSlope_",str_split(var, "_", n=2)[[1]][2]))), color="blue", size=1)+
  theme(legend.position = "none")+
  facet_wrap(~Site.Code)+geom_point( aes(shape = c(shapes))) + scale_shape_identity()+labs(y="Overall B-IBI Score")


ggsave("SiteScores_PSSB.png", width=20, height=10)
