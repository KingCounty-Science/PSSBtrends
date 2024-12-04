## Use this code only for King County samples that do not have sample replicates. For samples with replicates, refer to the 2022 report folder for City of Seattle contract to see how I calculate densities using replicates. BUT for doing taxa trends with sample replicates, MAKE SURE 'Surface.Area' in the PSSB data has been reported accurately and consistently! Does it reflect the surface area for that one sample, or the surface area of the combined replicates? I suspect different agencies (and different people within agencies) might enter this information differently.


setwd(here::here())
# setwd("./outputs")
# rawcounts<-read.csv("Collapsed_Coarse_Taxa.csv")


###Raw Data Processing####
library(plyr)
library(openxlsx)
library(lubridate)
library(stringr)
library(tidyverse)
library(dplyr)
getwd()


atts<-read.xlsx("./Inputs/2012_taxa_attributes.xlsx")
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

file.path="./Inputs/taxonomy_data/"
raw<-taxaBind(file.path)

length(unique(raw$Project))
length(unique(raw$Agency))
length(unique(raw$Site.Code))
length(unique(raw$WRIA.Number))
length(unique(raw$Stream.or.River))
length(unique(raw$Subbasin))


#####Taxa Mapping####
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

######Prepare the data####
OTU$Visit.Date<-format(as.Date(OTU$Visit.Date, "%m/%d/%Y"))
OTU$Year<-year(OTU$Visit.Date)

OTU<-subset(OTU, is.na(QC.Replicate.Of.Sample.Code)|QC.Replicate.Of.Sample.Code=="")##remove QC replicates
# OTU<-subset(OTU, Non.B.IBI=="False")##remove non-target organisms ###update 10-30-2024-- added rolling exclusion later on to better reflect PSSB process.
OTU[which(OTU$OTU=="DNI"),"OTU"]<-OTU[which(OTU$OTU=="DNI"),"Taxon"]###These are marked as "DNI" in BCG translation table, but they aren't on B-IBI exclusion list. Adding back in for now.
OTU[which(is.na(OTU$OTU)),"OTU"]<-OTU[which(is.na(OTU$OTU)),"Taxon"]###These are taxa without a translation in BCG translation table. Adding back in as-is for now.

OTU$Unique<-as.logical(OTU$Unique)

##keep sample.code in to keep counts separated by sample, for density conversions later.
OTU_collapsed<-ddply(OTU, .(Visit.ID, OTU, Sample.Code, WRIA.Number, Agency, Basin, Subbasin, Stream.or.River, Project, Visit.Date, Year, Latitude, Longitude, Lab.Name, Site.Code), summarize, Quantity_OTU = sum(Quantity), Unique_OTU=any(Unique))

OTU_collapsed$Visit.Date<-as.Date(OTU_collapsed$Visit.Date, "%Y-%m-%d")

######Add correct Taxonomic Hierarchy####
########Create lookup table for taxa hierarchy.
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

append[append$Taxon=="Parathyas","Subclass"]<-"Acari" ##the BCG attribute table is missing subclass for this taxa, need to let Sean know to fix in next version.

hier_combined<-rbind(PSSB_taxa, append) ###combine
any(duplicated(hier_combined$Taxon))

# ##Append correct hierarchy to taxa
OTU_collapsed<-merge(OTU_collapsed, hier_combined, by.x="OTU", by.y="Taxon", all.x=T)
any(is.na(OTU_collapsed$Phylum))
OTU_collapsed[which(is.na(OTU_collapsed$Phylum)),]
any(OTU_collapsed$Phylum=="")

#####Roll-up by broad rules ####
KC_taxa_coarse<-OTU_collapsed

KC_taxa_coarse$OTU_COARSE<-""
names(KC_taxa_coarse)

coarse_rules<-data.frame(taxa=c("Oligochaeta", "Acari", "Gastropoda","Dytiscidae", "Simuliidae", "Chironomidae", "Trichoptera"), ranktouse=c("Subclass", "Subclass", "Family", "Family", "Family", "Family", "Genus"), rank=c("Subclass", "Subclass", "Class", "Family", "Family", "Family", "Order"))

# fine_rules<-data.frame(taxa=c("Oligochaeta", "Acari", "Gastropoda","Dytiscidae", "Simuliidae", "Chironomidae", "Trichoptera"), ranktouse=c("Genus", "Genus", "Genus", "Genus", "Genus", "Species", "Species"), rank=c("Subclass", "Subclass", "Class", "Family", "Family", "Family", "Order"))
# coarse_rules<-fine_rules

for (j in 1:nrow(coarse_rules)){
  
  STE_rank<-coarse_rules[j, "ranktouse"]
  rank<-coarse_rules[j, "rank"]
  taxa<-coarse_rules[j, "taxa"]
  index<-which(names(KC_taxa_coarse)== STE_rank)
  rankindex<-which(names(KC_taxa_coarse)== rank)
  halt1<-which(names(KC_taxa_coarse)== "Phylum")
  halt2<-which(names(KC_taxa_coarse)== "Subspecies")##changed from Species to Subspecies
  
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
##For fine STE, need to adjust species names in OTU column so that they include genus
KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Species),"OTU_COARSE"]<-paste0(KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Species),"Genus"]," ", KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Species),"Species"])
##For fine STE, need to adjust subgenus names in OTU column so that they include genus
KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Subgenus),"OTU_COARSE"]<-paste0(KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Subgenus),"Genus"]," (", KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Subgenus),"Subgenus"], ")")

library(plyr)
OTU_collapsed2<-ddply(KC_taxa_coarse, .(Visit.ID, Sample.Code, OTU_COARSE,Agency, WRIA.Number, Basin, 
                                        Subbasin, Stream.or.River, Project, Visit.Date, 
                                        Year, Latitude, Longitude, Lab.Name, Site.Code
), summarize, Quantity_OTU = sum(Quantity_OTU), Unique_OTU=any(Unique_OTU))

######Append on correct hierarchy after rollup####
###some taxa have multiple unique hierarchies because different entries got rolled up to the same level, and taxonomists have been uneven in entering in 'infraorder', 'suborder' and 'infraclass'. Need to run some loops to clean this up by selecting the most complete hierarchy available
names(KC_taxa_coarse)
new_hierarchy<-unique(KC_taxa_coarse[,c(19:41)])
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

check<-check[,-which(names(check) %in% "sum")]
check<-unique(check)

new_hierarchy<-subset(new_hierarchy, ! OTU_COARSE %in% check$OTU_COARSE)
new_hierarchy<-rbind(new_hierarchy, check)
countlength<-ddply(new_hierarchy, .(OTU_COARSE), summarize, count=length(OTU_COARSE))
fixthese<-countlength[countlength$count>1,]

OTU_collapsed3<-merge(OTU_collapsed2, new_hierarchy, by.x="OTU_COARSE", by.y="OTU_COARSE", all.x=T)
any(is.na(OTU_collapsed3$Phylum))
unique(OTU_collapsed3[which(is.na(OTU_collapsed3$Phylum)),]$OTU_COARSE)
any(OTU_collapsed3$Phylum=="")

#####Add Taxa Attributes####
################read in PSSB attribute table, do rolling lookup between coarse taxa hierarchy and attribute table 
names(OTU_collapsed3)
missing_atts<-unique(subset(OTU_collapsed3, select=c(18:39,1)))
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
attribs<-unique(attribs2)
any(duplicated(attribs$Taxon))
attribs[(which(duplicated(attribs$Taxon))),]
missing<-unique(KC_taxa_coarse$OTU_COARSE)[!unique(KC_taxa_coarse$OTU_COARSE) %in% attribs$Taxon] ##These taxa do not have a match in the PSSB attribute table

##append attributes from the lookup table we made
OTU_collapsed3<-left_join(OTU_collapsed3, attribs, by=c("OTU_COARSE"="Taxon"))
any(OTU_collapsed3$Clinger=="")
any(is.na(OTU_collapsed3$Clinger))
##Fill in attributes as "FALSE" for the taxa with no attribute matches
OTU_collapsed3[which(is.na(OTU_collapsed3$Clinger)),"Clinger"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Intolerant)),"Intolerant"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Long.Lived)),"Long.Lived"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Predator)),"Predator"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Tolerant)),"Tolerant"]<-FALSE

#####Taxa Exclusion####
###make sure to exclude taxa before scoring but after mapping

exlude<-read.xlsx("./Inputs/PSSB_exclusions.xlsx")
exlude<-subset(exlude, select=c(Taxon.Name, Excluded))


##Perform a series of rolling lookups to find a match in taxa hierarchies in the Taxon Name column of the exclusion table, starting with direct match, then from lowest hierarchy to highest

exlude2<-data.frame(Taxon.Name=character(), Excluded=character(), OTU_COARSE=character())
names(OTU_collapsed3)
test<-unique(subset(OTU_collapsed3, select=c(18:39,1)))
for (i in 1:ncol(test)){
  k<-(ncol(test)+1)-i ##work backwards through the hierarchy columns
  attribs<-merge(exlude, test[, c(k, 23)], by.x="Taxon.Name", by.y=names(test[k]))
  # names(attribs)<-str_replace(names(attribs), "2012.", "")
  names(attribs)<-str_replace(names(attribs), "OTU_COARSE.1", "OTU_COARSE")
  exlude2<-rbind(attribs, exlude2)
  test<-subset(test, !OTU_COARSE %in% exlude2$OTU_COARSE)
}

exlude2<-subset(exlude2, Excluded==T)
exlude2<-unique(exlude2)
OTU_collapsed3<-merge(OTU_collapsed3, exlude2, by="OTU_COARSE", all.x=T)
OTU_collapsed3<-subset(OTU_collapsed3, is.na(Excluded)|Excluded==F)
OTU_collapsed3<-OTU_collapsed3[,c(1:44)]

write.csv(OTU_collapsed3, "Collapsed_Coarse_Taxa_taxa_trends.csv")
##########Convert to Density####
# missinggrids<-raw[raw$Sample.Code %in% missinggridinfo$Sample.Code,c("Sample.Code","Surface.Area", "Sampling.Grid.Squares.Counted", "Total.Sampling.Grid.Squares")]
# missinggrids[,c("Sample.Code","Surface.Area", "Sampling.Grid.Squares.Counted", "Total.Sampling.Grid.Squares")]<-missinggridinfo[match(missinggrids$Sample.Code, missinggridinfo$Sample.Code),c("Sample.Code","Surface.Area", "Sampling.Grid.Squares.Counted", "Total.Sampling.Grid.Squares")]
setwd(here::here())
setwd("./inputs")
# merged_OTU<-read.csv("Merged_data_OTUs.csv", header=T)
# raw<-read.csv("taxonomy_data.csv") ##need the raw data just to get sampling grid and sampling area info


missinggridinfo<-read.csv("missing grid info.csv")###missing data obtained from Sean. One ABR sample with count <500, assumed fully counted
setwd(here::here())

rawcounts<-OTU_collapsed3
rawcounts<-subset(rawcounts, Agency=="King County - DNRP") ## Use this code only for King County samples that do not have sample replicates. For samples with replicates, refer to the 2022 report folder for City of Seattle contract to see how I calculate densities using replicates.

##add the missing grid data to the raw data

raw[raw$Sample.Code %in% missinggridinfo$Sample.Code,c("Sample.Code","Surface.Area", "Sampling.Grid.Squares.Counted", "Total.Sampling.Grid.Squares")]<-missinggridinfo[match(raw[raw$Sample.Code %in% missinggridinfo$Sample.Code,"Sample.Code"], missinggridinfo$Sample.Code),c("Sample.Code","Surface.Area", "Sampling.Grid.Squares.Counted", "Total.Sampling.Grid.Squares")]
raw[raw$Sample.Code %in% missinggridinfo$Sample.Code,c("Sample.Code","Surface.Area", "Sampling.Grid.Squares.Counted", "Total.Sampling.Grid.Squares")]

###calculate the proportion of the sampling grid counted
prop_counted<-subset(raw, select=c(Sample.Code, Sampling.Grid.Squares.Counted, Total.Sampling.Grid.Squares, Surface.Area))
aggreg<-aggregate(Quantity~Sample.Code+Sampling.Grid.Squares.Counted+Total.Sampling.Grid.Squares+Surface.Area+Surface.Area.Units+Visit.Date, raw, FUN=sum)
# write.csv(aggreg, "sample_density.csv")
prop_counted$prop_counted<-prop_counted$Sampling.Grid.Squares.Counted/prop_counted$Total.Sampling.Grid.Squares
prop_counted<-unique(prop_counted)
library(dplyr)
counts<-left_join(rawcounts, prop_counted, by="Sample.Code")

counts$total_sample_count<-counts$Quantity_OTU/counts$prop_counted
subset(counts, select=c(total_sample_count, Quantity_OTU, prop_counted, Surface.Area))
counts$density<-counts$total_sample_count/counts$Surface.Area
subset(counts, select=c(total_sample_count, Quantity_OTU, prop_counted, Surface.Area, density))
sample_id<-unique(subset(counts, select=Sample.Code))
sample_id$Sample.ID<-seq(1, length(sample_id$Sample.Code), 1)
counts<-left_join(counts, sample_id, by="Sample.Code")
unique(counts[which(is.na(counts$density)),]$Sample.Code)##these are samples that are missing sampling grid info or sampling area info
unique(counts[which(is.na(counts$density)),]$Agency)
unique(subset(counts, Agency=="King County - DNRP"&Year>2001&is.na(density)&Project=="Ambient Monitoring")$Sample.Code)
write.csv(counts, "Collapsed_Coarse_Taxa_density.csv")
counts<-read.csv("Collapsed_Coarse_Taxa_density.csv")
counts<-subset(counts, Agency=="King County - DNRP"&Year>2001&Project!="Regulatory Effectiveness")
counts<-subset(counts, select=-X)
##Append _UNIQUE to unique coarse OTUs
counts$OTU_COARSE_Unique<-ifelse(counts$Unique_OTU, paste0(counts$OTU_COARSE, "_UNIQUE"), counts$OTU_COARSE)


##################DPAC Preparation####

names(counts)
counts2<-counts[, c(3,2, 1,52,15, 50,17)] ##get just minimum sample data needed: Sample.Code, Visit.ID, Taxon Name, Taxon Name _Unique, Site.Code, density, and Unique 
counts2[which(duplicated(counts2)),] #make sure no data is duplicated

cast<- reshape2::dcast(counts2, OTU_COARSE+OTU_COARSE_Unique+Unique_OTU~Sample.Code, sum, value.var="density") ##reshape the data to wide format for DPAC

##add on hierarchical info

hier<-counts[, c(1, 18:39)]
hier<-unique(hier)
names(hier)
hier<-hier[,c(1:16,18, 17,19, 21, 20, 22,23)] ##re-order hierarchy to be in correct cascading order
hier[hier$Species!=""&!is.na(hier$Species),"Species"]<-hier[hier$Species!=""&!is.na(hier$Species),"OTU_COARSE"]## make the hierarchies for species and subgenus reflect OTU_COARSE for better matching in code loop below.
hier[hier$Subgenus!=""&!is.na(hier$Subgenus),"Subgenus"]<-hier[hier$Subgenus!=""&!is.na(hier$Subgenus),"OTU_COARSE"]
cast<-merge(cast, hier, by="OTU_COARSE") ##append hierarchies
names(cast)

cast<-cast[cast$OTU_COARSE_Unique!="NA_UNIQUE",]##These are Vashon sites with taxa that were not in the lookup table...need to figure out what to do with these guys
names(cast)

####Append _UNIQUE to hierarchies that match each unique coarse OTU, so that the code doesn't distribute uniquely identified taxa
### Previously I've done this step in excel using complex formula involving match function, but the loop below does the same thing. 
### BUT it will end in error if there's no match in the hierarchy for what is listed in the OTU_COARSE column.
### Errors should be investigated in OTU_collapse.R and the lookup table used in that code.
# cast[1,(ncol(cast)-21):ncol(cast)]##chec that these columns are the taxonomic hierarchy columns
for ( i in 1:nrow(cast)){
  if(cast[i,"Unique_OTU"]==T){
    # which(str_detect( cast$OTU_COARSE[i],t(cast[i,(ncol(cast)-21):ncol(cast)])))
    # str_detect(cast$OTU_COARSE[i],cast[i,(ncol(cast)-1)])
    idx<-match(cast$OTU_COARSE[i], cast[i,(ncol(cast)-21):ncol(cast)])## match OTU_COARSE to position in hierarchical columns
    cast[i,(ncol(cast)-21):ncol(cast)][,idx]<-paste0(cast[i,(ncol(cast)-21):ncol(cast)][,idx], "_UNIQUE") ##matching hierarchical position will get _UNIQUE added to it
    # print( cast[i,524:545][,idx])
  }
}

# cast[i,] ##if error in above loop, this will tell you where the error occurred

cast[cast==""]<-NA ##Any blank hierarchies convert to NA

comm=cast[, c(2,4:(ncol(cast)-22))] ##select columns OTU_COARSE_Unique, and all Visit.IDs
hier=cast[,c(2,(ncol(cast)-21):ncol(cast))] ##select columns OTU_COARSE_Unique and all hierarchical data
taxa.var="OTU_COARSE_Unique" ##identify column used for identification (OTU_COARSE_Unique)
out<-list(comm=comm, hier=hier, taxa.var=taxa.var)
class(out) <- "wide_class"

mod_dpac_s<-function (x, value.var.cols = NULL) 
{
  
  if (class(x) != "wide_class") 
    stop("Need an object of class 'wide_class'!")
  if (is.null(value.var.cols)) 
    stop("Must specify value.var!")
  
  comm <- x[["comm"]]
  hier <- x[["hier"]]
  taxa.var <- x[["taxa.var"]]
  temp<-data.frame("ID"=character(), taxa.var=character(), Quantity_old=numeric(), Quantity_new=numeric(), ambp=logical(), assigned=character())
  names(temp)[names(temp)=="taxa.var"]<-taxa.var
  for(value.var in names(comm[value.var.cols])){
    print(value.var)
    if (!value.var %in% names(comm)) 
      stop("value.var not found in data")
    if (any(is.na(comm[, value.var]))) 
      stop("No NAs in value.var allowed!")
    
    keep <- apply(hier, 2, function(x) any(is.na(x)))
    keep[rle(keep)$lengths[1]] <- TRUE
    keep[taxa.var] <- TRUE
    hier <- hier[, keep]
    lev <- rev(names(hier))
    lev <- lev[!lev %in% taxa.var]
    foo <- function(y, value.var) {
      # print(y[,which(names(y)==i)])
      if(all(is.na(y[,which(names(y)==i)]))){
        return(y)
      }
      else{
        if((which(names(y) == i) + 1)==(ncol(hier)-1)){
          childs <- !is.na(y[, which(names(y) == i) + 1])
        }
        else{
          childs <- apply(!is.na(y[, c((which(names(y) == i) + 1):(ncol(hier)-1))]), 1, any)
        }
        y[childs, taxa.var]
        parent <- !childs
        y[parent, taxa.var]
        if (sum(y[childs, value.var]) == 0 | all(childs)) {
          return(y)
        }
        else {
          w <- y[childs, value.var]/sum(y[childs, value.var])
          y[childs, "Quantity_new"] <- y[childs, value.var] + y[parent, value.var] * w
          y[parent, "Quantity_new"] <- 0
          y$ambp[parent] <- TRUE
          y$assigned[parent]<-paste(c(y[childs, taxa.var]), collapse=", ")
          return(y)
        }
      }}
    wdf <- cbind(hier, comm)
    wdf$ambp <- FALSE
    wdf$Quantity_new<-wdf[,value.var]
    wdf$ID<-value.var
    wdf$assigned<-""
    wdf<-subset(wdf, get(value.var)>0)
    for (i in lev[-1]) {
      wdf <- ddply(wdf, i, foo, value.var)
    }
    
    names(wdf)[which(names(wdf)==value.var)]<-"Quantity_old"
    commout <- wdf[, c("ID", taxa.var, "Quantity_old", "Quantity_new",  "ambp", "assigned")]
    
    temp<-rbind(temp, commout)
  }
  return(temp)
  
}

####DPAC#####
##This will take a several hours to run for a large dataset
dpac_out<-mod_dpac_s(x=out, value.var.cols=c(2:ncol(comm))) ##value.var.cols are the Visit.ID columns
dpac_out
counts2[!paste0(counts2$Sample.Code, counts2$OTU_COARSE_Unique) %in% paste0(dpac_out$ID, dpac_out$OTU_COARSE_Unique),]###these are missing taxa OTUs in regulatory effectiveness project. Don't worry about for now, not part of trends.

write.csv(dpac_out, "dpac_out.csv") #02012023 version, fixed Tipulidae s.l. taxonomy issue

dpac_out<-read.csv("dpac_out.csv")

#####################taxa tremds

library(rkt)
library(EnvStats)
library(tidyverse)
library(stringr)


# source('C:/Users/esosik/OneDrive - King County/Documents/biostats.R', encoding = 'UTF-8')
counts<-read.csv("Collapsed_Coarse_Taxa_density.csv")
lookup<-subset(counts, select=c(Site.Code, Sample.Code,Visit.ID, Project, Stream.or.River, Subbasin, WRIA.Number, Visit.Date, Year, Basin))
lookup<-unique(lookup)

# data<-read.csv("dpac_out.csv", header=T)
data<-dpac_out
data<-subset(data, Quantity_new!=0)### remove ambiguous taxa that have been re-assigned
data[which(duplicated(data[,c(1,2)])),]
data$OTU_COARSE_Unique<-str_replace(data$OTU_COARSE_Unique, "_UNIQUE", "") ##Remove "_UNIQUE"
data[which(duplicated(data[,c(1,2)])),]
excludevisits<-unique(subset(counts, Project=="Regulatory Effectiveness")$Sample.Code) ## remove reg effectiveness samples, since they use three replicates
data<-subset(data, !ID %in% excludevisits)
replicates<-unique(subset(raw, QC.Replicate.Of.Sample.Code!="")$Sample.Code)
excludereplicates<-unique(subset(counts, Sample.Code %in% replicates)$Sample.Code) ### remove all replicates and QC samples
data<-subset(data, !ID %in% excludereplicates)
data[which(duplicated(data[,c(1,2)])),]


keep<-ddply(rawcounts, .(Sample.Code, Visit.ID), summarize, sumorgs=sum(Quantity_OTU))
keep<-subset(keep, sumorgs>=450)
data<-subset(data, ID %in% keep$Sample.Code)


##clean up species data##
data<-reshape2::dcast(data, ID~OTU_COARSE_Unique, value.var = "Quantity_new")
data<-merge(data, lookup, by.x="ID", by.y="Sample.Code")
data$month <- format(as.Date(data$Visit.Date, '%Y-%m-%d'), "%m") # add month column
data$year <- format(as.Date(data$Visit.Date, '%Y-%m-%d'), "%Y") # add year column
data <- droplevels(data[data$month %in% c('07', '08', '09', '10'),]) # exclude samples collected outside of July - Oct
data <- droplevels(data[data$year > 2001,]) # subset to sites since 2002
# data<-subset(data, year<=2020) ##subset thru 2020
fixDUW<-unique(counts[counts$Sample.Code=="09DUW0277/05",]$Visit.ID)
data[data$ID==fixDUW,"Site.Code"] <-"09DUW0277" ##Correct SiteID
data[data$ID==fixDUW,"Stream.or.River"] <-"Riverton Creek (003D)" ##Correct SiteID
data[data$ID==fixDUW,"Subbasin"] <-"Duwamish River Subbasin" ##Correct SiteID
data[data$ID==fixDUW,"Basin"] <-"Duwamish - Green River Basin" ##Correct SiteID
data[data$ID==fixDUW,"WRIA.Number"] <-9 ##Correct SiteID
library(plyr)
data <- ddply(data, 'Site.Code', mutate, por = length(unique(year))) # add on period of record col
data<-droplevels(data[data$por > 9,])# look only at sites with >9 years of data


data[is.na(data)] <- 0
colnames(data)<-str_replace_all(colnames(data), " ", ".")
colnames(data)<-str_replace_all(colnames(data), "/", "_")
colnames(data)<-str_replace_all(colnames(data), "\\(", "_.")
colnames(data)<-str_replace_all(colnames(data), "\\)", "._")
colnames(data)<-str_replace_all(colnames(data), "#", "num")
data[which(duplicated(subset(data, select=c(Site.Code, Visit.Date, ID)))),]
data[duplicated(subset(data, select=c(Stream.or.River, ID, year))),]
data[duplicated(subset(data, select=c(Stream.or.River, ID, year)), fromLast=TRUE),]

names(data)
# data$Sample.Code <- paste(data$Site.Code, format(as.Date(data$Visit.Date, '%Y-%m-%d'), "%Y"), sep = "_") # add sample identifier
data$Sample.Code <- data$ID # add sample identifier





###2/9/2021 update fix some taxa that should be rolled up###
data$Leuctridae<-rowSums(data[, c("Moselia", "Despaxia.augusta", "Leuctridae")], na.rm=T)
data$Clinocerinae<-rowSums(data[, c("Wiedemannia", "Trichoclinocera", "Roederiodes", "Empididae.Genus.B", "Clinocera")], na.rm=T)
data$Hemerodromiinae<-rowSums(data[, c("Hemerodromiinae", "Hemerodromia", "Chelifera_Metachela", "Neoplasta")], na.rm=T)
data$Sphaeriidae<-rowSums(data[,c("Pisidium","Sphaerium","Sphaeriidae")], na.rm=T)
data$Glossiphoniidae<-rowSums(data[,c("Glossiphoniidae","Glossiphonia.elegans","Helobdella.stagnalis")], na.rm=T)
data$Hydrophilidae<-rowSums(data[,c("Hydrophilidae","Ametor","Anacaena", "Laccobius")], na.rm=T) ###,"Laccobius"
data$Gomphidae<-rowSums(data[,c("Gomphidae","Octogomphus.specularis")], na.rm=T)
data$Kogotus.Rickera.Cultus<-rowSums(data[,c("Kogotus_Rickera","Perlodidae","Cultus")], na.rm=T)
data$Hydroptila<-rowSums(data[,c("Hydroptilidae","Hydroptila")], na.rm=T)
data$Uenoidae..s.l.<-data$Uenoidae.s.l.##"Thremmatidae"
data$Trepaxonemata<-rowSums(data[,c("Trepaxonemata","Polycelis")], na.rm=T)
data$Optioservus<-rowSums(data[,c("Optioservus","Elmidae")], na.rm=T) ##Rotate through adding Elmidae to these taxa to see if it effects trends: Cleptelmis addenda, Heterlimnius.corpulentus, Lara, Narpus, Optioservus, Zaitzevia


data<-subset(data, select=-c(Sphaerium, Pisidium, Glossiphonia.elegans,Helobdella.stagnalis,
                             Ametor,Anacaena,Octogomphus.specularis,Kogotus_Rickera,Perlodidae,
                             Cultus,Hydroptilidae,Uenoidae.s.l.,Polycelis,Elmidae, Moselia, Laccobius, Despaxia.augusta, 
                             Wiedemannia, Trichoclinocera, Roederiodes, Empididae.Genus.B, Clinocera, Hemerodromia, Chelifera_Metachela, Neoplasta))#Laccobius, Thremmatidae

names(data)
names(data[c(1:238, 252:254, 239:251)])
data<-data[c(1:238, 252:254, 239:251)]### make sure all the taxa count data is before the metadata




write.csv(data, "KC_noAmbig_noReps_Rolledup_density.csv")
data<-read.csv("KC_noAmbig_noReps_Rolledup_density.csv")
##### end 9/21/2021 update####

################### Regional trends  ####
### convert to frequency data##

data_freq<-data
names(data_freq)
rownames(data_freq)<-data_freq$ID
data_freq<-data_freq[3:(ncol(data_freq)-13)]###select only taxa count data
data_freq[data_freq!=0]<-1 ##if present, mark as 1

###get rid of taxa that don't appear at least once in the dataset##
y1<-data_freq
min.cv=0
min.po=0
min.fo=1
max.po=100
max.fo=nrow(data_freq)
pct.missing=100
#statistical functions
cova<<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100
freq.occur<<-function(x,na.rm) sum(!x==0,na.rm=TRUE)
pct.occur<<-function(x,na.rm) sum(!x==0,na.rm=TRUE)/length(x)*100
fpct.missing<<-function(x,na.rm) sum(is.na(x))/length(x)*100

#delete offending variables
z<-as.matrix(apply(y1,2,cova)<min.cv |
               apply(y1,2,pct.occur)>max.po |
               apply(y1,2,pct.occur)<min.po |
               apply(y1,2,freq.occur)>max.fo |
               apply(y1,2,freq.occur)<min.fo |
               apply(y1,2,fpct.missing)>pct.missing)
z<-y1[,z[,1]==FALSE]


data_freq<-z

data_freq$ID<-rownames(data_freq)
data_freq<-merge(data[c(2, (ncol(data)-12):ncol(data))], data_freq, by="ID") ###add metadata back in

data_freq2<-data_freq %>% group_by(year) %>% summarise_if(is.numeric, max) #region-wide presence/absence, year to year (i.e. 1 if detected anywhere, 0 if not)
data_freq3<-data_freq %>% group_by(year) %>% summarise_if(is.numeric, sum) #frequency of occurence across sites, year to year (i.e. 23 if detected at 23 sites within a year)
names(data_freq3)
nsamp<-ddply(data, .(year), summarize, nsamp=length(unique(ID)))
data_freq3<-merge(data_freq3, nsamp, by.x="year", by.y="year")
data_freq3<-data_freq3 %>% mutate(across(Acari:Uenoidae..s.l., function(x) x/nsamp*100))
data_freq3<-subset(data_freq3, select=-nsamp)

# test<-reshape2::melt(data_freq3, id.vars=c("year", "Visit.ID", "WRIA.Number", "Year", "month", "por"))
# library(plyr)
# ties<-subset(ddply(test, .(variable, value), summarize, ties=length(value)), ties>1)
# FOQ<-read.csv("freq_of_occurrence_density_03132023_2002-2021.csv")
# FOQ_none<-subset(FOQ, mk.trend=="none")
# ties[ties$variable %in% FOQ_none$taxon,]
# ddply(test, .(variable), summarize, nyears=length(value))
# no_trends<-test[test$variable %in% FOQ_none$taxon,]
# ggplot(no_trends, aes(x=year, y=value))+geom_point()+facet_wrap(.~variable, scales="free")

#### region-wide pres/abs using logistic regression
### A positive trend indicates the taxa is occurring more frequently over time.

data.freq_trend<-data.frame(taxon=character(), slope=double(), inter=double(),pval=double(),trend=character(), samples=double(), sum.taxon=double(),stringsAsFactors = F)
names(data_freq2)
for (i in 7:ncol(data_freq2)) {
  data.lr<-NULL
  data.lr$taxon<-as.character(colnames(data_freq2[i]))
  data.lr$slope<-glm(as.formula(paste0(colnames(data_freq2[i]), " ~ ", "as.numeric(year)")), data=data_freq2, family="binomial")$coefficients[2]
  data.lr$inter<-glm(as.formula(paste0(colnames(data_freq2[i]), " ~ ", "as.numeric(year)")), data=data_freq2,family="binomial")$coefficients[1]
  data.lr$pval<-summary(glm(as.formula(paste0(colnames(data_freq2[i]), " ~ ", "as.numeric(year)")), data=data_freq2, family="binomial"))$coefficients[8]
  data.lr$trend <- as.character(ifelse(data.lr$pval < 0.05 & data.lr$slope < 0, "negative", ifelse(data.lr$pval < 0.05 & data.lr$slope > 0, "positive", "none"))) # add trend col
  data.lr$samples<-nrow(subset(data_freq2, data_freq2[i]!=0))
  data.lr$sum.taxon<-sum(subset(data_freq2, data_freq2[i]!=0)[i])
  
  data.lr<- as.data.frame(data.lr)
  data.freq_trend<- rbind(data.lr, data.freq_trend)
  data.freq_trend %>% mutate_if(is.factor, as.character) -> data.freq_trend
  
}
nrow(subset(data.freq_trend, trend!="none"))
subset(data.freq_trend, trend!="none")
subset(data.freq_trend, pval<.1)
unique(data.freq_trend$taxon)
write.csv(data.freq_trend, "region_pres-abs_density_03132023_2002-2021.csv")


####Look for trends in frequency of occurrence across streams, year to year. 
### A positive trend indicates the taxa is occurring at more streams over time.

data.freq_trend2<-data.frame(taxon=character(), mk.tau=double(), mk.Sen=double(),mk.pval=double(), mk.inter=double(), mk.trend=character(),samples=double(), sum.taxon=double(), stringsAsFactors = F)
names(data_freq3)
for (i in 7:ncol(data_freq3)) {
  data.lr<-NULL
  data.lr$taxon<-as.character(colnames(data_freq3[i]))
  data.lr$mk.tau<-kendallTrendTest(as.formula(paste0(colnames(data_freq3[i]), " ~ ", "as.numeric(year)")), data=data_freq3)$estimate[1]
  data.lr$mk.Sen<-kendallTrendTest(as.formula(paste0(colnames(data_freq3[i]), " ~ ", "as.numeric(year)")), data=data_freq3)$estimate[2]
  data.lr$mk.pval<-kendallTrendTest(as.formula(paste0(colnames(data_freq3[i]), " ~ ", "as.numeric(year)")), data=data_freq3)$p.value
  data.lr$mk.inter<-kendallTrendTest(as.formula(paste0(colnames(data_freq3[i]), " ~ ", "as.numeric(year)")), data=data_freq3)$estimate[3]
  data.lr$mk.trend <- as.character(ifelse(data.lr$mk.pval < 0.05 & data.lr$mk.tau < 0, "negative", ifelse(data.lr$mk.pval < 0.05 & data.lr$mk.tau > 0, "positive", "none")))
  data.lr$samples<-nrow(subset(data_freq3, data_freq3[i]!=0))
  data.lr$sum.taxon<-sum(subset(data_freq3, data_freq3[i]!=0)[i])
  data.freq_trend2<- rbind(data.lr, data.freq_trend2)
  data.freq_trend2 %>% mutate_if(is.factor, as.character) -> data.freq_trend2
  
}

# ggplot(data_freq3, aes(y=Cinygmula, x=year))+geom_point()
# ggplot(data_freq3, aes(y=Glutops, x=year))+geom_point()
# ggplot(data_freq3, aes(y=Psychodidae, x=year))+geom_point()
# ggplot(data_freq3, aes(y=Heterlimnius.corpulentus, x=year))+geom_point()
# ggplot(data_freq3, aes(y=Narpus, x=year))+geom_point()

nrow(subset(data.freq_trend2, mk.trend!="none"))
subset(data.freq_trend2, mk.trend!="none")
subset(data.freq_trend2, mk.pval<.1)
write.csv(data.freq_trend2, "freq_of_occurrence_density_10122023_2002-2021.csv")

test<-reshape2::melt(data_freq3, id.vars=c("year", "Visit.ID", "WRIA.Number", "Year", "month", "por"))
library(plyr)
ties<-subset(ddply(test, .(variable, value), summarize, ties=length(value)), ties>1)
FOQ<-read.csv("freq_of_occurrence_density_10122023_2002-2021.csv")
FOQ_none<-subset(FOQ, mk.trend=="none")
ties[ties$variable %in% FOQ_none$taxon,]
ddply(test, .(variable), summarize, nyears=length(value))
no_trends<-test[test$variable %in% FOQ_none$taxon,]
ggplot(no_trends, aes(x=year, y=value))+geom_point()+facet_wrap(.~variable, scales="free")


##Look for regional trends in taxa abundance, year to year.
### note: there's an extremely slight difference in rounding that causes this R project to have extremely slightly different abundance trend (region, subbasin and site) results for Skweltsa, Skwala, Parapsyche, Optioservus, Neophylax, Hemerodromiinae, Glossosoma, Dixa, Dicranota and Cinygma the **extremely** slight difference in densities (out to the 5th decimal place) is likely because of a different version of R being used.

names(data)
l=unique(c(as.character(data$Site.Code)))
# data <- data[complete.cases(data[, 1]), ]
str(data)
data$Site.Number<-as.numeric(factor(data$Site.Code, levels=l))
data$year<-as.numeric(data$year)

metrics<-colnames(data)[3:(ncol(data)-14)]
data.lr1<-data.frame(taxon=character(), rkt.pval=double(),rkt.Sen=double(), rkt.tau=double(), RKTtrend=character(),  samples=double(), sum.taxon=double(), stringsAsFactors = F)
# each<-metrics[248]
# each<-"Heterlimnius.corpulentus"
for (each in metrics) {
  data.lr<-NULL
  data.lr$taxon<-paste0(each)
  data.lr$rkt.pval <- with(data, rkt(date = year, y = get(each), block = Site.Number, correct = FALSE)[1])
  data.lr$rkt.Sen <- with(data, rkt(date = year, y = get(each), block = Site.Number, correct = FALSE)[3])
  data.lr$rkt.tau <- with(data, rkt(date = year, y = get(each), block = Site.Number, correct = FALSE)[12])
  data.lr$RKTtrend <- with(data.lr, ifelse(rkt.pval < 0.05 & rkt.tau < 0, "negative", ifelse(rkt.pval < 0.05 & rkt.tau > 0, "positive", "none")))
  data.lr$samples<-nrow(subset(data, get(each)!=0))
  data.lr$sum.taxon<-sum(subset(data, select=get(each)))
  data.lr<-unlist(data.lr)
  data.lr<-as.data.frame(t(data.lr))
  colnames(data.lr)[ncol(data.lr)-5]<-paste0("rkt.pval")
  colnames(data.lr)[ncol(data.lr)-4]<-paste0("rkt.Sen")
  colnames(data.lr)[ncol(data.lr)-3]<-paste0("rkt.tau")
  colnames(data.lr)[ncol(data.lr)-2]<-paste0("RKTtrend")
  data.lr1<- rbind(data.lr, data.lr1)
  print(each)
}

# ggplot(data, aes(y=Cinygmula, x=year, group=year))+geom_point()+geom_boxplot()
# ggplot(data, aes(y=Glutops, x=year, group=year))+geom_point()+geom_boxplot()
# ggplot(data, aes(y=Psychodidae, x=year, group=year))+geom_point()+geom_boxplot()
# ggplot(data, aes(y=Heterlimnius.corpulentus, x=year, group=year))+geom_point()+geom_boxplot()
# ggplot(data, aes(y=Narpus, x=year, group=year))+geom_point()+geom_boxplot()

nrow(subset(data.lr1, RKTtrend!="none"))
subset(data.lr1, RKTtrend!="none")
subset(data.lr1, rkt.pval<.1)
data.lr1$rkt.pval
str(data.lr)
write.csv(data.lr1, "regionwide_abundance_trends_RKT_density_03132023_2002-2021.csv")


#### Site/Basin trends ####

####look for trends in presence/absence by site using logistic regression
### A positive trend indicates the taxa is occurring more frequently over time within a given site.

scal<-"Site.Code" # "Site.Code" or "Subbasin"

data_freq_scal<-data_freq %>% group_by(get(scal), year) %>% summarise_if(is.numeric, max) #presence/absence, year to year (i.e. 1 if detected anywhere, 0 if not)
names(data_freq_scal)[1]<-paste0(scal)

data.lr3<-data.frame(taxon=character(), site=character(), slope=double(), inter=double(),pval=double(), trend=character(), samples=double(), sum.taxon=double(),stringsAsFactors = F)
names(data_freq_scal)
form<-names(data_freq_scal)[8:(ncol(data_freq_scal))]
# var<-"Zaitzevia"
for (var in form) {
  print(var)
  data.lr<-NULL
  data.lr$taxon<-ddply(data_freq_scal, .(get(scal)), summarize, var)[2]
  data.lr$site<-ddply(data_freq_scal, .(get(scal)), summarize, var)[1]
  data.lr$slope<-ddply(data_freq_scal, .(get(scal)), summarize, slope= glm(get(paste0(var)) ~ as.numeric(year), family="binomial")$coefficients[2])
  data.lr$inter<-ddply(data_freq_scal, .(get(scal)), summarize, inter= glm(get(paste0(var)) ~ as.numeric(year), family="binomial")$coefficients[1])
  data.lr$pval<-ddply(data_freq_scal, .(get(scal)), summarize, pval= summary(glm(get(paste0(var)) ~ as.numeric(year), family="binomial"))$coefficients[8])
  data.lr$trend <- as.character(ifelse(data.lr$pval$pval < 0.05 & data.lr$slope$slope < 0, "negative", ifelse(data.lr$pval$pval < 0.05 & data.lr$slope$slope > 0, "positive", "none"))) # add trend col
  data.lr$samples<-ddply(data_freq_scal, .(get(scal)), summarize, samples=length(which(get(var)!=0)))
  data.lr$sum.taxon<-ddply(data_freq_scal, .(get(scal)), summarize, sum.taxon=sum(get(var)))
  data.lr<- as.data.frame(data.lr)
  data.lr<-subset(data.lr, select=-c(slope.get.scal., inter.get.scal., pval.get.scal., samples.get.scal., sum.taxon.get.scal.))
  colnames(data.lr)<-c("taxon","site", "slope", "inter", "pval", "trend", "samples", "sum.taxon")
  data.lr3<- rbind(data.lr, data.lr3)
  data.lr3 %>% mutate_if(is.factor, as.character) -> data.lr3
  print(var)
}

nrow(subset(data.lr3, trend!="none"))
nrow(subset(data.lr3, pval<.1))
subset(data.lr3, pval<.1)
write.csv(data.lr3, paste0(scal,"_pres-abs_density_03132023_2002-2021.csv"))


##### Look for Abundance trends by site ### 
data.lr2<-data.frame(taxon=character(), site=character(),mk.tau=double(),mk.Sen=double(),mk.pval=double(), mk.inter=double(),  samples=double(), sum.taxon=double(),mk.trend=character(), stringsAsFactors = F)

# ddply(data,.(Site.Code),summarize, summ=(length(get(paste0(var)) )))

form<-names(data)[3:(ncol(data)-14)]
for (var in form){
  data.lr<-NULL
  data.lr$taxon<-ddply(data, .(Site.Code), summarize, var)[2]
  data.lr$site<-ddply(data,.(Site.Code),summarize, mk.Sen=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[2])[1]
  data.lr$mk.tau<-ddply(data,.(Site.Code),summarize, mk.tau=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[1])[2]
  data.lr$mk.Sen<-ddply(data,.(Site.Code),summarize, mk.Sen=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[2])[2]
  data.lr$mk.pval<-ddply(data,.(Site.Code),summarize, mk.pval=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$p.value)[2]
  data.lr$mk.inter<-ddply(data,.(Site.Code),summarize, mk.inter=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[3])[2]
  data.lr$mk.trend <- ifelse(data.lr$mk.pval < 0.05 & data.lr$mk.tau < 0, "negative", ifelse(data.lr$mk.pval < 0.05 & data.lr$mk.tau > 0, "positive", "none"))
  data.lr$samples<-ddply(data, .(Site.Code), summarize, samples=length(which(get(var)!=0)))
  data.lr$sum.taxon<-ddply(data, .(Site.Code), summarize, sum.taxon=sum(get(var)))
  data.lr<- as.data.frame(data.lr)
  data.lr<-subset(data.lr, select=-c(samples.Site.Code, sum.taxon.Site.Code))
  colnames(data.lr)<-c("taxon","site", "mk.tau", "mk.Sen", "mk.pval", "mk.inter", "trend","samples", "sum.taxon")
  data.lr2<- rbind(data.lr, data.lr2)
  
  data.lr2 %>% mutate_if(is.factor, as.character) -> data.lr2
  
}

subset(data.lr2, trend!="none")
subset(data.lr2, mk.pval<.1)
write.csv(data.lr2, "Site.Code_abundance_trends_density_03132023_2002-2021.csv")

##### Look for Spread trends by Basin ### 

data_freq_scal<-data_freq %>% group_by(Subbasin, year) %>% summarise_if(is.numeric, sum) #frequency of occurence across sites, year to year (i.e. 23 if detected at 23 sites within a year)

data.lr2<-data.frame(taxon=character(), site=character(),mk.tau=double(),mk.Sen=double(),mk.pval=double(), mk.inter=double(),  samples=double(), sum.taxon=double(),mk.trend=character(), stringsAsFactors = F)

form<-names(data_freq_scal)[8:ncol(data_freq_scal)]
for (var in form){
  data.lr<-NULL
  data.lr$taxon<-ddply(data_freq_scal, .(Subbasin), summarize, var)[2]
  data.lr$site<-ddply(data_freq_scal,.(Subbasin),summarize, mk.Sen=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[2])[1]
  data.lr$mk.tau<-ddply(data_freq_scal,.(Subbasin),summarize, mk.tau=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[1])[2]
  data.lr$mk.Sen<-ddply(data_freq_scal,.(Subbasin),summarize, mk.Sen=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[2])[2]
  data.lr$mk.pval<-ddply(data_freq_scal,.(Subbasin),summarize, mk.pval=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$p.value)[2]
  data.lr$mk.inter<-ddply(data_freq_scal,.(Subbasin),summarize, mk.inter=kendallTrendTest(get(paste0(var)) ~ as.numeric(year))$estimate[3])[2]
  data.lr$mk.trend <- ifelse(data.lr$mk.pval < 0.05 & data.lr$mk.tau < 0, "negative", ifelse(data.lr$mk.pval < 0.05 & data.lr$mk.tau > 0, "positive", "none"))
  data.lr$samples<-ddply(data_freq_scal, .(Subbasin), summarize, samples=length(which(get(var)!=0)))
  data.lr$sum.taxon<-ddply(data_freq_scal, .(Subbasin), summarize, sum.taxon=sum(get(var)))
  data.lr<- as.data.frame(data.lr)
  data.lr<-subset(data.lr, select=-c(samples.Subbasin, sum.taxon.Subbasin))
  colnames(data.lr)<-c("taxon","site", "mk.tau", "mk.Sen", "mk.pval", "mk.inter", "trend","samples", "sum.taxon")
  data.lr2<- rbind(data.lr, data.lr2)
  
  data.lr2 %>% mutate_if(is.factor, as.character) -> data.lr2
  
}


subset(data.lr2, trend!="none")
subset(data.lr2, mk.pval<.1)
write.csv(data.lr2,  "Subbasin_Spread_trends_density_03132023_2002-2021.csv")


#Look for RKT abundance trends by basin ####
names(data)
l=unique(c(as.character(data$Site.Code)))
# data <- data[complete.cases(data[, 1]), ]
str(data)
data$Site.Number<-as.numeric(factor(data$Site.Code, levels=l))
data$year<-as.numeric(data$year)

data.lr2<-data.frame(taxon=character(), site=character(), rkt.pval=double(),rkt.Sen=double(), rkt.tau=double(), RKTtrend=character(),  samples=double(), sum.taxon=double(), stringsAsFactors = F)

form<-colnames(data)[3:(ncol(data)-14)]
# each<-"Cleptelmis.addenda"
# unique(data$Basin  )
# data.frame(Basin=unique(data$Basin), Taxon=each)

detach(package:dplyr, unload=T) ##sometimes dplyr and plyr packages don't play nicely and we get errors with get() in the loop below. Toggle dply off and on to fix
library(dplyr)
paste0(each)
for (each in form){
  data.lr<-NULL
  data.lr_taxon<-NULL
  data.lr_rkt.pval<-NULL
  data.lr_rkt.Sen<-NULL
  data.lr_rkt.tau<-NULL
  data.lr_taxon<-ddply(data, .(Subbasin), summarize, each)[c(1,2)]
  data.lr_rkt.pval<-as.data.frame(ddply(data,.(Subbasin),summarize, rkt.pval=rkt(date = year, y = get(each), block=Site.Number, correct=FALSE)[1]))
  data.lr_rkt.pval$rkt.pval<-unlist(data.lr_rkt.pval[,2])
  data.lr_rkt.Sen<-as.data.frame(ddply(data,.(Subbasin),summarize, rkt.Sen=rkt(date = year, y = get(each), block=Site.Number, correct=FALSE)[3]))
  data.lr_rkt.Sen$rkt.Sen<-unlist(data.lr_rkt.Sen[,2])
  data.lr_rkt.tau<-as.data.frame(ddply(data,.(Subbasin),summarize, rkt.tau=rkt(date = year, y = get(each), block=Site.Number, correct=FALSE)[12]))
  data.lr_rkt.tau$rkt.tau<-unlist(data.lr_rkt.tau[,2])
  data.lr<-Reduce(function(x, y) merge(x, y), list(data.lr_taxon, data.lr_rkt.pval, data.lr_rkt.Sen, data.lr_rkt.tau))
  data.lr$RKTtrend <- ifelse(data.lr$rkt.pval < 0.05 & data.lr$rkt.tau < 0, "negative", ifelse(data.lr$rkt.pval < 0.05 & data.lr$rkt.tau > 0, "positive", "none"))
  data.lr$samples<-ddply(data, .(Subbasin), summarize, samples=length(which(get(each)!=0)))$samples
  data.lr$sum.taxon<-ddply(data, .(Subbasin), summarize, sum.taxon=sum(get(each)))$sum.taxon
  
  data.lr2<- rbind(data.lr, data.lr2)
  
  
  
}


nrow(subset(data.lr2, RKTtrend!="none"))
subset(data.lr2, rkt.pval<.05)
write.csv(data.lr2, "Subbasin_abundance_trends_RKT_density_03132023_2002-2021.csv")

