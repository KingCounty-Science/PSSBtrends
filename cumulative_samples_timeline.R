#=== === === === === === === ===
# Script started by Rebekah Stiling April 2025
# This creates figures and visualizations about the history of PSSB for a 4/9/2025 PSSB 2.0 workshop.
# I reference this websolution for help: https://stackoverflow.com/questions/65190609/cumulative-stacked-area-plot-for-counts-in-ggplot-with-r
#=== === === === === === === ===

# load relevant packages
library(tidyverse) #for data minipulation and plotting
library(viridis) # for color palet
library(RColorBrewer)

#create empty folder for graphs and csvs
output_folder <- paste0("output",Sys.Date()) 
dir.create(output_folder)

#although we don't need all the taxa data, the visit data is included in this framework.
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

#about how many samples are there?
samps <-raw |> select(Sample.ID) |> unique() 

#make the 
raw_2 <- raw |> mutate(vis_date = as.Date(raw$Visit.Date, format = '%m/%d/%Y')) 

df_sa<- raw_2 |> #_sa = Sample_Agency
  filter(vis_date < '2024-01-01') |> 
  select(Sample.ID, Agency, Visit.Date, vis_date) |> 
  unique()

df_sa <- df_sa |> mutate(year = year(vis_date))

#out of curiosity, how many agencies are there?
agencies <- df_sa |> select(Agency) |> unique()

# I need to make a dataframe with organization, year, count of events:
long_df <-df_sa |> 
  group_by(Agency, year) |> 
  summarise(count = n()) |> 
  ungroup()

range(long_df$year)

# "Complete" the data, fill in count = 0
final_dat <- long_df %>%
  complete(Agency, year = 1994:2023, fill = list(count = 0)) %>%
  group_by(Agency) %>%
  arrange(year) %>%  # Arrange by year so adding works
  mutate(aggcount = cumsum(count)) %>% 
  ungroup()

# Plot results
p<-final_dat %>%
  ggplot(aes(x = year, y = aggcount, fill = Agency)) +
  geom_area() +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  scale_x_continuous(limits = c(1993,2024),
                     breaks = seq(1995, 2025, 5),
                     minor_breaks = seq(1995, 2025, 5)) +
  labs(y = "Aggregate count of samples") + 
  theme(text=element_text(size=16), # Adjust X axis label size
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16))

p

ggsave(plot = p,
       filename = paste0(output_folder,"/agency_timeline_plot6x12_font_16-20.png"),
       units = "in", width = 12, height = 6)

  
  
