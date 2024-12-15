##################################
## libraries
##
library(dplyr)
library(tidyr)
library(tidyverse)


#######################
## read data
##
readcount_95.tbl <- read_delim("vsearch_readcounts_95.txt", delim = "\t") %>%  
  rename(OTU = `#OTU ID`)


##########################
## transpose the table
##
transposed_95_tbl <- readcount_95.tbl %>%
  pivot_longer(cols = -OTU, names_to = "SampleID", values_to = "value") %>%  
  pivot_wider(names_from = OTU, values_from = value)  
data.tbl_95 <- read_delim("metadata_16S.txt") %>% 
  select("SampleID", "nEQR", "Location", "Station") %>% 
  mutate(Loc_Stat = paste(Location, Station, sep = "_")) %>% 
  select(-Location, -Station) %>% 
  mutate(SampleID = as.character(SampleID))

#########################
## join tables
##
joined_tbl <- left_join(transposed_95_tbl, data.tbl_95, by = "SampleID")
data_tabell_OTU95 <- joined_tbl %>%
  select(SampleID, Loc_Stat, nEQR, everything()) %>% drop_na(nEQR)

data_tabell_OTU95$SampleID<- as.numeric(as.character(data_tabell_OTU95$SampleID))
data_tabell_OTU95 <- data_tabell_OTU95[order(data_tabell_OTU95$SampleID), ]

#######################################
### Anonymazation
###
location_col <- data_tabell_OTU95[[2]]
unique_locations <- unique(location_col)
anonymized_locations <- setNames(seq_along(unique_locations), unique_locations)
data_tabell_OTU95[[2]] <- anonymized_locations[location_col]

###########################
## save and load data
##
save(data_tabell_OTU95, file = "data_tabell_OTU95.RData")
load("data_tabell_OTU95.RData")
