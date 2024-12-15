############################
### Libearies
###
library(glmnet)
library(dplyr)
library(caret)
library(tidyverse)

#############################
### Reading the data
###
meta.tbl <- read_delim("metadata_AQUAeD96.txt")
MAG.coverage.tbl <- read_delim("coverage/MAG_coverages.tsv")
data.tbl <- read_delim("metadata_16S.txt") %>% 
  select("SampleID", "nEQR", "Location", "Station") %>% 
  mutate(Loc_Stat = paste(Location, Station, sep = "_")) %>% 
  select(-Location, -Station) %>% 
  mutate(SampleID = as.character(SampleID)) %>% drop_na(nEQR)

#####################################################
### Preparation of data and put into same table
###
MAG.cov.tbl <- MAG.coverage.tbl %>% column_to_rownames(var = "MAG_id")
MAG.cov.tbl_t <- as.data.frame(t(MAG.cov.tbl))
all(meta.tbl$SampleID == rownames(MAG.cov.tbl_t)) 

meta.tbl_selected <- meta.tbl %>% select(SampleID, nEQR)
meta.tbl_selected$SampleID <- as.numeric(meta.tbl_selected$SampleID)

####################################################
### Merge the data
###
merged_data <- MAG.cov.tbl_t %>%
  rownames_to_column(var = "SampleID") 
merged_data$SampleID <- as.numeric(merged_data$SampleID)
merged_data <- inner_join(meta.tbl_selected, merged_data, by = "SampleID") %>% drop_na(nEQR)
data.tbl$SampleID <- as.numeric(data.tbl$SampleID)
data.tbl <- data.tbl %>% select(-nEQR)
joined_tbl <- full_join(merged_data, data.tbl, by = "SampleID") %>% drop_na(nEQR)
data_tabell_WGS_coverage <- joined_tbl %>%
  select(SampleID, Loc_Stat, nEQR, everything())
data_tabell_WGS_coverage$SampleID<- as.numeric(as.character(data_tabell_WGS_coverage$SampleID))
data_tabell_WGS_coverage <- data_tabell_WGS_coverage[order(data_tabell_WGS_coverage$SampleID), ]

#########################
### Anonymazation
###
location_col <- data_tabell_WGS_coverage[[2]]
unique_locations <- unique(location_col)
anonymized_locations <- setNames(seq_along(unique_locations), unique_locations)
data_tabell_WGS_coverage[[2]] <- anonymized_locations[location_col]

###########################
## save and load data
##
save(data_tabell_WGS_coverage, file = "data_tabell_WGS_coverage.RData")
load("data_tabell_WGS_coverage.RData")

