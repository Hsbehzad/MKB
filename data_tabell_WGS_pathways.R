###############################
## libraries
##
library(dplyr)
library(tidyverse)
library(glmnet)
library(tidyr)

#############################
### Reading the data
###
MAG.coverage.tbl <- read_delim("coverage/MAG_coverages.tsv")
product.tbl <- read_delim("annotations/distills/product.tsv")


##########################################################
### Merge the two tables
###
colnames(MAG.coverage.tbl)[1] <- "genome"
colnames(product.tbl)[1] <- "genome"
WGS_merged_table <- merge(MAG.coverage.tbl, product.tbl, by = "genome")

######################################################
## Generate functional profiles using dot product
## calculate the dot product for each sample using loop
##
sample_columns <- colnames(MAG.coverage.tbl)[-1]
pathway_columns <- colnames(product.tbl)[-1]
functional_profiles <- list()
for (SampleID in sample_columns) {
  profile <- data.frame(Pathway = pathway_columns, Total_Contribution = numeric(length(pathway_columns)))
  for (pathway in pathway_columns) {
    profile[profile$Pathway == pathway, "Total_Contribution"] <- sum(WGS_merged_table[[SampleID]] * WGS_merged_table[[pathway]])
  }
  functional_profiles[[SampleID]] <- profile
}
for (SampleID in sample_columns) {
  print(paste("Functional profile for sample", SampleID))
  print(functional_profiles[[SampleID]])
}


###############################################################
## Combine the list of data frames into a single data frame 
## with an additional column for SampleID
##
functional_profiles_df <- bind_rows(functional_profiles, .id = "SampleID")
functional_profiles_wide <- functional_profiles_df %>%
  pivot_wider(names_from = Pathway, values_from = Total_Contribution)

######################################################
## Merge tables
##
data.tbl <- read_delim("metadata_16S.txt") %>% 
  select("SampleID", "nEQR", "Location", "Station") %>% filter(SampleID != 1087) %>% 
  mutate(Loc_Stat = paste(Location, Station, sep = "_")) %>% 
  select(-Location, -Station) %>% 
  mutate(SampleID = as.character(SampleID))
joined_pathway_tbl <- full_join(functional_profiles_wide, data.tbl, by = "SampleID") %>% drop_na(nEQR)
data_tabell_WGS_pathways <- joined_pathway_tbl %>%
  select(SampleID, Loc_Stat, nEQR, everything())
data_tabell_WGS_pathways$SampleID<- as.numeric(as.character(data_tabell_WGS_pathways$SampleID))
data_tabell_WGS_pathways <- data_tabell_WGS_pathways[order(data_tabell_WGS_pathways$SampleID), ]

#########################
### Anonymazation
###
location_col <- data_tabell_WGS_pathways[[2]]
unique_locations <- unique(location_col)
anonymized_locations <- setNames(seq_along(unique_locations), unique_locations)
data_tabell_WGS_pathways[[2]] <- anonymized_locations[location_col]

###########################
## save and load data
##
save(data_tabell_WGS_pathways, file = "data_tabell_WGS_pathways.RData")
load("data_tabell_WGS_pathways.RData")

