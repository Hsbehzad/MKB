########################################
### Libraries
###
library(tidyverse)
library(microseq)
library(dplyr)
library(ggplot2)

#########################################
### same codes will run on:
###
# vsearch_centroids_95.fasta
# vsearch_centroids_97.fasta
# vsearch_centroids_99.fasta

###########################################
### The sequence table
###
sequence_95.tbl <- readFasta("vsearch_centroids_95.fasta") %>%   
  mutate(clustered_size = str_remove(Header, "OTU[0-9]+;size=")) %>%
  mutate(clustered_size = as.numeric(clustered_size)) %>%
  mutate(Header = str_remove(Header, ";size=[0-9]++")) %>%
  arrange(desc(clustered_size))
View(sequence_95.tbl)

##########################################

demulti_95.tbl <- readFasta("/mnt/SCRATCH/helenbeh/vsearch/all_derep_minsize_95.fasta") %>%  
  mutate(centroid_size = str_remove(Header, ".+;size=")) %>%
  mutate(centroid_size = as.numeric(centroid_size)) %>%
  arrange(desc(centroid_size))
sequence_95.tbl <- demulti_95.tbl %>%
  select(-Header) %>%
  right_join(sequence_95.tbl, by = "Sequence")
ggplot(sequence_95.tbl) +
  geom_text(aes(x = centroid_size, y = clustered_size, label = Header)) +
  labs(x = "Centroid sizes", y = "Clustered sizes")

sample_95.tbl <- read_delim("metadata_16S.txt")
readcount_95.tbl <- read_delim("vsearch_readcounts_95.txt", delim = "\t") %>%  
  rename(OTU = `#OTU ID`)
sum_95.tbl <- readcount_95.tbl %>%
  select(-OTU) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "SampleID", values_to = "vsearch_readpairs")

ordered_sum_95.tbl <- sum_95.tbl[order(sum_95.tbl$vsearch_readpairs), ]
sample_95.tbl_selected <- sample_95.tbl %>% select(SampleID, nEQR)
ordered_sum_95.tbl <- ordered_sum_95.tbl %>%
  mutate(SampleID = as.numeric(SampleID))
sample_95.tbl_selected <- sample_95.tbl_selected %>%
  mutate(SampleID = as.numeric(SampleID))
ordered_sum_95.tbl <- ordered_sum_95.tbl %>%
  left_join(sample_95.tbl_selected, by = "SampleID") %>%
  drop_na(nEQR)

#######################################
### scatter plot
###

ggplot(ordered_sum_95.tbl, aes(x = SampleID, y = vsearch_readpairs, color = nEQR)) +
  geom_point(size = 3) + # Adjust point size as needed
  scale_color_gradient(low = "red", high = "blue") + 
  labs(
    title = "0.95 OTU identity",
    x = "Sample ID",
    y = "number of readpairs",
    color = "nEQR"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##################################
### bar plot
###
ggplot(ordered_sum_95.tbl, aes(x = reorder(factor(SampleID), nEQR), y = vsearch_readpairs, fill = nEQR)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") + 
  labs(
    title = "0.95 OTU identity",
    x = "SampleID",
    y = "number of readpairs",
    fill = "nEQR"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

########################################
### result table
###
# 
# sample_95.tbl <- left_join(sample_95.tbl, sum_95.tbl, by = "SampleID") %>%
# write_delim(sample_95.tbl, delim = "\t", file = "metadata_16_95.txt")

