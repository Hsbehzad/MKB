library(tidyverse)
library(dplyr)

meta.tbl <- read_delim("metadata_AQUAeD96.txt")
MAG.tbl <- read_delim("MAGs/MAG_table.tsv")

MAG.coverage.tbl <- read_delim("coverage/MAG_coverages.tsv")

product.tbl <- read_delim("annotations/distills/product.tsv")

annot.tbl <- read_delim(str_c("annotations/annotations.tsv"))
