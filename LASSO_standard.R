############################
### Libraries
###

library(tidyverse)
library(glmnet)
library(midiv)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

################################
### Loading the data
### 
# load("data_tabell_OTU95.RData") 
 load("data_tabell_OTU97.RData")
# load("data_tabell_OTU99.RData")
# load("data_tabell_WGS_coverage.RData")
# load("data_tabell_WGS_pathways.RData")

#############################################################
### Choosing a standard data table for running LASSO
### for OTU data sets, removing sample 1087 that is missing on WGS data
###
data_tabell <- data_tabell_OTU97 %>% 
  filter(SampleID != "1087")
data_tabell <- data_tabell_OTU97
save_file_path <- "data_tabell_OTU97.RData"

 
#############################################################
### Preparation and CLR transform
### 
y <- data_tabell$nEQR
segments <- as.numeric(factor(data_tabell$Loc_Stat)) 
X <- data_tabell |> 
  select(-SampleID, -nEQR, -Loc_Stat) |>  
  as.matrix()
X.clr <- t(apply(X, 1, clr, n.pseudo = mean(X) / 1000)) 
# X.clr <- na.omit(X.clr)
rownames(X.clr) <- data_tabell$SampleID


####################################
### Run LASSO
###
cv.obj <- cv.glmnet(x = X.clr, y = y, foldid = segments, keep = T)


##########################################
### Calculate prediction error and create plot
###
plot(cv.obj)
title(main=paste0("0.97 OTU"), 
      outer=TRUE, line=-1)

y.pred <- cv.obj$fit.preval[, cv.obj$lambda == cv.obj$lambda.1se] 

pred.tbl <- data.frame(y_obs = y, y_pred = as.numeric(y.pred))
View(pred.tbl)


#############################################################
### Total antall variabler og antall variabler valgt av LASSO
###
n_features <- ncol(X.clr)
print(paste("Total number of features:", n_features))

selected_features <- coef(cv.obj, s = "lambda.1se") 
selected_indices <- which(selected_features != 0)   
n_selected_features <- length(selected_indices) - 1 
print(paste("Number of features selected by LASSO:", n_selected_features))



########################################
### Calculate MSE (RMSE) and R2 and MAE
### 
mse <- mean((pred.tbl$y_obs - pred.tbl$y_pred)^2) 
mae <- mean(abs(pred.tbl$y_obs - pred.tbl$y_pred)) 

# Coefficient of Determination (R²)
# R² = 1 - (RSS / TSS)
rss <- sum((pred.tbl$y_obs - pred.tbl$y_pred)^2)  
tss <- sum((pred.tbl$y_obs - mean(pred.tbl$y_obs))^2)  
r_squared <- 1 - (rss / tss)

####################################
### Which MAG are selected
###

#############################################################
### Extract Selected MAGs and Coefficients
###
selected_features <- coef(cv.obj, s = "lambda.1se") 
selected_indices <- which(selected_features != 0)   
selected_coefficients <- as.numeric(selected_features[selected_indices]) 
selected_MAGs <- rownames(selected_features)[selected_indices] 
selected_MAGs <- selected_MAGs[selected_MAGs != "(Intercept)"] 

######################################
### MAGs and theri coefficeints
###

MAG_table <- data.frame(
   MAG_id = selected_MAGs,
   Coefficient = selected_coefficients[-1] 
 )


MAG_table <- MAG_table |>
  left_join(MAG.tbl, by = "MAG_id") |>
  arrange(desc(Coefficient))

print(MAG_table)
View(MAG_table)

write.csv(MAG_table,"MAG_table.csv")

#################################################
### Pathways and their coefficients
###

pathways_table <- data.frame(
  pathways = selected_MAGs,
  Coefficient = selected_coefficients[-1] 
)


MAG_ids <- product.tbl$genome
selected_pathways_in_product.tbl <- intersect(selected_MAGs, colnames(product.tbl)) 
selected_pathways_data <- product.tbl[, selected_pathways_in_product.tbl, drop = FALSE]
pathway_presence <- selected_pathways_data != 0
MAGs_with_pathways <- lapply(selected_pathways_in_product.tbl, function(pathway) {
  MAG_ids[pathway_presence[, pathway]]
})
pathway_MAG_summary <- data.frame(
  pathways = selected_pathways_in_product.tbl,
  MAGs = sapply(MAGs_with_pathways, paste, collapse = ", ") 
)

pathway_MAG_summary <- left_join(pathways_table, pathway_MAG_summary, by = "pathways")
View(pathway_MAG_summary)
write.csv(pathway_MAG_summary, "MAGs_with_selected_pathways.csv", row.names = FALSE)



final_pathway <- pathway_MAG_summary
final_pathway$MAGs <- as.character(final_pathway$MAGs)
final_pathway$MAG_count <- sapply(strsplit(final_pathway$MAGs, ","), length)
write.csv(final_pathway, "final_pathway.csv", row.names = FALSE)



pathway_MAG_summary_long <- pathway_MAG_summary %>%
  mutate(RowID = row_number()) %>%  
  separate_rows(MAGs, sep = ",")  
pathway_MAG_summary_long <- pathway_MAG_summary_long %>%
  mutate(MAGs = str_trim(MAGs)) 

pathway_MAG_summary_joined <- pathway_MAG_summary_long %>%
  left_join(MAG.tbl, by = c("MAGs" = "MAG_id"))
write.csv(pathway_MAG_summary_joined, "pathway_MAG_summary_joined.csv", row.names = FALSE)


####################################################
### Split the taxonomy column into separate parts
### and remove prefixes
###

taxonomy_split <- strsplit(as.character(pathway_MAG_summary_joined$GTDB_taxonomy), ";")
taxonomy_df <- do.call(rbind, lapply(taxonomy_split, function(x) {
  length(x) <- 7
  x
}))

taxonomy_clean <- apply(taxonomy_df, 2, function(column) {
  gsub("^[a-z]__", "", column)
})

taxonomy_df <- as.data.frame(taxonomy_clean, stringsAsFactors = FALSE)
colnames(taxonomy_df) <- c("domain", "phylum", "class", "order", "family", "genus", "species")

pathway_MAG_summary_joined_final <- cbind(pathway_MAG_summary_joined, taxonomy_df)

head(pathway_MAG_summary_joined_final)
phylum_counts <- table(pathway_MAG_summary_joined_final$phylum)
phylum_counts_df <- as.data.frame(phylum_counts)
colnames(phylum_counts_df) <- c("Phylum", "Count")
print(phylum_counts_df)
write.csv(phylum_counts_df, "phylum_counts_df.csv", row.names = FALSE)

###############################################
### Group by pathway and phylum, count the MAGs
###

pathway_MAG_summary_joined_final <- pathway_MAG_summary_joined_final %>%
  mutate(phylum = ifelse(is.na(phylum) | phylum == "", "Unknown", phylum))

phylum_count <- pathway_MAG_summary_joined_final %>%
  group_by(pathways, phylum) %>%
  summarise(MAG_count = n(), .groups = 'drop')

phylum_summary_table <- phylum_count %>%
  pivot_wider(names_from = phylum, values_from = MAG_count, values_fill = 0)

View(phylum_summary_table)
write.csv(phylum_summary_table, "phylum_summary_table.csv", row.names = FALSE)

####################################################
### barplot for pathways
###

ggplot(data = phylum_count, aes(x = pathways, y = MAG_count, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") +  
  labs(
    title = "Phylum Distribution Across Pathways",
    x = "Pathways",
    y = "Count of MAGs",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    legend.position = "bottom",  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  scale_fill_brewer(palette = "Set3")  

ggplot(data = phylum_count, aes(x = pathways, y = MAG_count, fill = phylum)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  theme_minimal() +
  theme(
    legend.position = "bottom",  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


###############################################
### Create a mapping of pathways to numbers for pathways
###

pathway_labels <- c(
  "1" = "Sulfur metabolism: thiosulfate => sulfite",
  "2" = "SCFA and alcohol conversions: Butyrate, pt 2",
  "3" = "Photosynthesis: Photosystem I",
  "4" = "Methanogenesis and methanotrophy: acetate => methane, pt 3",
  "5" = "Methanogenesis and methanotrophy: acetate => methane, pt 2",
  "6" = "Complex V: F-type ATPase, eukaryotes"
)

phylum_count$pathway_num <- as.factor(as.numeric(as.factor(phylum_count$pathways)))
ggplot(data = phylum_count, aes(x = pathway_num, y = MAG_count, fill = phylum)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  theme_minimal() +
  theme(
    legend.position = "bottom",  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    axis.text.y = element_text(size = 14)  
  ) +
  labs(x = "Pathways", y = "MAG Count", fill = "Phylum")


################################################
### Print the results
###
cat("Mean Squared Error (MSE):", mse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("R-squared (R²):", r_squared, "\n")

############################################
### Save results
###
save(cv.obj, pred.tbl, n_features, n_selected_features, mse, mae, r_squared, file = save_file_path)

