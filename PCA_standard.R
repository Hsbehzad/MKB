#########################################
### Libraries
###
library(ggplot2)
library(ggfortify)
library(magrittr)
library(midiv)

###################################
### Loading the data table
### 
load("data_tabell_OTU95.RData")
load("data_tabell_OTU97.RData")
load("data_tabell_OTU99.RData")
load("data_tabell_WGS_coverage.RData")
load("data_tabell_WGS_pathways.RData")

#############################################################################
### Choosing a standard data table for running PCA
### for OTU data sets, removing sample 1087 that is missing on WGS data
### 
data_tabell <- data_tabell_OTU95 %>% 
   filter(SampleID != "1087")
 
data_tabell <- data_tabell_WGS_pathways
save_file_path <- "data_tabell_WGS_pathways.RData"


#######################################################
### Relative abundance and CLR-transformed values for OTUs
###
y <- data_tabell$nEQR
names(y) <- data_tabell$SampleID
X <- data_tabell %>%
  select(starts_with("OTU")) %>%
  as.matrix()
rownames(X) <- data_tabell$SampleID
X.clr <- t(apply(X, 1, clr))
PCA_table <- cbind(nEQR = data_tabell$nEQR, X.clr)

#######################################################
### Relative abundance and CLR-transformed values for WGS_coverage
###
y <- data_tabell$nEQR
names(y) <- data_tabell$SampleID
X <- data_tabell %>%
  select(starts_with("AQU")) %>%
  as.matrix()
rownames(X) <- data_tabell$SampleID
X.clr <- t(apply(X, 1, clr))
PCA_table <- cbind(nEQR = data_tabell$nEQR, X.clr)

#######################################################
### Relative abundance and CLR-transformed values for WGS_pathways
###
y <- data_tabell$nEQR
names(y) <- data_tabell$SampleID
X <- data_tabell %>%
  select(4:ncol(data_tabell)) %>%
  as.matrix()
rownames(X) <- data_tabell$SampleID
X.clr <- t(apply(X, 1, clr))
PCA_table <- cbind(nEQR = data_tabell$nEQR, X.clr)

#######################################################
### Running PCA and Score plot
###
PCA <- prcomp(X.clr, scale = TRUE)
evar <- (PCA$sdev^2)/sum(PCA$sdev^2)
summary(PCA)
scores <- as.data.frame(PCA$x)
scores$nEQR <- y
ggplot(scores, aes(x = PC1, y = PC2, color = nEQR)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "red", high = "blue") +
  xlab(paste0("PC1 (", format(round(evar[1] * 100, 2), nsmall = 2), "% variance)")) +
  ylab(paste0("PC2 (", format(round(evar[2] * 100, 2), nsmall = 2), "% variance)")) +
  ggtitle("WGS_pathways") +
  theme(plot.title = element_text(hjust = 0.5))
  theme_minimal() +
  coord_fixed()

########################################################
### Extract the eigenvalues from the PCA object
###
 eigenvalues <- PCA$sdev^2

############################################################ 
### Create a scree plot, Proportion of variance explained
###
screeplot_PCA <- plot(evar, type = "b",
                         xlab = "Principal Component",
                         ylab = "Eigenvalue", main = "OTU 0.97")


##############################################
### Proportion of variance explained for first 20 components
###
plot(eigenvalues[1:20]/sum(eigenvalues), type = "b",
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained", main = "WGS_pathways, first 20 components")


############################################
### Save results
###
save(PCA, evar, file = save_file_path)



