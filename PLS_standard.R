##################################
### Libraries
###
library(tidyverse)
library(pls)
library(midiv)
library(compositions)

###########################
### Loading data
### 
# load("data_tabell_OTU95.RData")
# load("data_tabell_OTU97.RData")
 load("data_tabell_OTU99.RData")
# load("data_tabell_WGS_coverage.RData")
# load("data_tabell_WGS_pathways.RData")

#############################################################
### Choosing a standard data table for running LASSO
### for OTU data sets, removing sample 1087 that is missing on WGS data
###

data.tbl <- data_tabell_OTU99 %>% 
  filter(SampleID != "1087")
save_file_path <- "data_tabell_OTU99.RData"

#########################
### The response vector
###
y <- data.tbl$nEQR
names(y) <- data.tbl$SampleID

#############################################
### The cross-validation segments (subsets)
###
ulocstat <- unique(data.tbl$Loc_Stat)
segment.lst <- c()
for(i in 1:length(ulocstat)){
  segment.lst <- c(segment.lst, list(which(data.tbl$Loc_Stat == ulocstat[i])))
}

#################################
### The predictors (features)
###

X <- data.tbl %>%
  select(-c(SampleID, Loc_Stat, nEQR)) %>%  
  as.matrix()

# Set row names and apply CLR transformation
rownames(X) <- data.tbl$SampleID
X.clr <- t(apply(X, 1, clr))


###########################
### Fitting the PLS model
###
pls.obj <- plsr(y ~ X.clr, ncomp = 20, validation = "CV", segments = segment.lst)

######################################################
### Optimal number of components
### Find for which component we have minimum PRESS
### PRESS = Predicted Residual Sum of Squares
###
opt.comp <- which(pls.obj$validation$PRESS == min(pls.obj$validation$PRESS))
cat("Found best model at", opt.comp, "components\n")

#########################################
### The predictions for the best model
###
pls.preds <- pls.obj$validation$pred[,1,opt.comp]

###############################################################
### Extract the number of components and PRESS values
###
num_components <- 1:length(pls.obj$validation$PRESS)
cv_values <- pls.obj$validation$PRESS

########################################################
### Plot CV value (PRESS) vs. Number of Components
###
plot(num_components, cv_values, type = "b", pch = 19, col = "blue",
     xlab = "Number of Components", ylab = "CV Error (PRESS)",
     main = "0.99 OTU")
abline(v = opt.comp, col = "red", lty = 2) 
text(opt.comp, cv_values[opt.comp] + 0.5, labels = paste("Optimal:", opt.comp), col = "red")


########################################
### Calculate MSE, R2 og MAE
### 
mse <- mean((y - pls.preds)^2)
mae <- mean(abs(y - pls.preds))
r_squared <- 1 - sum((y - pls.preds)^2) / sum((y - mean(y))^2)

################################################
### Print the results
###
cat("Mean Squared Error (MSE):", mse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("R-squared (R²):", r_squared, "\n")

############################################
### Save results
###
save(pls.obj, opt.comp, pls.preds, mse, mae, r_squared, file = save_file_path)




