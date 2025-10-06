'
Get weighted score using cv.glmnet

Author: Yang Sui(ysui1@mgh.harvard.edu)

Notes:
    1. Read data from stdin and output results to stdout.
    2. The input format should be "IID,outcome,score1,score2......scoreall"
    2. Rows with NA (missing value) will be ignored.
    3. A glm model was fit, supported "family": https://www.statmethods.net/advstats/glm.html

Usage:
    WeightedScore.R -m family_link_function [-a alpha]

Options:
    -m string       Family regression type that glm supports: binomial,gaussian,poisson ect.
    -a string       alpha = 0.5 for Elastic Net, 0 for Ridge, 1 for Lasso.
    --version       Show version.
' -> doc

## Auto-detect and install needed packages.
## Written by someone
options (warn = -1)
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table,pROC,dplyr,rsq)
library(glmnet)
opts <- docopt(doc,version="V1_2024.07.17")

#options(digits=5)


#########input datasets
df <- na.omit(read.table(file("stdin"),header = T,check.names=F))

##################add one in each time
#########get the weighted score
# weight score by Elastic Net, Ridge, and Lasso
# Adjust the selection to include all given PGS columns

PGS_columns <- colnames(df)[c(-1,-2)]
all_in_one <- df 
all_in_one <- na.omit(all_in_one)

# Initialize a dataframe to store weighted scores
weighted_scores_df <- data.frame(IID = all_in_one[,1], outcome = all_in_one[,2])

# Loop through each combination of PGS scores
for (i in 2:length(PGS_columns)) {
  selected_columns <- PGS_columns[1:i]
  x <- as.matrix(all_in_one[, selected_columns]) # Standardize scores (remove scale() if scores are already standardized)
  y <- all_in_one[,2] # Outcome
  
  fit <- cv.glmnet(x, y, alpha = as.numeric(opts$a), family = opts$m)  # alpha = 0.5 for Elastic Net, 0 for Ridge, 1 for Lasso

  # Extract the coefficients
  betas <- coef(fit, s = "lambda.min")[-1] # Extracts the coefficients corresponding to the lambda that gives the minimum cross-validated error. [-1] removes the intercept term from the coefficients. This is necessary because the intercept term is not part of the weighted score calculation.
  # Print betas to standard error stream
  cat(c(betas,'\n'), file = stderr())
  # Compute the weighted score
  weighted_score <- rowSums(sweep(x, 2, betas, `*`)) * (i / sum(betas)) # This multiplies each column of the matrix x by the corresponding element in the betas vector
  
  # Add the weighted score to the dataframe
  score_column_name <- paste0("weighted_score_", i)
  weighted_scores_df[[score_column_name]] <- weighted_score
}


write.table(weighted_scores_df ,"",row.names=F,col.names=T,quote=F)



