'
Regression analysis using glm.

updates:2025.03.16 by Yang Sui(ysui@broadinstute.org)
to fix the AUC, now the output of AUC_95%_CI is Low_AUC_95%_CI, AUC and High_AUC_95%_CI

Also this version can output the real P not only P<2e-16

Aim:
   1. Fix the model
   2. Get the AUC of the model
   3. Get the Coefficients
   4. Get the Incremental_Model_Rsq
   5. Get the Pearson_Correlation


Notes:
    1. Read data from stdin and output results to stdout.
    2. Rows with NA (missing value) will be ignored.
    3. A glm model was fit, supported "family": https://www.statmethods.net/advstats/glm.html

Usage:
    GlmReression.R -f regression_formula -m family_link_function [-n null_model_formula] [-a AUC] [-i beta_95%CI] [-r ModelFitRsquare] [-p PearsonCorrelation] [-t PhenotypeName]

Options:
    -f string       Regression formular of Full model.eg:CAD~PRS+age+sex+PC.
    -m string       Family regression type that glm supports: binomial,gaussian,poisson ect.
    -n string       Regression formular of Null model.eg:CAD~age+sex+PC
    -a string       Calculate model AUC, can use if predict binary variable outcome.
    -t string       Phenotype name, necessary if -a is provided.
    -i string       Calculate coefficients 95% CI.
    -r string       Output the model fit r square. Provide anything if "yes".
    -p string       Output the pearson correlation between raw phenotype and predicted phenotype. Provide anything if "yes".
    -h --help       Show this screen.
    --version       Show version.
' -> doc

## Auto-detect and install needed packages.
## Written by someone
options (warn = -1)
.libPaths(c("rpackage_score",.libPaths()))
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table,pROC,dplyr,rsq)

opts <- docopt(doc,version="V1_2023.11.02")

#options(digits=5)

outAUC <- opts$a
outCI <- opts$i
outR <- opts$r
outP <- opts$p
N <- opts$n
t <- opts$t

df <- na.omit(read.table(file("stdin"),header = T,check.names=F))

reg.formula <- opts$f

print("Regression formula:")
print(reg.formula)

glm.fit <- glm(as.formula(reg.formula),family=opts$m,data=df)
null.glm.fit <- glm(as.formula(N),family=opts$m,data=df)

# Output regression results  
#summary(glm.fit)
#can output real p value.
summary(glm.fit)$coefficients

# Calculate 95% CI
if(is.null(outCI) == F){
    x= confint(glm.fit)
    y = as.data.frame(x)
    y$CI='95CI'
    y
}

# Calculate AUC
if(is.null(outAUC) == F){
    pred_val <- predict(glm.fit, type='response')
    #roc_obj <- roc(df[[colnames(df)[1]]], pred_val)
    roc_obj <- roc(df[[t]], pred_val) ###get the phenotype by their name
    myauc <- auc(roc_obj)
    ci_auc <- ci.auc(roc_obj)
    cat(c('AUC: ', myauc, '\n'),sep='\t')
    cat('AUC_95%_CI:', ci_auc[1],'-',ci_auc[2],'-',ci_auc[3], '\n')
}

# Calculate model fit correlation
if(is.null(outR) == F){
    model_r2 <- rsq(glm.fit,adj=F)
    model_r2 <- formatC(model_r2, digits = 4, format = "g")
    cat(c('Full_Model_Rsq: ', model_r2, '\n'),sep='\t')
}

# Calculate R2 of R2(full_model) - R2(null_model)
if(is.null(N) == F){
    model_r2_null <- rsq(null.glm.fit,adj=F)
    model_r2_null <- formatC(model_r2_null, digits = 4, format = "g")
    cat(c('Null_Model_Rsq: ', model_r2_null, '\n'),sep='\t')
    diff <- as.numeric(model_r2)-as.numeric(model_r2_null)
    diff <- formatC(diff, digits = 6, format = "g")
    cat(c('Incremental_Model_Rsq: ', diff, '\n'),sep='\t')
}

# Calculate correlation between fit phenotype and raw phenotype.
obs <- glm.fit$y
pred <- glm.fit$fitted.values

if(is.null(outP) == F){
     pr <- cor.test(obs,pred,method="pearson")$estimate
     pr <- formatC(pr, digits = 4, format = "g")
     cat(c('Pearson_Correlation:', pr, '\n'),sep='\t')
  # 95% CI of pearson correlation:
     CI <- cor.test(obs,pred,method="pearson")$conf.int[1:2]
     cat(c('95%_CI_PearsonR2:', CI, '\n'),sep='\t')
}
