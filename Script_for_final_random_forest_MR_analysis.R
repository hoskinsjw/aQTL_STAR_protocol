#!/usr/bin/env Rscript

# Will need input arguments for RFCV workspace, column label for phenotype being analyzed, the number of 
# MRs to infer and the output file prefix
new_args = commandArgs(trailingOnly=TRUE)

# Test for at least 4 input arguments: if not, return an error
if (length(new_args)<4) {
  stop("Need at least 4 arguments with the file name for the random forest cross-validation (RFCV) analysis workspace, the column label for the phenotype being analyzed (in quotes), the number of MRs being inferred (based on cross-validation random forest analysis), and the output file prefix (in quotes).", call.=FALSE)
}

# Send standard output to a log file.
sink(paste(new_args[4],"_final_rf_MR_std_out.log",sep = ""),type = "output",split = T)

# Install and load required packages
if (!requireNamespace("snowfall", quietly = TRUE))
  install.packages("snowfall")
if (!requireNamespace("doParallel", quietly = TRUE))
  install.packages("doParallel")
if (!requireNamespace("foreach", quietly = TRUE))
  install.packages("foreach")
if (!requireNamespace("future", quietly = TRUE))
  install.packages("future")
if (!requireNamespace("randomForest", quietly = TRUE))
  install.packages("randomForest")
if (!requireNamespace("matrixStats", quietly = TRUE))
  install.packages("matrixStats")
if (!requireNamespace("caret", quietly = TRUE))
  install.packages("caret")

library(snowfall)
library(doParallel)
library(foreach)
library(future)
library(randomForest)
library(matrixStats)
library(caret)

# Read in workspace
print("Reading in workspace from random forest cross-validation analysis")
load(new_args[1])

# Train a random forest with the top phenotype-associated regulators to pick the top regulators (by importance type 1) 
# as representative phenotypic MRs. The number of desired MRs should have been determined based on the results of the 
# cross-validation random forest analysis in the previous step, and is fed into this analysis as the 3rd argument.
print(paste("Identifying",new_args[3],"representative", new_args[2],"master regulators (MRs)"))
suppressWarnings(RNGkind(sample.kind = "Rounding")) # This is now necessary after R v3.6.0 for consistent results with pre-3.6.0 scripts because the default sampler changed.
set.seed(135)
if(bin_pheno){ # If phenotype is binary...
  init_rndFor=randomForest(x=t(train_zvip),y=as.factor(train_phenos),ntree = 1000,keep.forest = T,importance = T)
} else{ # If phenotype is continuous...
  init_rndFor=randomForest(x=t(train_zvip),y=train_phenos,ntree = 1000,keep.forest = T,importance = T)
}
init_import=importance(init_rndFor,type = 1,scale = F)
init_import=init_import[order(init_import,decreasing = T),,drop=F]
mrs=rownames(init_import)[1:new_args[3]]
train_mrs=train_zvip[mrs,]
test_mrs=test_zvip[mrs,]

# Train the final random forest model using only the representative MRs just identified above.
print("Training final random forest model with only the representative MRs")
set.seed(314)
if(bin_pheno){ # If phenotype is binary...
  final_rndFor=randomForest(x=t(train_mrs),y=as.factor(train_phenos),xtest = t(test_mrs),ytest = as.factor(test_phenos),
                            ntree = 1000,keep.forest = T,importance = T)
} else{ # If phenotype is continuous...
  final_rndFor=randomForest(x=t(train_mrs),y=train_phenos,xtest = t(test_mrs),ytest = test_phenos,
                            ntree = 1000,keep.forest = T,importance = T)
}

# Output MR importance plot from final random forest model
pdf(file=paste(new_args[4],new_args[2],"MRs_importance_plot.pdf",sep = "_"),width = 6,height = 10)
varImpPlot(final_rndFor,type=1,scale = F,n.var = new_args[3],cex = 0.5,main=paste("Final random forest model representative",new_args[2],"MRs"))
dev.off()
    
# Check linear models or confusion matrix statistics for phenotype prediction accuracy
print(paste("Testing linear models associating final MR random forest-predicted",new_args[2],"to actual",new_args[2],"in both the training and test sets"))
if(bin_pheno){ # If phenotype is binary...
  sink(paste(new_args[4],new_args[2],"actual_vs_predicted_statistics.txt",sep = "_"),append = F,split = T)
  cat(paste("Training set","\n",sep = ""))
  print(confusionMatrix(final_rndFor$predicted,as.factor(train_phenos)))
  cat(paste("\n","Test set","\n",sep = ""))
  print(confusionMatrix(final_rndFor$test$predicted,as.factor(test_phenos)))
  sink()
} else{ # If phenotype is continuous...
  train_lm=lm(train_phenos~final_rndFor$predicted)
  test_lm=lm(test_phenos~final_rndFor$test$predicted)
  sink(paste(new_args[4],new_args[2],"actual_vs_predicted_linear_models.txt",sep = "_"),append = F,split = T)
  cat(paste("Training set linear model for actual~predicted",new_args[2],"\n"))
  summary(train_lm)$coefficients
  cat(paste("\nadjusted r-squared =",summary(train_lm)$adj.r.squared),"\n")
  cat(paste("\n\nTest set linear model for actual~predicted",new_args[2],"\n"))
  summary(test_lm)$coefficients
  cat(paste("\nadjusted r-squared =",summary(test_lm)$adj.r.squared),"\n")
  sink()
}

# Output the representative phenotypic MR list and update the workspace in case the final 
# random forest model is needed for any further validation against independent data sets.
write.table(mrs, file = paste(new_args[4],"representative",new_args[2],"MRs_list.txt",sep = "_"),quote = F,col.names = F,row.names = F)
save.image(paste(new_args[4],"_final_MR_random_forest_workspace.RData",sep = ""))

print("Master regulator analysis is complete!")















