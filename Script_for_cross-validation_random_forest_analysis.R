#!/usr/bin/env Rscript

# Will need input arguments for activity data, phenotype data, column label for phenotype being analyzed, 
# the output file prefix, and optional training/test sets sampling seed
args = commandArgs(trailingOnly=TRUE)

# Test for at least 4 input arguments: if not, return an error
if (length(args)<4) {
  stop("Need at least 4 arguments with file names for activity data, phenotype data, the column label for the phenotype being analyzed (in quotes), and the output file prefix (in quotes). A seed number for the random sampling for training/test sets is an optional 5th argument", call.=FALSE)
} else if (length(args)==4) {
  # default seed
  args[5] = 123
}

# Send standard output to a log file.
sink(paste(args[4],"_rfcv_std_out.log",sep = ""),type = "output",split = T)

# Install and load required packages
print("Installing missing packages")
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

print("Loading required packages")
library(snowfall)
library(doParallel)
library(foreach)
library(future)
library(randomForest)
library(matrixStats)

# Read in data
print("Reading in data")
vip=read.table(args[1],sep = "\t",header = T,row.names = 1)
phenos=read.table(args[2],sep="\t",header = T,row.names = 1)

# Test for binary vs continuous phenotype. If neither, stop analysis.
num_pheno=is.numeric(phenos[,args[3]])
bin_pheno=(length(levels(as.factor(phenos[,args[3]])))==2)
if(!num_pheno && !bin_pheno){
  print("Phenotype data must be either continuous numeric or binary (numeric or categorical)")
  q("n")
}

# Split data 70/30 for training/test sets, but make sure the phenotype distribution is comparable between sets.
print("Splitting samples into training and test sets")
suppressWarnings(RNGkind(sample.kind = "Rounding")) # This is now necessary after R v3.6.0 for consistent results with pre-3.6.0 scripts because the default sampler changed.
set.seed(args[5])
rnd_samples=sample(colnames(vip),size = ceiling(dim(vip)[2]*0.7)) # 70% of samples rounded up to nearest integer
train_vip=vip[,rnd_samples]
train_phenos=phenos[rnd_samples,args[3]]
test_vip=vip[,!(colnames(vip) %in% rnd_samples)]
test_phenos=phenos[!(rownames(phenos) %in% rnd_samples),args[3]]

# These plots can indicate concerning deviations in the distributions of phenotype between the training and test sets
print("Plot phenotype distributions for training and test sets. If concerning difference are observed, consider either re-running this analysis with a different sampling seed number (as a 5th argument), or performing outlier sample removal prior to re-running the analysis.")

if(bin_pheno){ # If phenotype is binary...
  png(paste(args[4],"_distributions_in_training_and_test_sets.png",sep=""),width = 800,height = 600)
  par(mfrow=c(1,1))
  barplot(rbind(proportions(table(train_phenos)),proportions(table(test_phenos))),
          col = c("blue","red"),beside = T,legend = c("Train","Test"),
          xlab=args[3],ylab="Proportion",main="Phenotype distributions")
  dev.off()
} else{ # If phenotype is continuous...
  png(paste(args[4],"_distributions_in_training_and_test_sets.png",sep=""),width = 800,height = 600)
  par(mfrow=c(1,1))
  plot(density(train_phenos),xlab=args[3],main="Phenotype distributions",col="blue")
  lines(density(test_phenos),col="red")
  legend(x="topright",legend=c("Training","Test"),lty=par("lty"),col=c("blue","red"),text.col=c("blue","red"))
  dev.off()
}

# Identify the regulators whose expression or activities are best associated with the phenotype
print("Identifying the regulators whose activities are most significantly associated with the phenotype")
if(bin_pheno){ # If phenotype is binary...
  a_pheno=list()
  for(i in 1:dim(train_vip)[1]){
    a_pheno[[i]]=summary(glm((as.numeric(as.factor(train_phenos))-1)~as.numeric(train_vip[i,]),family = binomial())) # Have to coerce the phenotype data in this way to ensure it is interpreted as a binary input
  }
} else{ # If phenotype is continuous...
  a_pheno=list()
  for(i in 1:dim(train_vip)[1]){
    a_pheno[[i]]=summary(lm(train_phenos~as.numeric(train_vip[i,])))
  }
}
names(a_pheno)=rownames(train_vip)
sum_a_pheno=as.data.frame(matrix(nrow = length(a_pheno),ncol = 2))
for(i in 1:length(a_pheno)){
  sum_a_pheno[i,1]=coef(a_pheno[[i]])[2,1]
  sum_a_pheno[i,2]=coef(a_pheno[[i]])[2,4]
}
colnames(sum_a_pheno)=c("Beta","P")
rownames(sum_a_pheno)=names(a_pheno)
sum_a_pheno=sum_a_pheno[order(sum_a_pheno$P),] # re-order regulators by significance of association w/ phenotype
sum_a_pheno$BonferroniP=sum_a_pheno$P*dim(sum_a_pheno)[1]
bonf_sig=sum(sum_a_pheno$BonferroniP<0.05)
print(paste(bonf_sig,"regulators associated with phenotype with Bonferroni-adjusted P<0.05",sep = " "))

# Select regulators most significantly associated with phenotype for cross-validation random forest with a max of 500
if(bonf_sig<500){
  cvrf_regs=rownames(sum_a_pheno)[1:bonf_sig]
} else {
  cvrf_regs=rownames(sum_a_pheno)[1:500]
}
train_vip_top=train_vip[cvrf_regs,]
test_vip_top=test_vip[cvrf_regs,]

# Let's Z-transform the expression values and activities prior to RF modeling, as I did originally
cl=makeCluster(6)
registerDoParallel(cl)
train_zvip=list()
train_zvip=foreach(i=1:dim(train_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_vip_top[i,])))
train_zvip=as.data.frame(t(structure(train_zvip, row.names = c(NA, -length(train_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_zvip)=colnames(train_vip_top)
rownames(train_zvip)=rownames(train_vip_top)

cl=makeCluster(6)
registerDoParallel(cl)
test_zvip=list()
test_zvip=foreach(i=1:dim(test_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_vip_top[i,])))
test_zvip=as.data.frame(t(structure(test_zvip, row.names = c(NA, -length(test_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_zvip)=colnames(test_vip_top)
rownames(test_zvip)=rownames(test_vip_top)

# This is a wrapper function for running rfcv() with different RNG seeds
rfcv_wrapper=function(seed, train_x, train_y){
  set.seed(seed)
  temp=rfcv(trainx=train_x, trainy=train_y, scale = F, step = -5, cv.fold = 5)
  return(temp)
}

# Run cross-validation random forest to identify the number of features required to plateau the error 
# using the above identified associated regulators
print("Running random forest cross-validation analyses across different feature counts with 12 different sampling seeds")
cv_rndFor=list()
cl=makeCluster(future::availableCores())
registerDoParallel(cl)
if(bin_pheno){ # If phenotype is binary...
  cv_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar%
    rfcv_wrapper(seed=i, train_x=t(train_zvip), train_y=as.factor(train_phenos))
} else{ # If phenotype is continuous...
  cv_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar%
    rfcv_wrapper(seed=i, train_x=t(train_zvip), train_y=train_phenos)
}
stopCluster(cl)

# Average across 12 different seeds for the predictors counts in the rfcv() analysis, and then plot it.
cv_rndFor_error=cv_rndFor[[1]]$error.cv
for(i in 2:12){
  cv_rndFor_error=cbind(cv_rndFor_error,cv_rndFor[[i]]$error.cv)
}
cv_avg=rowMeans(cv_rndFor_error)
cv_sd=rowSds(cv_rndFor_error)

png(paste(args[4],"_RF_cross-validation_error_across_feature_counts.png",sep=""),width = 800,height = 600)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_avg)),
     y=cv_avg,
     ylim=range(c(cv_avg-cv_sd,cv_avg+cv_sd)),
     pch=19,
     xlab="Feature count",
     ylab="Mean cross-validation error",
     main=paste("Average",args[3],"rfcv() results over 12 different sampling seeds",sep = " "))
arrows(x0=as.numeric(names(cv_avg)),
       y0=cv_avg-cv_sd,
       x1=as.numeric(names(cv_avg)),
       y1=cv_avg+cv_sd,
       length = 0.05,
       angle = 90,
       code = 3) # This draws the error bars
dev.off()

# Save workspace for easy loading in the final MR analysis
print("Saving workspace that will be used in the final master regulator (MR) analysis")
save.image(paste(args[4],"_rfcv_workspace.RData",sep = ""))

# RFCV analysis finished
print("Random forest cross-validation analysis is complete. Look at the plot of mean cross-validation error across feature counts to determine how many features are needed to minimize prediction error while keeping the feature count as low as possible. This feature count will be used in the following master regulator analysis to pick out that number of putative master regulators (MRs) for your phenotype of interest.")
























