#!/usr/bin/env Rscript

# Will need input arguments for expression data, covariables data, ARACNe network, and output file prefix
args = commandArgs(trailingOnly=TRUE)

# Test for at least 4 input arguments: if not, return an error
if (length(args)<4) {
  stop("Need at least 4 arguments! The first 3 are file names for expression data, covariables and ARACNe network (in that order). Output file prefix (in quotes) should be the 4th argument.", call.=FALSE)
}

# Send standard output to a log file.
sink(paste(args[4],"_VIPER_and_expression_normalization_std_out.log",sep = ""),type = "output",split = T)

# Install and load required packages
print("Installing missing packages")
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.15")
  }
if (!requireNamespace("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")
if (!requireNamespace("aracne.networks", quietly = TRUE))
  BiocManager::install("aracne.networks")
if (!requireNamespace("viper", quietly = TRUE))
  BiocManager::install("viper")
if (!requireNamespace("snowfall", quietly = TRUE))
  install.packages("snowfall")
if (!requireNamespace("doParallel", quietly = TRUE))
  install.packages("doParallel")
if (!requireNamespace("foreach", quietly = TRUE))
  install.packages("foreach")
if (!requireNamespace("future", quietly = TRUE))
  install.packages("future")
if (!requireNamespace("RNOmni", quietly = TRUE))
  install.packages("RNOmni")

print("Loading required packages")
library(preprocessCore)
library(aracne.networks)
library(viper)
library(snowfall)
library(doParallel)
library(foreach)
library(future)
library(RNOmni)

# Change the sample.kind = "Rounding" with RNGkind to ensure reproducibility of results generated with 
# scripts prior to R v3.6.0 that use RNG seeds. The default since v3.6.0 is "Mersenne-Twister".
suppressWarnings(RNGkind(sample.kind = "Rounding"))

# My modified write.regulon() function that is much faster than the one in aracne.networks. 
# NOTE: this function requires the "aracne.networks" package to be loaded.
write.regulon3<-function(regulon,file="",sep="\t",header=TRUE,n=Inf,regulator=NULL,cpus=future::availableCores(),toScreen=F){
  fullTab<-data.frame()
  if(is.null(regulator)){
    cl=makeCluster(cpus)
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    tab=list()
    tab=foreach(tf=names(regulon)) %dopar%
      cbind(rep(tf,length(names(regulon[[tf]]$tfmode))),names(regulon[[tf]]$tfmode),regulon[[tf]]$tfmode,regulon[[tf]]$likelihood)
    fullTab=as.data.frame(do.call("rbind",tab)) # Very rapidly convert a list of data.frames to data.frame
  } else {
    tf<-regulator
    x<-regulon[[tf]]
    targets<-names(x$tfmode)
    moas<-x$tfmode
    likelihoods<-x$likelihood
    tab<-cbind(rep(tf,length(targets)),targets,moas,likelihoods)
    fullTab<-rbind(fullTab,tab)
  }
  colnames(fullTab)=c("Regulator","Target","MoA","likelihood")
  rownames(fullTab)=NULL
  
  if(file!=""){
    write.table(fullTab,file = file,sep = sep, col.names = header, row.names = F,quote = F)
  }
  
  if(toScreen){
    print(fullTab)
  }
  
  return(fullTab)
}

# Read in the expression values used to infer the ARACNe networks, and the covars file that has the samples 
# in the final order for the eQTL and aQTL analyses
print("Reading in data sets")
genes=read.table(args[1],sep = "\t",header = T, row.names = 1)
genes_set=ExpressionSet(assayData=as.matrix(genes))
covars=as.data.frame(t(read.table(args[2],sep="\t",header = T,row.names = 1)))

# Scale (i.e. Z-transform) genes in the full log2(TPM) dataset.
print("Z-scaling expression data")
cl=makeCluster(future::availableCores())
registerDoParallel(cl)
genes_z=list()
genes_z=foreach(i=1:dim(genes)[1]) %dopar%
  as.numeric(scale(as.numeric(genes[i,])))
genes_z=as.data.frame(t(structure(genes_z, row.names = c(NA, -length(genes_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(genes_z)=colnames(genes)
rownames(genes_z)=rownames(genes)
z_set=ExpressionSet(assayData=as.matrix(genes_z))

# Convert the ARACNe network to interactome that includes mode of action (MoA) based on both the network and expression data
print("Generating interactome from ARACNe network and expression data")
regulon=aracne2regulon(args[3],genes_set,format="3col") # NOTE: The aracne network file cannot have a header!!!
regTable=write.regulon3(regulon,file=paste(args[4],"_interactome.txt",sep=""))

reg_count=sum(!duplicated(regTable$Regulator))
print(paste(reg_count,"distinct regulators in the interactome",sep = " "))

# Run VIPER to infer activities of regulators using the inferred interactome. 
vip=viper(z_set,regulon)
vip_results=exprs(vip)  # exprs() extracts the results matrix from the ExpressionSet object that is generated

# Filter viper activities and expression data to the same samples in the same order as in covars
filt_exp=genes[,na.omit(match(rownames(covars),colnames(genes)))]
filt_vip=as.data.frame(vip_results[,na.omit(match(rownames(covars),colnames(vip_results)))])

# Also, filter filt_exp for regulators with inferred activities in filt_vip
filt_regs=genes[na.omit(match(rownames(vip_results),rownames(genes))),na.omit(match(rownames(covars),colnames(genes)))]

# Check to see how well correlated the activities are to their respective expression values.
print("Checking regulator expression-activity correlations")
cl=makeCluster(future::availableCores())
registerDoParallel(cl)
expActCor=list()
startTime=Sys.time()
expActCor=foreach(i=1:dim(filt_regs)[1]) %dopar%
  cor(t(filt_regs)[,i],t(filt_vip)[,i])
expActCor=unlist(expActCor) # Very rapidly convert a list to data.frame
stopCluster(cl)
Sys.time()-startTime
names(expActCor)=rownames(filt_regs)

png(paste(args[4],"_activity_expression_correlation.png",sep=""),width = 800,height = 600)
plot(density(expActCor),main="Activity-Expression Correlation",xlab="Pearson r")
dev.off()

# Quantile normalize all expressed genes for final sample set
print("Quantile normalizing and inverse normal transforming expression data for downstream eQTL analyses")
qn_exp=as.data.frame(normalize.quantiles(as.matrix(filt_exp)))
rownames(qn_exp)=rownames(filt_exp)
colnames(qn_exp)=colnames(filt_exp)

# Inverse normal transform the quantile normalized filt_exp.
cl=makeCluster(future::availableCores())
registerDoParallel(cl)
int_exp=list()
int_exp=foreach(i=1:dim(qn_exp)[1],.packages = 'RNOmni') %dopar%
  RankNorm(as.numeric(qn_exp[i,]))
int_exp=as.data.frame(t(structure(int_exp, row.names = c(NA, -length(int_exp[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(int_exp)=colnames(qn_exp)
rownames(int_exp)=rownames(qn_exp)

# Filter the int_exp for regulators with inferred activities.
filt_int_exp=int_exp[na.omit(match(rownames(vip_results),rownames(int_exp))),]

# Write the filtered, normalized and transformed data to file
print("Writing output data to file")
write.table(cbind("Gene_name"=rownames(filt_int_exp),filt_int_exp),
            paste(args[4],"_regulators_QNorm_INT_expression.txt",sep=""),
            sep="\t",quote = F,row.names = F)
write.table(cbind("Gene_name"=rownames(int_exp),int_exp),
            paste(args[4],"_expressed_genes_QNorm_INT_expression.txt",sep=""),
            sep="\t",quote = F,row.names = F)
write.table(cbind("Gene_name"=rownames(filt_vip),filt_vip),
            paste(args[4],"_regulators_activities.txt",sep=""),
            sep="\t",quote = F,row.names = F)
