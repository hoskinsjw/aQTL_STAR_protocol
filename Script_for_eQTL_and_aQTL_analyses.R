#!/usr/bin/env Rscript

# Will need input arguments for input SNP dosage file, SNP locations file, gene expression/activity file,
# gene locations file, covariables file, cis P threshold, trans P threshold, and output prefix (in quotes)
args = commandArgs(trailingOnly=TRUE)

# Test for at least 8 input arguments: if not, return an error
if (length(args)<8) {
  stop("Need at least 8 arguments with the file paths for the SNP dosages, SNP locations, gene expression/activity, gene locations, and covariables, followed by the cis-QTL P threshold, trans-QTL P threshold, and the output file prefix (in quotes).", call.=FALSE)
}

# Send standard output to a log file.
sink(paste(args[8],"_std_out.log",sep = ""),split = T)

# Install and load required packages
if (!requireNamespace("MatrixEQTL", quietly = TRUE))
  install.packages("MatrixEQTL")

library(MatrixEQTL)

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = args[1];
snps_location_file_name = args[2];

# Gene expression file name
expression_file_name = args[3];
gene_location_file_name = args[4];

# Covariates file name
# Set to character() for no covariates
covariates_file_name = args[5];

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = as.numeric(args[6]);
pvOutputThreshold_tra = as.numeric(args[7]);

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Normal quantile transformation of gene expression data

for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cis_eqtls<-me$cis$eqtls
cat('Detected local QTLs:', dim(cis_eqtls)[1],'\n');
trans_eqtls<-me$trans$eqtls
cat('Detected distant QTLs:', dim(trans_eqtls)[1],'\n');

## Plot the Q-Q plot of local and distant p-values
jpeg(paste(args[8],"_QQ_plot.jpg",sep = ""))
plot(me)
dev.off()

write.table(cis_eqtls,paste(args[8],"_cis_results_",args[6],"_P_threshold.txt",sep = ""),sep="\t",quote = FALSE,row.names=FALSE)
write.table(trans_eqtls,paste(args[8],"_trans_results_",args[7],"_P_threshold.txt",sep = ""),sep="\t",quote = FALSE,row.names=FALSE)

