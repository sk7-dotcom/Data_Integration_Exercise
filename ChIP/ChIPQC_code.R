ChIPQC
library(BiocManager)
install("ChIPQC")
library(ChIPQC)

samples <- read.csv('ChIPQC/ChIPQC.csv')

chipObj <- ChIPQC(samples, annotation="hg19") 

ChIPQCreport(chipObj, reportName="ChIP QC report: Nanog and Pou5f1", reportFolder="ChIPQCreport")