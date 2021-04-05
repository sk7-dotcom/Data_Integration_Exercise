#running ChIPQC
library(BiocManager)
install("ChIPQC")
library(ChIPQC)

samples <- read.csv('ChIPQC/ChIPQC.csv')

chipObj <- ChIPQC(samples) 

ChIPQCreport(chipObj, reportName="ChIP QC report", reportFolder="~/Desktop/ChIPQC/")
