##ChIP-RNA-Seq integration: Bio722 Workshop


#align to S.venezuelae genome with BowTie2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rbowtie2")
BiocManager::install("Rsamtools")
BiocManager::install("geneXtendeR")


library(Rbowtie2)
library(Rsamtools)
library(geneXtendeR)


bowtie2_build(references="ChIP/ref/vnz_chromosome.fna", bt2Index = "ChIP/index/venezuelae_bowtie", overwrite = TRUE)

bowtie2(bt2Index = "ChIP/index/venezuelae_bowtie", 
        samOutput = "ChIP/sam_out/Sam_Output_flagged_v1", 
        seq1 = "ChIP/raw_data/lsr2ChIP_flagged_v1.fastq", 
        seq2 = NULL, 
        overwrite = TRUE)

bowtie2(bt2Index = "ChIP/index/venezuelae_bowtie", 
        samOutput = "ChIP/sam_out/Sam_Output_flagged_v2", 
        seq1 = "ChIP/raw_data/lsr2ChIP_flagged_v2.fastq", 
        seq2 = NULL)

bowtie2(bt2Index = "ChIP/index/venezuelae_bowtie", 
        samOutput = "ChIP/sam_out/Sam_Output_unflagged_v1", 
        seq1 = "ChIP/raw_data/lsr2ChIP_unflagged_control1.fastq", 
        seq2 = NULL)

bowtie2(bt2Index = "ChIP/index/venezuelae_bowtie", 
        samOutput = "ChIP/sam_out/Sam_Output_unflagged_v2", 
        seq1 = "ChIP/raw_data/lsr2ChIP_unflagged_control2.fastq", 
        seq2 = NULL)

#convert SAM to BAM

asBam(file = "ChIP/sam_out/Sam_Output_flagged_v1", destination = "ChIP/bam_out/flagged_v1_asBAM", overwrite = TRUE)
asBam(file = "ChIP/sam_out/Sam_Output_flagged_v2", destination = "ChIP/bam_out/flagged_v2_asBAM", overwrite = TRUE)
asBam(file = "ChIP/sam_out/Sam_Output_unflagged_v1", destination = "ChIP/bam_out/unflagged_v1_asBAM", overwrite = TRUE)
asBam(file = "ChIP/sam_out/Sam_Output_unflagged_v2", destination = "ChIP/bam_out/unflagged_v2_asBAM", overwrite = TRUE)

#sort BAM alignment files by genomic coordinates (instead of by name)

sortBam("ChIP/bam_out/flagged_v1_asBAM.bam", "ChIP/sorted_bam/flagged_v1_asBAM")
sortBam("ChIP/bam_out/flagged_v2_asBAM.bam", "ChIP/sorted_bam/flagged_v2_asBAM")
sortBam("ChIP/bam_out/unflagged_v1_asBAM.bam", "ChIP/sorted_bam/unflagged_v1_asBAM")
sortBam("ChIP/bam_out/unflagged_v2_asBAM.bam", "ChIP/sorted_bam/unflagged_v2_asBAM")

#index BAM files

indexBam("ChIP/sorted_bam/flagged_v1_asBAM.bam", "ChIP/index_bam/flagged_v1_asBAM")
indexBam("ChIP/sorted_bam/flagged_v2_asBAM.bam", "ChIP/index_bam/flagged_v2_asBAM")
indexBam("ChIP/sorted_bam/unflagged_v1_asBAM.bam", "ChIP/index_bam/unflagged_v1_asBAM")
indexBam("ChIP/sorted_bam/unflagged_v2_asBAM.bam", "ChIP/index_bam/unflagged_v2_asBAM")

#Filter BAM files 

filterBam("ChIP/sorted_bam/flagged_v1_asBAM.bam", "ChIP/filtered_bam/flagged_v1_asBAM", index = "ChIP/index_bam/flagged_v1_asBAM.bam.bai")
filterBam("ChIP/sorted_bam/flagged_v2_asBAM.bam", "ChIP/filtered_bam/flagged_v2_asBAM", index = "ChIP/index_bam/flagged_v2_asBAM.bam.bai" )
filterBam("ChIP/sorted_bam/unflagged_v1_asBAM.bam", "ChIP/filtered_bam/unflagged_v1_asBAM", index = "ChIP/index_bam/unflagged_v1_asBAM.bam.bai" )
filterBam("ChIP/sorted_bam/unflagged_v2_asBAM.bam", "ChIP/filtered_bam/unflagged_v2_asBAM", index = "ChIP/index_bam/unflagged_v2_asBAM.bam.bai" )

#Peak Calling 

