Welcome to workshop portion of today’s tutorial!

This is part 3 of 3 sections towards integrating our data: Data
Integration

For more context and details on each step look under the slides folder
in the main repo.

We will do all the analysis in R.

------------------------------------------------------------------------

#### Package installation

To start this analysis, we will start by calling in all the packages we
will need. Install them using the code below if you have not already
done so. You also want to `setwd()` again to what it was in the
ChIP-analysis aspect of the workshop.

``` r
#Packages Used ----

#install.packages('tidyverse')
#install.packages('VennDiagram')
#install.packages('data.table')
#install.packages('dplyr')
#install.packages('tidyr')
#install.packages('ggplot2')

library(tidyverse)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(data.table)
library(ggplot2)

# setwd("~/Desktop/Data_Int_Files") # <-- edit this based on the directory of your choice
```

#### Data cleanup and final check

To cleanup the output from the `Bedtools` function to just include the
GeneID and the distance to the nearest gene the following code will be
used. It will result in a file called ChIP\_data.csv that will we
continue with.

``` r
#Cleaning up ChIP-seq output ----

peaks_closest_features <- read.delim("C:/Users/sreed/OneDrive - McMaster University/Documents/BIOLOGY722/Data_Integration_Exercise/ChIP/MACS3/peaks_closest_subset.txt", header = FALSE)
#the way the GFF file is constructed we have to clean it up a little to get a table we can use for further analyses

tag <- c("tag1", "tag2")
peaks_features_split <- separate(peaks_closest_features, 1, tag, sep = ";")

colnames(peaks_features_split) <- c("ID", "name", "distance")

peaks_feature_trim <- select(peaks_features_split, ID,  distance)

cleaned_features <- peaks_feature_trim %>%
  mutate_if(is.character, stringr::str_replace_all, pattern = 'ID=SVEN_', replacement = '')

#now we have a nice table and we want to remove duplicates (if the same gene is ID'd by all four pairwise comparisons, just keep the one with greatest distance). 

dedup_features <- unique(setDT(cleaned_features)[order(ID, -distance)], by = "ID")

ChIP_dat <- dedup_features %>% mutate(ID = as.numeric(ID)) %>% drop_na() 
```

Next, we will clean up RNA-sequencing data. Here we are turning GeneID’s
into a numeric type and assigning biosynthetic gene clusters to
diffrentially expressed genes.

``` r
#Cleaning up RNA-seq output ----

#ChIP_dat <- read.csv('~/BIOLOGY722/Data_Integration_Exercise/ChIP/ChIP_table_NEW.csv')
RNA_dat <- read.csv('~/BIOLOGY722/Data_Integration_Exercise/RNA-seq/RNA_seq_data.csv')
#removing 'SVEN' from geneid
RNA_dat <- RNA_dat %>%
  mutate_if(is.character, stringr::str_replace_all, pattern = 'SVEN_', replacement = '') %>%
  mutate(ID = as.numeric(ID)) %>% drop_na() 

RNA_clusters <- RNA_dat %>%
  mutate(Cluster = case_when((ID >= 223 & ID <= 234) ~ "Ectoine",
                               (ID >= 261 & ID <= 306) ~ "Terpene1",
                               (ID >= 463 & ID <= 531) ~ "T1PKS-T3PKS-NRPS",
                               (ID >= 540 & ID <= 561) ~ "Lantipeptide/terpene" ,
                               (ID >= 612 & ID <= 630) ~ "Lantipeptide(venezuelin)",
                               (ID >= 755 & ID <= 772) ~ "Indole(acryriaflavin)",
                               (ID >= 913 & ID <= 928) ~ "Chloramphenicol",
                               (ID >= 1844 & ID <= 1884) ~ "Other1",
                               (ID >= 2566 & ID <= 2577) ~ "Siderophore(desferrioxamine-like)",
                               (ID >= 3103 & ID <= 3132) ~ "Lassopeptide",
                               (ID >= 4061 & ID <= 4110) ~ "Other2",
                               (ID >= 4179 & ID <= 4189) ~ "Butyrolactone(gaburedin)",
                               (ID >= 4620 & ID <= 4662) ~ "Melanin1",
                               (ID >= 5076 & ID <= 5111) ~ "Butyrolactone",
                               (ID >= 5119 & ID <= 5145) ~ "Thiopeptide",
                               (ID >= 5351 & ID <= 5383) ~ "T3PKS1",
                               (ID >= 5413 & ID <= 5426) ~ "Siderophore1",
                               (ID >= 5471 & ID <= 5482) ~ "Siderophore2",
                               (ID >= 5817 & ID <= 5840) ~ "Bacteriocin1",
                               (ID >= 5951 & ID <= 6002) ~ "Butyrolactone/T2PKS",
                               (ID >= 6112 & ID <= 6204) ~ "Other3",
                               (ID >= 6134 & ID <= 6282) ~ "NRPS-ladderane",
                               (ID >= 6436 & ID <= 6490) ~ "Terpene2",
                               (ID >= 6527 & ID <= 6535) ~ "Bacteriocin2",
                               (ID >= 6767 & ID <= 6814) ~ "T2PKS",
                               (ID >= 6833 & ID <= 6842) ~ "Melanin2",
                               (ID >= 7032 & ID <= 7080) ~ "NRPS",
                               (ID >= 7101 & ID <= 7119) ~ "Terpene3",
                               (ID >= 7223 & ID <= 7259) ~ "T3PKS2",
                               (ID >= 7417 & ID <= 7452) ~ "Terpene-NRPS"))
```

Before proceeding to the any analysis we want to make sure that there
are no missing values in the tables we just created, so the following
code will assess this.

``` r
#Data Check ----

RNA_data<- RNA_dat %>% filter_all(any_vars(is.na(.))); RNA_data
CHIP_data<- ChIP_dat %>% filter_all(any_vars(is.na(.))); CHIP_data
```

#### Venn Diagram

Now that we have verified there is no missing data, we will start making
some plots.

To start our analysis we will make some Venn Diagrams that look at how
many genes are both differentialy expressed and bound at or near Lsr2.

``` r
#Venn diagram ----

Cluster_bound = RNA_clusters %>% filter(Cluster != is.na(Cluster)) %>% select(ID)

closest_100 <- ChIP_dat %>% filter(distance <= 100)

directly_bound <- ChIP_dat %>% filter(distance == 0)

venn.diagram(
  x = list(RNA_clusters$ID, directly_bound$ID),
  category.names = c("RNA_seq", "ChIP_seq"),
  filename = 'directly_bound.png',
  fill= c ("light blue", "pink"),
  output=TRUE
)

venn.diagram(
  x = list(RNA_clusters$ID, closest_100$ID),
  category.names = c("RNA_seq" , "ChIP_seq"),
  filename = 'closest_100.png',
  fill= c ("light blue", "pink"),
  output=TRUE)

venn.diagram(
  x = list(Cluster_bound[,1], RNA_clusters$ID, closest_100$ID),
  category.names = c("BGC", "RNA_seq" , "ChIP_seq"),
  filename = 'All_overlap.png',
  fill= c ("mediumorchid", "light blue", "pink"),
  output=TRUE)
```

#### Effect of Lsr2 binding on antibiotic production

Finally, we are going to do what we came here to do, integrate!

1.  To start we are going to assign either `Bound` or `Not Bound` status
    to the genes that are differentialy expressed.

2.  We are then plotting this binding status against the log2 Fold
    Change of the genes and assigning the color variable to clusters.

``` r
#Lsr2_bound or not

chip_id <- ChIP_dat %>% pull(ID)

lsr2_bound <- RNA_clusters %>% 
  mutate(Lsr2_bound = case_when(ID %in% chip_id ~ "Bound", 
                                TRUE ~ "Not Bound"))

pdf('Lsr2_bind_vs_log2FC.pdf')
ggplot(lsr2_bound, aes(Lsr2_bound, log2FoldChange)) + geom_violin() + geom_jitter(aes(color = Cluster), size = 2, alpha = 0.5) + theme_bw() + xlab('Lsr2 Binding Status') + ylab('log2 Fold Change') + geom_hline(yintercept = 0, alpha = 0.7, color = 'red') + ylim(min = -10, max = 10) 
dev.off()
```

Finally, generate `sessionInfo()` for a report on all the tools used
today.

``` r
sessionInfo()
```
