#Data Integration

# Packages to install ----

#BiocManager::install("igvR")

#Packages used: ----
library(tidyverse)
library(VennDiagram)
#library(igvR)


#Datasets: ----
RNA_dat <- read.csv('mutant_versus_wildtype_gene_expression_data_BETTER.csv')
ChIP_dat <- read.csv('ChIP/ChIP_table.csv')
cluster_chip_data <- read.csv('ChIP/ChIP_table_clusters.csv')
cluster_data <- read.csv('RNA_table.csv')

#Data Checks: ----

#Matching columns across data
ChIP_dat <- ChIP_dat %>% rename(ID = geneid)
cluster_data <- cluster_data %>% rename(geneid = ID)

#checking for NAs in datatables 
RNA_data<- RNA_dat%>% filter_all(any_vars(is.na(.))) %>% arrange(desc(baseMean))
CHIP_data<- ChIP_dat%>% filter_all(any_vars(is.na(.)))

#removing 'SVEN' from geneid
RNA_dat <- RNA_dat %>%
  mutate_if(is.character, stringr::str_replace_all, pattern = 'SVEN_', replacement = '') %>%
  mutate(geneid = as.numeric(X)) %>% drop_na() 

ChIP_dat <- ChIP_dat %>% mutate(geneid = as.numeric(geneid)) %>% drop_na(geneid) 

#joining data based on CHIP table 
join_try <- ChIP_dat %>% inner_join(RNA_dat, by = 'geneid'); join_try

#checking joint table for NAs
join<- join_try%>% filter_all(any_vars(is.na(.))) %>% arrange(desc(geneid))

#Venn attempt ----

ChIP <- ChIP_dat %>% select(geneid); RNA <- RNA_dat_2 %>% select(geneid)

fit <- euler(c(RNA = RNA[,1], CHIP = ChIP[,1]),
             shape = "ellipse") # <--- keeps aborting r session so I am not doing this anymore 

#clearly all the gene from ChIP are in the RNA set, so thats nice. This can be used for future 
#union efforts. 



venn.diagram(
  x = list(RNA[,1], ChIP[,1]),
  category.names = c("RNA_seq" , "CHIP_seq"),
  filename = 'venn_diagramm.png',
  output=TRUE
)

#Lsr2_bound or not
chip_id <- ChIP_dat[,1]

plot_try <- cluster_data %>% 
  mutate(Lsr2_bound = case_when(ID %in% chip_id ~ "Bound", TRUE ~ "Not Bound"))
  
#Preliminary plots ----

ggplot(join_try, aes(distance, log2FoldChange)) + geom_point()

plot_1 <- ggplot(plot_try, aes(ClusterID, log2FC, color = Lsr2_bound)) + 
  geom_jitter() + coord_flip() + theme_bw() + geom_hline(yintercept = 0) +
  ylim(min = -10, max = 10)








