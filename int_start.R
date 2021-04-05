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
fold_enrich <- read.csv('ChIP/chip_fe.csv')
bound_genes <- read.csv('ChIP/bound_genes.csv')
EL_data <- read.csv('ChIP/RNA_plot_EL.csv')
EL_chip <- read.csv('ChIP/chip_EL.csv')

#Data Checks: ----

#Matching columns across data
ChIP_dat <- ChIP_dat %>% rename(geneid = ID)
cluster_data <- cluster_data %>% rename(geneid = ID)
chip_fe <- fold_enrich %>% rename(geneid = ID)
RNA_dat <- RNA_dat %>% rename(ID = X)

#checking for NAs in datatables 
RNA_data<- RNA_dat%>% filter_all(any_vars(is.na(.))) %>% arrange(desc(baseMean))
CHIP_data<- ChIP_dat%>% filter_all(any_vars(is.na(.)))

#removing 'SVEN' from geneid
RNA_dat <- RNA_dat %>%
  mutate_if(is.character, stringr::str_replace_all, pattern = 'SVEN_', replacement = '') %>%
  mutate(ID = as.numeric(ID)) %>% drop_na() 

ChIP_dat <- ChIP_dat %>% mutate(geneid = as.numeric(geneid)) %>% drop_na(geneid) 



#Venn attempt ----

Cluster_bound = RNA_clusters %>% filter(Cluster != is.na(Cluster)) %>% select(geneid)

venn.diagram(
  x = list(EL_data$geneid, EL_chip[,1], Cluster_bound[,1]),
  category.names = c("BGC", "RNA_seq" , "ChIP_seq"),
  filename = 'venn_diagramm.jpg',
  output=TRUE
)

#Lsr2_bound or not

chip_id <- ChIP_dat[,1]

plot_try <- RNA_clusters %>% 
  mutate(Lsr2_bound = case_when(geneid %in% chip_id ~ "Bound", TRUE ~ "Not Bound"))

pdf('ChIP/Pictures/Lsr2_bind_vs_log2FC.pdf')

ggplot(plot_try, aes(Lsr2_bound, log2FoldChange, color = Cluster)) + geom_jitter(size = 2, alpha = 0.5) + 
  theme_bw() + xlab('Lsr2 Binding Status') + ylab('log2 Fold Change') + geom_hline(yintercept = 0, alpha = 0.7, color = 'red') + 
  ylim(min = -10, max = 10)

dev.off()
  
#Preliminary plots ----



# Modeling

#joining data based on CHIP table 
join_try <- fold_enrich %>% inner_join(RNA_clusters, by = 'geneid'); join_try
distance <- ChIP_dat %>% inner_join(RNA_clusters, by = 'geneid'); distance


ggplot(distance, aes(distance, log2FoldChange, color = cluster)) + geom_point()

#checking joint table for NAs

fold_enrich <- fold_enrich %>%
    mutate(geneid = as.numeric(geneid)) %>% drop_na() 

cluster_data <- cluster_data %>%
  mutate(geneid = as.numeric(geneid)) %>% drop_na() 


join<- join_try%>% filter_all(any_vars(is.na(.))) %>% arrange(desc(geneid))

hist(join_try$score)

m1 <- lm(log2FoldChange ~ signal, data = join_try); summary(m1)
plot(m1)

m2 <- glm(log2FoldChange ~ score, data = join_try, family=gaussian(link="identity"))

library(ggplot2)

ggplot(join_try, aes(signal, log2FoldChange, color = cluster)) + geom_point()


###########################

RNA_clusters <- RNA_dat %>%
  mutate(Cluster = case_when((ID >= 223 & ID <= 234) ~ "Ectoine",
                               (ID >= 261 & ID <= 306) ~ "Terpene1",
                               (ID >= 463 & ID <= 531) ~ "T1PKS-T3PKS-NRPS",
                               (ID >= 540 & ID <= 561) ~ "Lantipeptide/terpene"	,
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

RNA_clusters <- RNA_clusters %>% rename(geneid = ID)



