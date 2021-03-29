library(dplyr)

RNA_dat <- read.csv('RNA-seq/mutant_versus_wildtype_gene_expression_data.csv')

RNA_dat_2 <- RNA_dat %>%
  mutate_if(is.character, stringr::str_replace_all, pattern = 'SVEN_', replacement = '') %>%
  mutate(geneid = as.numeric(geneid))
  
RNA_dat %>% 
  rowwise %>%
  summarise(NA_per_row = sum(is.na(.)))
  
RNA_dat_2[is.na(RNA_dat_2$padj), ]

map(RNA_dat_2, ~sum(is.na(.)))
