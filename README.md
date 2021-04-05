# CHIP-seq and RNA-seq integration tutorial

#### Presented by: 

Meagan Archer, Meghan Pepler, Sreedevi Kesavan, and Stephanie Ali Fairbairn 
______________________________________________________
April 7th 2021

Welcome to the second tutorial for Biology 722! 

We are going to be spending this lesson:

1. Introducing **Streptomyces venezuelae** and the biological question we want to answer today.
2. Discussing the ChIP-seq analysis pipeline and how integrating these results into the RNA-seq data can give greater biological clarity to data. 
3. Workshop on ChIP-seq pipeline and integration.
4. Results and Discussion on biological question.
______________________________________________________

### Streptomyces venezuelae

Streptomyces are soil-dwelling Gram-positive bacteria that have a complex multi-stage lifecycle and secondary metabolism. They specialize in antibiotic production and are the organism for todays analysis. 

![](https://github.com/sk7-dotcom/Data_Integration_Exercise/blob/main/ChIP/Pictures/S_Ven.JPG)

#### Biological question for today

Antibiotic production in this bacteria is goverened by ~30 biosynthetic gene clusters. The products of these genes coordinate together to produce a specific antibiotic. Each cluster will encode: enzymes responsible for the catalysis and modification of the core antibiotic backbone, one or more gene product(s) to regulate the activity of the cluster, and resistance factors or efflux pumps to prevent cell suicide when the cluster is active. Lsr2 is a nucleoid associated protein that behaves similarly to Histone proteins in eucaryotic systems. It has the ability to bend and contort DNA to regulate gene expression. Here we investigate whether Lsr2 binding has an repressive effect on the biosynthetic clusters. 

![](https://github.com/sk7-dotcom/Data_Integration_Exercise/blob/main/ChIP/Pictures/Clusters.JPG)

### Chromatin Immunoprecipitation Sequencing (ChIP-Seq)

![](https://github.com/sk7-dotcom/Data_Integration_Exercise/blob/main/ChIP/Pictures/ChIP_pipeline.jpg)
[CHIP-seq pipeline](https://doi.org/10.1016/j.ymeth.2020.03.005)



### How to ChIP

![](https://github.com/sk7-dotcom/Data_Integration_Exercise/blob/main/ChIP/Pictures/chip_workflow.jpg)
[This](https://github.com/hbctraining/Intro-to-ChIPseq) and is a page for going into more detail on the ChIP-seq analysis process. This page details the process our experiment took, but for more detail on file types and even pipelines that are run entirely in the cluster go to HBC training website above.  

### How to Integrate 

