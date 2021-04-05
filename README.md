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

#### [What is ChIP-Seq?](https://www.illumina.com/techniques/sequencing/dna-sequencing/chip-seq.html)

By combining chromatin immunoprecipitation (ChIP) assays with sequencing, ChIP sequencing (ChIP-Seq) is a powerful method for identifying genome-wide DNA binding sites for transcription factors and other proteins. Following ChIP protocols, DNA-bound protein is immunoprecipitated using a specific antibody. The bound DNA is then coprecipitated, purified, and sequenced.

![](https://github.com/sk7-dotcom/Data_Integration_Exercise/blob/main/ChIP/Pictures/ChIP_pipeline.jpg) Figure Can be found [here](https://www.nature.com/articles/sdata201510#Fig1)

The application of next-generation sequencing (NGS) to ChIP has revealed insights into gene regulation events that play a role in various diseases and biological pathways, such as development and cancer progression. ChIP-Seq enables thorough examination of the interactions between proteins and nucleic acids on a genome-wide scale.

### How to ChIP

![](https://github.com/sk7-dotcom/Data_Integration_Exercise/blob/main/ChIP/Pictures/chip_workflow.jpg)

[This](https://github.com/hbctraining/Intro-to-ChIPseq) is a page for going into more detail on the ChIP-seq analysis process. This page details the process our experiment took, but for more detail on file types and even pipelines that are run entirely in the cluster go to HBC training website above.  

#### Tools used in this tutorial (does not incluse packages for RNA-seq pipeline)

Most tools have been run with [Bioconductor](https://bioconductor.org/) tools but can just as easily be run with shell scripting. Look into the HBC training page for alternative code to several tools listed below:

1. [Rbowtie2](https://bioconductor.org/packages/release/bioc/manuals/Rbowtie2/man/Rbowtie2.pdf)
2. [Rsamtools](https://bioconductor.org/packages/release/bioc/manuals/Rsamtools/man/Rsamtools.pdf)
3. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. [MACS3](https://pypi.org/project/MACS3/) ([MACSr](https://52.71.54.154/packages/devel/bioc/manuals/MACSr/man/MACSr.pdf) is supposed to come out soon for the Bioconductor environment.)
5. [BedTools](https://bedtools.readthedocs.io/en/latest/) or [bedr](https://cran.r-project.org/web/packages/bedr/bedr.pdf)
