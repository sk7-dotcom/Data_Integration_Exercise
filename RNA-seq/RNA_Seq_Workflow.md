# Data_Integration_Exercise
## RNA-sequencing Differential Gene Expression Analysis
### Comparing Differentially Expressed Genes in Mutant and Wildtype *Streptomyces venezuelae* using R
##### Authors: Meagan Archer and Stephanie Ali Fairbairn

##### Date: 07/04/2021

##### Output: DESeq2 Differential Expression Results as .csv File

---

Welcome to workshop portion of today's tutorial!

This is part 1 of 3 sections towards integrating our data: RNA-seq analysis

For more context and details on each step look under the slides folder in the main repo. 

_________________

The raw reads have been previously quality checked, cleaned and aligned to the reference genome for *Streptomyces venezuelae* by a past member of the Elliot Lab. The next steps 
in the RNA-sequencing differential analysis pipeline are to (1) count the reads mapping to genomic features and (2) use the counts generated to test for differentially expressed 
genes in the mutant (T3d3225) compared to the wildtype (T3WT). We will be using the ```featureCounts()``` function of ```RSubread```, a package within the ```Bioconductor```
project package repository, to first assign the mapped sequencing reads to the genomic features contained in the *S. venezuelae* genome annotation (.gff) file. We will then
use the counts generated as input for the```DESeq2``` package, to test for differentially expressed genes. The workflow below assumes that ```Bioconductor```  and ```DESeq2```
have already been installed and that the necessary sorted .bam (containing the aligned reads) and genome annotation (.gff) files have been downloaded and saved to your working directory.

---

### (1) Use featureCounts() to assign mapped reads to genomic features


First, install the ```RSubread``` and ```ape``` packages. The ```ape``` package is needed to read in the .gff file. Please see here for the complete documentation for the [RSubread](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf) and 
[ape](https://cran.r-project.org/web/packages/ape/ape.pdf) packages. Load the associated libraries.
  
  
```r
# Install "RSubread" package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")

# Install "ape" package
install.packages("ape")

# Load the library for "RSubread" package
library(Rsubread)

# Load the library for "ape" package
library(ape)
```
  
  
Next, set your working directory using the ```setwd()``` function. Check that you are in the correct working directory using the ```getwd()``` function.
  
  
```r
# Set your working directory 
setwd("C:/Diff_Gene_Analysis")

# Check that you are in the correct working directory
getwd()
```
  
  
Once you have set up and are currently in your working directory, read the genome annotation (.gff) file into the environment as a data frame using the ```read.gff()``` function of the ```ape```
package and save the data frame into the variable, "s_venezuelae". This will allow you to check the structure of the genome annotation file. Pass the following arguments, ```file```, ```na.strings``` and ```GFF3``` to ```read.gff()```  as: ```"s_venezuelae.gff"```,
```c(".", "?")``` and ```TRUE```, respectively. The ```file``` argument specifies the annotation file. The ```na.strings``` argument specifies which characters will be turned into "NAs" in the final
format. The ```GFF3``` argument is a logical that specifies whether the annotation file is in .gff3 (```TRUE```) or .gff2 (```FALSE```) format. Perform a "quick check" by printing the first six rows of
the data frame using the ```head()``` function. 
  
  
```r
# Read the genome annotation file into the environment as a data frame and save it to the variable "s_venezuelae"
s_venezuelae <- read.gff("s_venezuelae.gff", na.strings = c(".", "?"), GFF3 = TRUE)

# Check the data frame structure

# Read the first six rows of the "s_venezuelae" data frame
head(s_venezuelae)
```


Depending on the annotation file, the structure of the data may differ slightly from the default parameters for ```featureCounts()```. ```featureCounts()``` expects the ```GTF.featureType``` to be ```"exon"```  and ```GTF.attrType``` to be ```"gene_id"``` by default.
```GTF.featureType``` denotes the feature type (in this case ```type```) used to select the rows of the genome annotation file for read summarization. ```GTF.attrType``` denotes
the attributes that will be used to group features into meta-features. 


By examining the structure of the "s_venezuelae" data frame, we can see that the ```GTF.attrType``` corresponds to ```"ID"```. This will need to be specified to ```featureCounts()``` in place of ```"gene_id"```.
  
  
Example Output ```head(s_venezuelae)```
```
seqid  source type start  end score strand phase
1 FR845719 GenBank gene    67  657    NA      -     1
2 FR845719 GenBank gene   166  396    NA      +     1
3 FR845719 GenBank gene   478  696    NA      +     1
4 FR845719 GenBank gene  1064 1678    NA      +     1
5 FR845719 GenBank gene  1987 2235    NA      +     1
6 FR845719 GenBank gene  2368 2784    NA      +     1
                   attributes
1 ID=SVEN_0001;Name=SVEN_0001
2 ID=SVEN_0002;Name=SVEN_0002
3 ID=SVEN_0003;Name=SVEN_0003
4 ID=SVEN_0004;Name=SVEN_0004
5 ID=SVEN_0005;Name=SVEN_0005
6 ID=SVEN_0006;Name=SVEN_0006
```


Since the ```GTF.featureType``` data in the annotation may differ, it is important to examine its corresponding column in the data frame. Use the ```unique()``` function
to search the data frame "s_venezuelae", by passing it as the first argument to the function ```unique()```. This will omit any duplicate rows in the ```type``` column, as specified by ```$```. 


```r
# Check the data in the "type" column of the "s_venezuelae" data frame
unique(s_venezuelae$type)
```


In this case, the only ```type``` data corresponds to ```gene``` and therefore there is no ```exon``` data. Therefore, ```GTF.featureType``` must be specified as ```gene```
to the ```featureCounts()``` function.


Example Output ```unique(s_venezuelae$type)```
```
[1] gene
Levels: gene
```
  
We will be counting the reads that are mapped to the genomic features of *S. venezuelae* for four previously sorted .bam files, with two sorted .bam files per genotype (mutant and wildtype), corresponding to two
biological replicates for each. First, specify the ```files``` argument to ```featureCounts()```. The ```files``` argument specifies the input files as a character vector. Aggregate the four sorted .bam
files using the ```dir()``` function. ```dir()``` specifies all files in a directory. To specify the current working directory, use, ```"."```, followed by ```"bam$"```, to select only the files ending in .bam.
Next, specify the annotation file to the ```annot.ext``` argument as ```"s_venezuelae.gff"```. Set the ```isGTFAnnotationFile``` to ```TRUE```, ```GTF.featureType``` to ```"gene"``` and
```GTF.attrType``` to ```"ID"```. Since we are not interested in summarizing the reads according to meta-features, set ```useMetaFeatures``` to ```FALSE```.

Check the first six rows of the outputed "sorted_counts" using ```head(sorted_counts)```.  
    
    
```r
# Use featureCounts to count the reads that are mapped to the genomic features and save the results into the variable "sorted_counts"
sorted_counts <- featureCounts(files=dir(".", "bam$"), annot.ext="s_venezuelae.gff", isGTFAnnotationFile=TRUE, GTF.featureType="gene", GTF.attrType="ID", useMetaFeatures = FALSE)

# Check the structure of "sorted_counts". Read the first six rows.
head(sorted_counts)
```
  
  
### (2) Use DESeq2 to perform differential gene expression analysis

#### A. Count Matrix Input

First, load the library for the ```DESeq2``` package.

```r
# Load the library for DESeq2
library(DESeq2)
```

To generate the ```DESeqDataSet``` for differential gene expression analysis, the ```DESeqDataSetFromMatrix()``` function requires the counts matrix (```countData```), the columns of the count matrix or sample information (```colData```) and the experimental design or design formula (```design```).


First, we must clean-up the "sorted_counts" output from ```featureCounts()```. We need to extract the ```$counts``` information from the "sorted_counts" output. To do this we will
subset the counts table from "sorted_counts" as a new data frame called "counts".


```r
# Subset the counts data from "sorted_counts"
counts <- sorted_counts$counts

# Check the first 6 rows of the new "counts" table.
head(counts)
```


Next, we must prepare the information corresponding to the columns of the count matrix or the sample information (```colData```). To do this we will create a data frame with
the names of each .bam file in one column and their corresponding genotype, "mutant" or "wildtype", in the second column.


```r
# Preparation of "colData"

# Create a vector called, "files", containing the names of the .bam files (will be column one)
files <- c("T3d3225v1.sorted.bam", "T3d3225v2.sorted.bam", "T3WTv1.sorted.bam", "T3WTv2.sorted.bam")

# Create a vector called, "genotype", containing the corresponding genotypes of the above .bam files (will be column two) 
genotype <- c("mutant", "mutant", "wildtype", "wildtype")

# Create the data frame "coldata"
coldata <- data.frame(files, genotype)
```


Now that we have prepared the input for the ```countData``` and ```colData``` arguments of the ```DESeqDataSetFromMatrix()``` function, we can generate the ```DESeqDataSet```.
Specifiy ```counts``` as the input for ```countData```, ```coldata``` as the input for ```colData``` and ```~genotype``` as the input for ```design```, as we are interested in
comparing differentially expressed genes between the mutant and wildtype.


```r
# Generate the DESeqDataSet, "deseqdata", for comparing differentially expressed genes in the mutant and wildtype
deseqdata <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~genotype)
```

#### B. Differential Expression Analysis


We will be using the standard ```DeSeq()``` function and associated functions to perform differential expression analysis. Since the last variable in our design formula above is
a factor, we must specify that the "wildtype" should be used to as the reference level for analysis. By default, ```DeSeq2``` will choose the reference level based on
whichever factor comes first in the alphabet. To assign the "wildtype" as the reference level, relevel the factors using the ```factor()``` function. Specify the data that we are
changing to the ```factor()``` function as the ```genotype``` in the ```DESeqDataSet``` using ```$```. Specify the levels of the factor by setting the ```levels``` argument
equal to a vector containing "wildtype" first and "mutant" second. Save this as "deseqdata$genotype".


```r
# Relevel the factors in "genotype" so that wildtype is the reference level for analysis
deseqdata$genotype <- factor(deseqdata$genotype, levels = c("wildtype","mutant"))
```


Next, we will pass the ```DESeqDataSet``` to ```DeSeq2``` for differential gene analysis. Output to the variable "deseqdata". Generate the results table using the function ```results()```. The results table will
contain the associated log2 fold changes, p-values and adjusted p-values.  Note that p-values are computed using the Wald test and that the default cut-off for the adjusted p-value is
alpha = 0.1.


```r
# Perform DeSeq() analysis and output to results table
deseqdata <- DESeq(deseqdata)

# Generate the results table
results <- results(deseqdata)

# View results table in the console
results
```


We will then organize the results according to the lowest adjusted p-value.


```r
# Reorders the results table such that the lowest adjusted p-value comes first and saves it to a variable called "resOrdered"
resultsOrdered <- results[order(results$padj),]

# View reordered results in the console
resultsOrdered
```


Finally, we will write the count data to a .csv file.


```r
# Write reordered results to a .csv file
write.csv(as.data.frame(resultsOrdered), 
          file="mutant_versus_wildtype_gene_expression_data.csv")
```
