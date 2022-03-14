# hello-world [FinalDraft.pdf](https://github.com/leylyx/hello-world/files/8241596/FinalDraft.pdf)
Click on the FinalDraftLink to see the code fully run with graphs and final visualization data
Copepod gene expression in Santa Cruz vs. San Diego. 
There are 4 groups, Control Santa Cruz, Control San Diego, Heat Stress Santa Cruz, Heat Stress San Diego
## Question: How will seeing the changes in gene expression of the heat shock proteins between San Diego population of copepods and Santa Cruz populations of copepods under acute thermal stress show us how the heat shock protein contributes to the adaptation of SD copepods to higher water temperatures?

## This data set was taken from a paper that sequenced the RNA of two different populations of fish living in different water temperatures, one population in Santa Cruz and one population in San Diego and tested mRNA expression of genes in control temperature and heat stress temperature to test how the populations react differently to heat (Schoville et al., 2012). The data in this coding was taken from the mRNA sequences that they sequenced from the fish that they inserted into NCBI and GEO dataset of the differential expression in RPKM of the significant genes they found.The heat shock protein is a key protein in preventing protein degradation because they help refold proteins and can be upregulated when needed under stressful situations (Mayer et al.,2005). There's a wide range of functions for heat shock proteins, including preventing cell aggreggation and preventing denaturation which can occur when a species undergoes heat stress (Mayer et al., 2005). These proteins are a gene of interest because San Diego populations of copepods are acclimated to higher water temperatures than Santa Cruz fish which could possibly be due to differences in their heat shock protein gene expressions. 

## The main hypothesis going into this paper is then: if copepods in San Diego have higher heat tolerance than copepods in Santa Cruz, then comparing gene expression changes of their heat shock proteins while they both undergo acute thermal stress will be able to show key differences in gene expression that account for their increased heat tolerance. 

## Pairwise sequence alignment of the initial mRNA of the heat shock protein sequences from control SD and control SC groups were compared to see if there differences in base sequence. This information was taken from NCBI. Then differential expression analysis was done with the GEO data withe RPKM reads of the mapping data of the 4 groups tested. The volcano plot, log expression changes and heat map was used to visualize the results of the differential expression analysis, to see the significantly changed genes between the 4 groups.   
```

```{r}

```{r}
## These libraries read the RNA / DNA off the fasta files and analyze them as DNA strings instead of just letters.
library(Biostrings)
 library(seqinr)
##Annotate library helps get data out of the files we are going to use 
BiocManager::install("annotate")
library(annotate)
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
## Readr allows .txt files to be loaded into the R studio
#read. delim function brings in the txt file
library(readr)
library(tidyverse)
## Libraries required for differential expression analysis 
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(DESeq2)
if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 
 BiocManager::install("SummarizedExperiment")

```


```{r}
## Reads the fasta file, skipping the first line with the name and looking at the sequence and analyzing it 
sc_RNA <- read.fasta(file = "~/Downloads/sequence (4).fasta")[[1]]
```
```{r}
## RNA stringset function reads the sequence as RNA and analyzes it instead of just reading random letters in a file. it reads them as amino acids
sc_RNA_Stringset <- readRNAStringSet("~/Downloads/sequence (4).fasta", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE, with.qualities=FALSE)
```


```{r}
## Description of method: Pairwise alignment looks at the nucleotide sequences of separate strings and assesses how similar they are in term of sequence by aligning them and finding areas of similar sequence. 
## Read RNA string set, reads the code as RNA instead of reading the letters, and assigns it as an RNA sequence. 
sd_RNA <- read.fasta(file = "~/Downloads/sequence (5).fasta")[[1]]
sd_RNA_stringset <- readRNAStringSet("~/Downloads/sequence (4).fasta", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE, with.qualities=FALSE)
```

```{r}
## Dot plot analyzes the differences between the sequences and shows the similarity in sequence
dotPlot(sd_RNA, sc_RNA, wsize=3, wstep=5, nmatch=3)
```
```{r}
## Biostrings library helps to load in the substitution matrix and Pairwise alignment code. Biostrings is a subset of a Bioconductor function, Nucleotide substitution matrix assigns scores to alignments of the nucleotides so that we can get a quantifiable score on how aligned or misaligned it is. 
library(Biostrings)
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
sigma 

```
```{r}
## Aligning the mRNA sequences of SD heatshock mRNA with SC heat shock mRNA with pairWise alignment functioning, and using the sigma Nucleotide substitution matrix we created. 
globalAligns1s2 <- pairwiseAlignment(sc_RNA_stringset, sd_RNA_stringset, substitutionMatrix = sigma, gapOpening = -2,
gapExtension = -8, scoreOnly = FALSE)

globalAligns1s2
```
```{r}
## There isn't a lot of changes in the mRNA of the heat shock proteins itself so there must be changes in level of gene expression. 


```{r}
## I'm reading in the mapping files that I got from the GEO database. The grep control is to filter out all of the datapoints that are not heatshock proteins by applying a filter for "heat shock" in the predicted products column of the gene mapping. 
SD_control_mapping <- read.delim("~/Downloads/GSE38546_RAW/GSM945350_SD_control_mapping 3.txt" )
SD_control_mapping
```
```{r}
SD_control_mapping <- SD_control_mapping[grep("heat shock",SD_control_mapping$Predicted.product),]
SD_control_mapping
```
SC_control_mapping <- read.delim("~/Downloads/GSE38546_RAW/GSM945352_SCN_control_mapping.txt")
SC_control_mapping[grep("heat shock",SC_control_mapping$Predicted.product),]
SC_control_mapping <- SC_control_mapping[grep("heat shock",SC_control_mapping$Predicted.product),]
```{r}
SD_heatshock_mapping <- read.delim("~/Downloads/GSE38546_RAW/GSM945351_SD_heatshock_mapping.txt")
SD_heatshock_mapping <- SD_heatshock_mapping[grep("heat shock",SD_heatshock_mapping$Predicted.product),]
```
```{r}
SD_heatshock_mapping
```
```{r}
SC_heatshock_mapping <- read.delim("~/Downloads/GSE38546_RAW/GSM945353_SCN_heatshock_mapping.txt")
SC_heatshock_mapping <- SC_heatshock_mapping[grep("heat shock",SC_heatshock_mapping$Predicted.product),]
# Taking the log of the RPKM to normalize the data of each data set and then plotting the controls and the heatshock against each other to see the different expression analysis and statistically up regulated or down regulated genes in each location. 
log_SC_heatshock_mapping <- log(SC_heatshock_mapping[,4])
```
```{r}
log_SC_control_mapping <- log(SC_control_mapping[,4])
```
```{r}
plot(log_SC_control_mapping, log_SC_heatshock_mapping)
```
log_SD_heatshock_mapping <- log(SD_heatshock_mapping[,4])
log(SD_heatshock_mapping[,4])
log_SD_control_mapping <-  log(SD_control_mapping[,4])
plot(log_SD_control_mapping, log_SD_heatshock_mapping)
## When we compare the graphs, we see that the SD graph has significantly more differentially expressed outliers that venture away from the group cluster whilst the SC graphs have a more equal expression when undergoing heat versus. control. This shows that in the San Diego fish, their heat shock protein is more differentially expressed especially when undergoing heat stress indicating that their heat shock proteins function differently and respond differently to heat. 


## DESEQ Analysis description: We have RPKM or Reads per kilobase of transcript, per million mapped reads. DESeq analysis of these values are going to find the proteins that are more statistically significant, find count outliers, and helps with the visualization process.    
```{r}
## The format of the dataset is not ideal for using DESeq, so compile all the row 4's of the four datasets that have the RPKM values into one matrix
data.set <- cbind(SD_control_mapping[,4], SD_heatshock_mapping[,4] , SC_control_mapping[,4], SC_heatshock_mapping[,4])

colnames(data.set) <- c("SDControl" , "SDHeatshock" , "SCControl" , "SCHeatshock")
```


```{r}
## Setting the rownames to the be the name of the protein product to align with the RPKM of each fish group and assigning it into a data matrix
rownames(data.set) <- (SC_control_mapping$Predicted.product)
data.matrix <- as.matrix(data.set)
```


```{r}
## Flipping the orientation of the data matrix so that ncol(countdata) == nrow(ColDATA)
countdata <- t(data.matrix)
countdata 
```
## Assigning the condition of each dataset to either "heatshock" or "control" based on the treatment of the fish 
metadata_SD_heatshock <- 
 SD_heatshock_mapping %>% 
  mutate_at("GenBank.Accession", str_replace, "waiting assignment", "heatshock")
metadata_SD_heatshock

metadata_SD_control <- 
 SD_control_mapping %>% 
  mutate_at("GenBank.Accession", str_replace, "waiting assignment", "control")
metadata_SD_control

```

```{r}
metadata_SC_heatshock <- 
 SC_heatshock_mapping %>% 
  mutate_at("GenBank.Accession", str_replace, "waiting assignment", "heatshock")
metadata_SC_heatshock
```
```{r}
SC_control_mapping
metadata_SC_control <- 
 SC_control_mapping %>% 
  mutate_at("GenBank.Acccession", str_replace, "waiting assignment", "control")
metadata_SC_control
##Re-arranging the row order so that col names match row names of the metadata and count data
metadata_SD_heatshock <-metadata_SD_heatshock[,c(3,1,2,4)]
metadata_SD_heatshock
```
metadata_SC_heatshock <-metadata_SC_heatshock[,c(3,1,2,4)]
metadata_SC_heatshock
```

```{r}
metadata_SD_control <-metadata_SD_control[,c(3,1,2,4)]
metadata_SD_control
```

```{r}
meta_data <- as.matrix(metadata_SC_control)
```

```{r}
## Counts data and meta_data match 
ncol(countdata)
nrow(meta_data)
```
##Creating the DESeq object 
colnames <- as.matrix(meta_data)
dds <- DESeqDataSetFromMatrix(countData = round(countdata), colData = metadata_SC_control, design = ~ Predicted.product , tidy=FALSE)
```


```{r}
## Running Differential Expression analysis on the DESeq matrix
dds <- DESeq(dds)
```
## Displaying the results 
res <- results(dds)
head(results(dds, tidy=FALSE))


```
## Summary of results
summary(res)
```
## Sort the summary list by P- value
res <- res[order(res$padj),]
head(res)

## Creating the volcano plot 
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
## Log2Fold Change vs -Log10P value shows plot points that are not highly differentially expressed but the most statistically significant of those differentially expressed data points. 
##Normalize the count matrix
dds_wt <- estimateSizeFactors(dds)
sizeFactors(dds_wt)
## Unsupervised clustering 
vsd_wt <- varianceStabilizingTransformation(dds_wt)
## Extract variance stabilizing transformation matrix
vsd_mat_wt <- assay(vsd_wt)
## Pairwise correlation value 
vsd_cor_wt <- cor(vsd_mat_wt)
View(vsd_cor_wt)
```

```{r}
plotPCA(vsd_wt, intgroup="Predicted.product")
```
## Create the heat map
library("RColorBrewer")
heatmap(vsd_cor_wt)
library(pheatmap)

# Plot heatmap

pheatmap(vsd_cor_wt, annotation = select.list(metadata_SC_control))
## Result analysis: From the pairwise sequence alignment, there wasn't a lot of misalignment or changes in the mRNA of the heat shock protein sequences between the control SD and control SC copepod groups. After doing the log change of the RPKM we saw that SD populations had more change in gene expression than SC population after going through heat shock. After doing the DESeq to analyze the RPKM reads, the volcano plot showed that there only a couple of values that had a signficant large fold change that were also statistically significant.  The heat map shows that there a different patterns of expressiong amongst various heat shock proteins but some were not differentially expressed when put through heat stress in SD vs. Santa Cruz. The results of this study overall show that the heat shock proteins could be of significance in helping these different populations adapt to different water temperatures but a lot of them were not significantly expressed differently to be considered a contendor for helping them with this adaptation.  


