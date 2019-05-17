#set the working directory
setwd("~/Documents/Informatics/AIH/cupcakes")

#install packages 
install.packages("DESeq2", dependencies = TRUE)
# if this doesn't work, try the following:
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

install.packages("biomaRt", dependencies = TRUE)

#load packages
library( "DESeq2" ) #gene expression analysis
library( "biomaRt" ) #converting gene names
library(ggplot2)
library(plyr)
vignette("DESeq2") #useful documentation on DESeq2

#load in gene counts file
genecounts<-read.delim(as.matrix("Cupcakes_AIH_gene_counts_final.txt"), header = TRUE, row.names = 1)
View(genecounts)

#make an ERCC gene counts dataframe and then remove those gene counts
ERCCgenecounts <- genecounts[grep("ERCC", rownames(genecounts)), ]
genecounts <- genecounts[ !rownames(genecounts) %in% rownames(ERCCgenecounts), ]

#load in metadata
metadata<-read.csv("Cupcakes_AIH_metadata.csv", header=TRUE, row.names = 1)

#check if rownames equal column names and are in the same order
all(rownames(metadata)==colnames(genecounts))

#check if metadata dataframe is filled
all(!is.na.data.frame(metadata))

#check if metadata is a factor
sapply(metadata,is.factor)

#take out version number of ENSEMBL ID
genenames<- sapply( strsplit( rownames(genecounts), split="\\." ), "[", 1 )
rownames(genecounts)<-genenames

#ensembl uses HG38 as a reference
listMarts(host="www.ensembl.org")
#create an object of class MART which describes the biomart, the dataset you are pulling from and the host
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
#view what you can filter by (similar to keytype)
listFilters(mart)
#query the BioMart database to pull the following atttributes, and to select the correct filter
#values are like keys, they are of the type 'filter
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "gene_biotype", "description"),
                  filters = "ensembl_gene_id",
                  values = genenames,
                  mart)

#create an index for mapping ensmbl gene ids from genemap into same order as they are in genenames
idx <- match(genenames, genemap$ensembl_gene_id )

#pull the variable names based on the idx
entrez <- genemap$entrezgene[ idx ]
hgnc_symbol <- genemap$hgnc_symbol[ idx ]
description <- genemap$description[ idx ]
gene_biotype <- as.character(genemap$gene_biotype[ idx ] )
ensembl <- genemap$ensembl_gene_id[ idx ]

ga<-as.data.frame(cbind(hgnc_symbol,description,gene_biotype,ensembl,entrez))

#pull of the of the genes which are protein_coding
pc<-ga[ga$gene_biotype=="protein_coding",]

#drop anything that says NA
pc<-pc[!(is.na(pc$ensembl)),]

#make an index for the protein coding genes and use it to create protein coding gene counts dataframe
idxpc<- match(pc$ensembl,rownames(genecounts))
genecounts_pc<-genecounts[idxpc,]

#turn the rownames into unique HGNC symbols (make.unique will append sequence numbers to duplicate HGNC symbols)
rownames(genecounts_pc)<-make.unique(as.character(pc$hgnc_symbol))
write.csv( as.data.frame(genecounts_pc), file="Cupcakes_AIH_genecounts-protein-coding.csv")

#PREFILTERING HERE

#view your metadata columns to view what you can choose for study design
colnames(metadata)

#create a DESeqDataSet from a matrix
dds<-DESeqDataSetFromMatrix(countData=genecounts_pc, colData= metadata, design=~Sex+AIH.Status)
?DESeqDataSetFromMatrix

#use default DESeq methods to estimate size factors, estimate dispersion and negative binomial GLM fitting/WALD statistics
dds<-DESeq(dds)
?DESeq

#NORMALIZATION AND QC HERE

#UNSUPERVISED CLUSTERING HERE: PCA plot
data <- plotPCA(vst(dds), intgroup=c("Sex", "AIH.Status", "AASLD.ID"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(as.numeric(PC1), as.numeric(PC2), color=AIH.Status)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


#extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, 
#standard errors, test statistics, p-values and adjusted p-values;
resDiag<-results(dds, contrast=c("AIH.Status", "Case_AIH", "Control_Healthy"), alpha=0.2)
summary(resDiag)

#pull out just FDR<0.2
resDiag.subset <- resDiag[ which(resDiag$padj < 0.2), ]
write.csv(resDiag.subset, file="Cupcakes_AIHvControl_FDR0.2_Genes.csv")
