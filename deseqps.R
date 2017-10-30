# Differential Expression Problem Set - Jason Catolico

#Question 1

#Code to be run in the directory that contains the forward and reverse folders. 

# printf "\n%-27s %-9s %-8s %-10s %-9s\n" "File" "Direction" "Mapped" "No_feature" "% mapped"; for file in ./{forward,reverse}/mouse*; do printf "%-27s %-9s %-8s %-10s %-9s\n" ${file#./*/} ${file:2:7} $(cat $file | awk '/ENS/ {mapped+=$2} /no_feat/ {nofeat+=$2} END {print mapped,nofeat,100*mapped/(mapped+nofeat)}'); done
  
  # File                        Direction Mapped   No_feature % mapped
  # mouseFuramidine1.genecount  forward   1096717  23062795   4.53948
  # mouseFuramidine2.genecount  forward   1274810  25305061   4.79615
  # mouseFuramidine3.genecount  forward   1168628  22799936   4.87567
  # mouseHeptamidine1.genecount forward   1243208  26381265   4.50039
  # mouseHeptamidine2.genecount forward   1509915  32197122   4.47952
  # mouseHeptamidine3.genecount forward   992168   20239500   4.67306
  # mouseSaline1.genecount      forward   894086   18270692   4.66526
  # mouseSaline2.genecount      forward   1148623  23598904   4.64136
  # mouseSaline3.genecount      forward   816598   17358780   4.49288
  # mouseFuramidine1.genecount  reverse   18186803 5589572    76.4911
  # mouseFuramidine2.genecount  reverse   20351417 5761795    77.9353
  # mouseFuramidine3.genecount  reverse   18287042 5275395    77.611
  # mouseHeptamidine1.genecount reverse   20623191 6590358    75.7828
  # mouseHeptamidine2.genecount reverse   24986076 8216580    75.2532
  # mouseHeptamidine3.genecount reverse   16359172 4509673    78.3904
  # mouseSaline1.genecount      reverse   14442057 4397184    76.6594
  # mouseSaline2.genecount      reverse   18638226 5687597    76.6191
  # mouseSaline3.genecount      reverse   13521744 4371971    75.567
  
  # You'd want to use the reverse files since %mapped is way higher.
  
  #Question 2
  
  #install deseq2
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
  
  #begin pdf
  pdf(file="C:/cygwin64/home/Salmon/dea/reverse/results.pdf")  
  
  #get to the directory
  directory<-"C:/cygwin64/home/Salmon/dea/reverse"
  samplefiles<-grep("mouse",list.files(directory),value=TRUE) #files with mouse in the name
  list <- c("control1", "control2", "exp1", "exp2")
  samplecondition<-sub("mouse","\\1",samplefiles) #remove mouse from the name
  samplenames<-sub(".genecount","\\1",samplecondition) #remove .genecount from the name
  samplecondition<-sub("\\d.*","\\1",samplecondition) #remove the trial number from the condition
  sampletable <- data.frame(samplename=samplenames,filename=samplefiles,condition=samplecondition) #create table
  
  ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampletable,directory=directory,design= ~ condition) #import data 
  ddsHTSeq$condition <- factor(ddsHTSeq$condition,levels=c("Saline","Heptamidine","Furamidine")) #assign names to conditions
  dds<-DESeq(ddsHTSeq) #run deseq on the table
  
  #report changes from control
  resHept<-results(dds,contrast=c("condition","Heptamidine","Saline"), format="DataFrame") 
  resHeptOrdered <- resHept[order(resHept$padj),] 
  summary(resHept) #summary table
  plotMA(resHept, alpha=0.1, ylim=c(-3,3), main="Heptamidine vs Saline") #red stuff is above the significance level.
  plotCounts(dds, gene=which.min(resHept$padj), intgroup="condition", main = "Gene with minimum padj value under Heptamidine treatment")
  plotCounts(dds, gene=which.max(resHept$log2FoldChange), intgroup="condition", main = "Gene with greatest log2FoldChange under Heptamidine treatment")
  
  resFur<-results(dds,contrast=c("condition","Furamidine","Saline"), format="DataFrame") 
  resFurOrdered <- resFur[order(resFur$padj),] 
  summary(resFur) #summary table
  plotMA(resFur, alpha=0.1, ylim=c(-3,3), main="Furamidine vs Saline") #red stuff is above the significance level.
  plotCounts(dds, gene=which.min(resFur$padj), intgroup="condition", main = "Gene with minimum padj value under Furamidine treatment")
  plotCounts(dds, gene=which.max(resFur$log2FoldChange), intgroup="condition",main = "Gene with greatest log2FoldChange under Furamidine treatment")
  
  #create PCA plot
  rld <- rlog(dds, blind=FALSE)
  plotPCA(rld, intgroup="condition")
  
#The PCA plot shows that one Heptamidine sample isn't clustering with the rest. Possible outlier.
  
  #create heatmap
  sampleDists <- dist(t(assay(rld)))
  library("RColorBrewer") 
  install.packages("pheatmap")
  library(pheatmap)
  library(ggplot2)
  sampleDistMatrix <- as.matrix(sampleDists) 
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
  colnames(sampleDistMatrix) <- NULL 
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

  
# The heatmap shows that one of the Heptamidine samples matches Furamidine sample profile more closely than the other samples. We must remove the sample.
  

  sampletable <- sampletable[-6,]

  ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampletable,directory=directory,design= ~ condition) #import data 
  ddsHTSeq$condition <- factor(ddsHTSeq$condition,levels=c("Saline","Heptamidine","Furamidine")) #assign names to conditions
  dds<-DESeq(ddsHTSeq) #run deseq on the table
  
  #report changes from control
  resHept<-results(dds,contrast=c("condition","Heptamidine","Saline"), format="DataFrame") 
  resHeptOrdered <- resHept[order(resHept$padj),] 
  summary(resHept) #summary table
  plotMA(resHept, alpha=0.1, ylim=c(-3,3), main="Heptamidine vs Saline") #red stuff is above the significance level.
  plotCounts(dds, gene=which.min(resHept$padj), intgroup="condition", main = "Gene with minimum padj value under Heptamidine treatment")
  plotCounts(dds, gene=which.max(resHept$log2FoldChange), intgroup="condition", main = "Gene with greatest log2FoldChange under Heptamidine treatment")
  
  resFur<-results(dds,contrast=c("condition","Furamidine","Saline"), format="DataFrame") 
  resFurOrdered <- resFur[order(resFur$padj),] 
  summary(resFur) #summary table
  plotMA(resFur, alpha=0.1, ylim=c(-3,3), main="Furamidine vs Saline") #red stuff is above the significance level.
  plotCounts(dds, gene=which.min(resFur$padj), intgroup="condition", main = "Gene with minimum padj value under Furamidine treatment")
  plotCounts(dds, gene=which.max(resFur$log2FoldChange), intgroup="condition",main = "Gene with greatest log2FoldChange under Furamidine treatment")

#create PCA plot
  rld <- rlog(dds, blind=FALSE)
  plotPCA(rld, intgroup="condition")

#We've removed the outlier and the PCA plot looks better.
  
#create heatmap
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists) 
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
  colnames(sampleDistMatrix) <- NULL 
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
#The heatmap looks better now.
  
  dev.off() 