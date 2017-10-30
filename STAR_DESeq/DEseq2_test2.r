meta <- read.table("/data/project/HTseq/htseq_test2/deseq2_inputsummary.txt", header =TRUE)
# this summary table was generated from HTSeq run and contain sampleID, HTSeq count output file names and experimental conditions of your datasets.
# Note: the directory of these files is assigned specifically as DEseqDataSetFromHTSeqCount variable. Don't put path info here in the meta table.


# assign attributes of sampleTable based on meta table
samplenames <- as.character(meta$sampleID)
filenames <- as.character(meta$filenames)
samplecondition <- as.character(meta$condition)

sampleTable <- data.frame(sampleName = samplenames, fileName = filenames, condition = samplecondition)
print(sampleTable)


library(DESeq2) # Call DESeq2 in R and input HTseq count

directory <- "/data/project/HTseq/htseq_test2/" # assign directory of HTseq count output files

ddsHTSeq <- DESeqDataSetFromHTSeqCount (sampleTable = sampleTable, directory = directory, design = ~ condition)
print(ddsHTSeq)

dds <- DESeq(ddsHTSeq)
res <- results(dds)  
res05 <- results(dds,alpha=0.05) # highlight DE genes with FDR <0.05
summary(res05) 
sum(res05$padj <0.05, na.rm=TRUE) # how many DE genes meet the threshhold

plotMA(res05, ylim=c(-10,10)) # MAplot
jpeg("DEseq_res05_MAplot2.jpg")
plotMA(res05, ylim=c(-10,10))
dev.off()

res05Ordered <- res05[order(res05$pvalue),]
write.csv(as.data.frame(res05Ordered), file="parental_brain_FDR05_results.csv")
