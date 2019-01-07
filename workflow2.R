#4.3 Sample Distances
#test the similarity between samples
#fit with hypothesis Euceldian distance 

sampleDists <- dist(t(assay(vsd)))
#t transposes the matrix
#vsd is a DESeq2 object with gene as rows and samples as columns

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists) 
rownames(sampleDistMatrix) <- paste (vsd$dex, vsd$cell, sep= " - ")
colnames(sampleDistMatrix) <- NULL
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_cols = sampleDists, clustering_distance_rows= sampleDists, col=colours)

#produces a heatmap of sample-to-sample distanes using rlog-transformed values
#row names of distance matrix ontains treatment tye nd patient ID
#could repeat using PoiClaClu package (see Love paper)

#4.4 PCA Plot
plotPCA(vsd, intgroup = c("dex", "cell"))
#intgroup arguement for labelling the samples

#plot pca using ggplot
 pcaData <-plotPCA(vsd, intgroup= c("dex", "cell"), returnData= TRUE)
 pcaData
 
 #round up the percentVar 
percentVar <- round(100*attr(pcaData, "percentVar"))
#attr function get or set specific attributes of an object

#use this data to build up a second plot, 
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

# 4.5 mds plot
# similar to PCA plot when we only have a matrix of distances not matrix of data (as usual)
#this can be repeated with poisson
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

#5 Differenial expression analysis
#running the differential expression pipeline
dds <- DESeq(dds)
#run the differential expression pipeline on the raw counts with a single call to the function DESeq:
#using pre-existing size factors, estimating dispersions, gene-wise dispersion estimates, mean-dispersion relationship, final dispersion estimates, fitting model and testing

#5.2 building the results table
#calling results extract estimated log2fold changes and p values for last variable in design formula
#comparison is printed at the top
res <- results(dds)

res <- results(dds, contrast=c("dex","trt","untrt"))
#makes lfcSE 

mcols(res, use.names = TRUE)
#first column, baseMean, is a just the average of the normalized count values, divided by the size factors, taken over all samples in the DESeqDataSet. The remaining four columns refer to a specific contrast, namely the comparison of the trt level over the untrt level for the factor variable dex.
#log2foldchange is base 2
#lrcSE (standard errr of log fold change) (is it false positve?)
#pvalue

summary(res)
#adj pvalue of less than 0.1 
#Note that there are many genes with differential expression due to dexamethasone treatment at the FDR level of 10%.
#This makes sense, as the smooth muscle cells of the airway are known to react to glucocorticoid steroids. 
#lower the false discovery rate threshold (the threshold on padj in the results table)
#raise the log2 fold change threshold from 0 using the lfcThreshold argument of results

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

#DESeq’s way of reporting that all counts for this gene were zero, 

#5.3 other comparisons
#find othe analysises - name f variablem name of level of numerator and the denmiator
results(dds, contrast = c("cell", "N061011", "N61311"))
#to upload onto the main comp
files <- c("workflow2.R", "Rplot01.jpeg")
scp_upload(session, files, to = "RNAseq_test")
