# experimental data found in airway package contains
# contains RNAseq experiment with airway smooth muscle cells
# treated with dexamethasone (an anti-inflammatory steroid)

# 1.1
# Four cell lines 
# Human Smooth Muscle (treated/untreated)

# 2.1 note that to prepare count matrices - all statistical-based

# methods (e.g. DESeq2) need input (i,j) matrix 
# equals # of fragments for paired end reads OR 
# # of reads (single) OR binding regions (CHIP)

# For DESeq2 do not normalise reads prior 
# for a quicker version - use tximport and kallisto (transcript abundance quantification)

#2.2 Not done here is aligning reads to a reference genome
#FASTQ files (nucleotide sequence and quality score)
#using STAR to align to human reference genome

#This is the form of the for loop used 
##for f in `cat files`; do STAR --genomeDir ../STAR/ENSEMBL.homo_sapiens.release-75 \
# --readFilesIn fastq/$f\_1.fastq fastq/$f\_2.fastq \
# --runThreadN 12 --outFileNamePrefix aligned/$f.; done
# for f in `cat files`; do samtools view -bS aligned/$f.Aligned.out.sam \
# -o aligned/$f.bam; done

# 2.3 Locating alignment files 
#only reads aligned to small region of chromosome 1
library("airway", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

# informations
library("airway")
data(airway)
as.data.frame(colData(airway))
summary(colSums(assay(airway))/1e6)
metadata(rowRanges(airway))

#shws where the package has been installed 
indir <- system.file("extdata", package="airway", mustWork=TRUE)

#Find 8 BAM files
list.files(indir)

#usually a Csv file containin a sample that links 
#FASTQ and BAM files 
csvfile <- file.path(indir, "sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)

#once the reads have been aligned 
#tools assigned to enomic features of each sample e.g. gft file

#2.4 DESeq2 import functions
#Generate count matrices through summarizeOverlaps
filenames <- file.path(indir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)

#indicate BAM files (Rn interface) only 2mill reads at a time
#chromosome names of genomic features in annotation
#SAME as what is used for read alignment e.g. chr1 and 1 (use seqinfo(bamfiles[1]))
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

#2.5 Defining gene models
#read in the gene model used for counting reads/fragments
#using an Ensembl GTF file

# A TxDb object is a database that can be used to generate a variety of range-based objects, such as exons, transcripts, and genes. We want to make a list of exons grouped by gene for counting read/fragments.

#Fr known genes can use UCSC Genome browser 
library("GenomicFeatures")
gtffile <- file.path(indir,"Homo_sapiens.GRCh37.75_subset.gtf")
#none of the sequences are circular indicated by cicular 0-length character vector
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())

#list of all the exons grouped by genes 
ebg <-exonsBy(txdb, by="gene")

#2.6 Read counting step
#single core-  Expect that the summarizeOverlaps call will take at least 30 minutes per file for a human RNA-seq file with 30 million aligned reads. By sending the files to separate cores, one can speed up the entire counting process.
library("GenomicAlignments")
library("BiocParallel")


register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
#mode arguement supplies the kind of reads overlapped
#fragments counted once to each gene, even
#if they overlap only once to each gene
#if they overlap multiple exons of a gene which may themselves be overlapping

#singleEnd=FALSE produce paired end reads count only a pair of genes

#fragments arguement can be used when se=F

#if the RNAseq exp is strand specific https://www.ecseq.com/support/ngs/how-do-strand-specific-sequencing-protocols-work
#this experiment isn't so set ignorestand to TRUE
#experimentally, strand-specific protocol is achieved through attaching different ligators to either ends of the cDNA (this indicates directionality FRTseq)

#2.7 se means summarised experiment
#assay means a matrix of counts, row ranges contains genomic ranges and col data contains sample info

#use to check
dim(se) 
assayNames(se)# expected output counts
head(assay(se), 3)
colSums(assay(se))
rowRanges(se)#helps show the metadata
str(metadata(rowRanges(se))) #shows metadata compactly
colData(se) #emputy and should contain all metadata
#Because we used a column of sampleTable to produce the bamfiles vector, we know the columns of se are in the same order as the rows of sampleTable. We can assign the sampleTable as the colData of the summarized experiment, by converting it into a DataFrame and using the assignment function:
colData(se) <- DataFrame(sampleTable)
colData(se)

#2.8Branching point 
# this point we have counted fragments which overlap the genes in the gene model specified
# could use lots of other packages tfor exploration and differential expression analysis
# e.g. the limma point seen by law et al. 
# How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? paper can help devide which package to use

# 3. The DESeqDataSet object 
# certain data classes such as SummarisedExperiment needed for each package
# the core Bioconductor classes provide useful functionality: for example, subsetting or reordering the rows or columns of a SummarizedExperiment automatically subsets or reorders the associated rowRanges and colData, which can help to prevent accidental sample swaps that would otherwise lead to spurious results
# DESeqDataSet for DESeq2
# difference between DESeqDataSet and SummarisedExperiment
# 1) assay slot accessed by counts accessor funtion and non negative intgers
#  2) design formulae for DESeqDataset how to treat the samples in the analysis (one exception is the size factor estimation, i.e., the adjustment for differing library sizes, which does not depend on the design formula). The design formula tells which columns in the sample information table (colData) specify the experimental design and how these factors should be used in the analysis.

#simplest design formulae for DE is ~condition that specifies wich of two samples belong to
#airway experiment, specific ~ cell + dex 
se$cell
se$dex
library("magrittr")
se$dex %<>% relevel("untrt") #pipe operator also in R preferred first factor is untrt
se$dex


#steps we used to produce this object were equivalent to those you worked through in the previous sections, except that we used all the reads and all the genes.
data("airway")
se <- airway
se$dex %<>% relevel("untrt")
se$dex
round( colSums(assay(se)) / 1e6, 1 )

colData(se)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)

#3.2 Starting from count matrices 
#building DESeqDataSet with only a count matrix and a table of sample information
#see the paper

#4 Exploratory analysis and visualization
#two pathways 1) tranformation of counts to explore sample relationships 2) original raw counts for stats testing 
#critical- statistical cou

nrow(dds)
dds <- dds[ rowSums(counts(dds))>1, ]
# remove rows with only zeros 
nrow(dds)

#4.2 variance stabilizing transformaion and the rlog
#homoskedastic- expected amount of variance is approximately the same across different mean values
#PCA usually displays plots with the highest variance as the ones with the same highest counts as they show the largest differences. However if you take log of normalised counts and pseudocount of one, can counteract it. 
#however a pseudocount also means lowest counts create noises
lambda <- 10^seq(from=-1, to=2, length=1000)
#seq is a function that generates regular sequences
cts <- matrix(rpois(1000*100,lambda), ncol=100)
#rpois Density, distribution function, quantile function and random generation for the Poisson distribution with parameter lambda.
library("vsn")
meanSdPlot(cts, ranks= FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks= FALSE)

#two different distubutions
