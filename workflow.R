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
