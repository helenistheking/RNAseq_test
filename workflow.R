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
#