### The dada2 side

## Getting set up

# Load the required packages:

library(dada2)
packageVersion("dada2")
#library(ShortRead)
#packageVersion("ShortRead")
#library(ggplot2)
#packageVersion("ggplot2")
#library(tidyverse)
#packageVersion("tidyverse")
#library(plyr)
#packageVersion("plyr")
#library("viridis")
#packageVersion("viridis")


## The path

# Define a path variable to check it is all OK for the work we are going to do:

path <- ("/path/where/data/is/found/")
path
fns <- list.files(path)
fns

## Collecting our data

# extract out our fastq sequences
fastqs <- fns[grepl(".fastq$", fns)]
fastqs

# sort them to ensure reads are in the same order
fastqs <- sort(fastqs)
fastqs

# make sub-lists for the forward and reverse reads
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

# get the sample.names, customised for these file names
sample.names <- sapply(strsplit(fnFs, "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
sample.names

# specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
fnFs
fnRs


## Examine the quality profiles of forward and reverse reads

# It is important to look at your data. We start by visualizing the quality profiles along the sequencing reads
# Visualize the quality profile of the forward reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


## Perform filtering and trimming

# make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filt_path

# make list of filtered names for later
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240), maxN=0, 
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
out


## Dereplication

# dereplicate the forward and reverse reads separately
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# rename the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
names(derepFs)


## Error rates 

# forward reads first and then look at the output
start_time <- Sys.time()
errF <- learnErrors(derepFs, multithread=TRUE)
end_time <- Sys.time()
end_time - start_time

start_time1 <- Sys.time()
errR <- learnErrors(derepRs, multithread=TRUE)
end_time1 <- Sys.time()
end_time1 - start_time1

# plot the errors as a trellis plot
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


## set our wd()

path <- ("/path/where//our/output/data/is/going/")

## Sample inference and merging paired reads

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, BAND_SIZE=2)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, BAND_SIZE=2)

# inspect the results in more detail
dadaFs[[1]]
dadaRs[[1]]

# Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


## Constructing the sequence table and removing chimaeras

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
head(seqtab)

# look at the top 2 x 2 only as a check
seqtab[1:2, 1:2]

# extract the data as a text file
Tseqtab <- t(seqtab)
write.table(seqtab, file = "simpleCountsCryptoBS2.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# get the column names
Tnames <- c("sequence", paste("samp", colnames(seqtab), sep = "_"))
write.table(t(Tnames), file = "simpleColumnNamesBS2.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Inspect the distribution of sequence lengths
test <- table(nchar(getSequences(seqtab)))
test1 <- as.data.frame(test)
colnames(test1) <- c('seqLength', 'frequency')

# this is not totally correct, asit shows the number of sequences at each
# length rather than their abundance, but does show the length distribution
plot(test)

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab)

## Checking our progress

# Remove chimeric sequences with some complicated code
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

## Sequence output ##

# Use code found at https://github.com/benjjneb/dada2/issues/48
# load our required packages
library(phyloseq)
packageVersion("phyloseq")

# starting to make a dataframe for the samples by getting their names
samples.out <- rownames(seqtab.nochim)
samples.out

# R code
seqs <- colnames(seqtab.nochim)
otab <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
colnames(otab) <- paste0("seq", seq(ncol(otab)))
otab = t(otab)
write.table(seqs, "dada_seqs.txt", quote=FALSE)
write.table(otab, "dada_table.txt", quote=FALSE,sep="\t")
