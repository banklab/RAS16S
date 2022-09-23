# Process PacBio raw reads (fastq.gz) using the DADA2 pipeline and do some basic plots

# Written by Adamandia Kapopoulou
# adamantia.kapopoulou@unibe.ch

# R 4.1.0 2021/05/19

# Install DADA2
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.13")

# Load libraries
library(dada2);packageVersion("dada2") # 1.20.0
library(Biostrings); packageVersion("Biostrings") # 2.60.1
library(ShortRead); packageVersion("ShortRead") # 1.50.0
library(ggplot2); packageVersion("ggplot2") # 3.3.4
library(reshape2); packageVersion("reshape2") # 1.4.4
library(gridExtra); packageVersion("gridExtra") # 2.3
library(phyloseq); packageVersion("phyloseq") # 1.36.0

# Set working directory where are stored all fastq.gz files
setwd("samples/")

# Create output folders
path.out <- "Figures/"
path.rds <- "RDS/"

fns <- list.files(pattern="fastq.gz", full.names=TRUE) # 11
lens.fn <- lapply(fns, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_all.", sep="", out.type="pdf"))
summary(lens)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   93    2357    2392    2326    2539    6692 
lens.fn <- lapply(fns[1], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_19H.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 341    1937    2564    2318    2706    4181 
lens.fn <- lapply(fns[2], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_20H.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 242    2106    2424    2309    2682    3582
lens.fn <- lapply(fns[3], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_4B.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 417    1721    2197    2118    2587    3702 
lens.fn <- lapply(fns[4], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_5B.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 170    2270    2422    2298    2443    5158 
lens.fn <- lapply(fns[5], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_5C.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 269    2172    2418    2278    2508    5514 
lens.fn <- lapply(fns[6], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_Weeks1.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 133    2378    2494    2361    2558    3773 
lens.fn <- lapply(fns[7], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_Weeks4.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 171    1804    2222    2142    2632    4336 
lens.fn <- lapply(fns[8], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_Mock.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 115    2378    2382    2399    2534    5174 
lens.fn <- lapply(fns[9], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_CtrlNeg.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 383    2429    2659    2562    2670    4712 
lens.fn <- lapply(fns[10], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_PLFilter1.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 93    1869    2388    2200    2613    5115 
lens.fn <- lapply(fns[11], function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_UNK.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 187    1749    2245    2111    2477    6692

# Primers used
F515 <- "AGRRTTYGATYHTDGYTYAG"
R806 <- "YCNTTCCYTYDYRGTACT"
# reverseComplement that handles the character vector DNA sequences generated in the dada2 package
rc <- dada2:::rc
# The classic dark-on-light ggplot2 theme. May work better for presentations displayed with a projector
theme_set(theme_bw())
# Remove primers and orient reads
path <- "."
nops <- file.path(path, "noprimers", basename(fns))
# Create output directory: noprimers/
prim <- removePrimers(fns, nops, primer.fwd=F515, primer.rev=dada2:::rc(R806), orient=TRUE)
# Inspect length distribution
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
dev.copy2pdf(file=paste("hist_lengths_noPrimers.", sep="", out.type="pdf"))
summary(lens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 53    2273    2342    2381    2486    5068 
filts <- file.path(path, "noprimers", "filtered", basename(fns))
# Create output directory: filtered/
# track <- filterAndTrim(nops, filts, minQ=3, minLen=500, maxLen=5500, maxN=0, rm.phix=TRUE, maxEE=2)
track
track <- filterAndTrim(nops, filts, minQ=2, minLen=500, maxLen=5500, maxN=0, rm.phix=TRUE, maxEE=3)
write.table(track, file="filterAndTrim_stats2.txt", quote=FALSE, sep="\t")
drp <- derepFastq(filts, verbose=TRUE)
err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
# 127791410 total bases in 54570 reads from 8 samples will be used for learning the error rates.
# Create output directory: RDS/
saveRDS(err, file.path(path.rds, "err.rds"))
plotErrors(err)
dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
saveRDS(dd, file.path(path.rds, "dd.rds"))
cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sapply(dd, function(x) sum(x$denoised)))
stats = cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sapply(dd, function(x) sum(x$denoised)))
write.table(stats, file="All_stats2.txt", quote=FALSE, sep="\t")
st <- makeSequenceTable(dd); dim(st) # 11 104
tax <- assignTaxonomy(st, "silva_nr_v128_train_set.fa.gz", multithread=TRUE)
tax[,"Genus"] <- gsub("Escherichia/Shigella", "Escherichia", tax[,"Genus"])
write.table(tax, file="seqs_tax2.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
head(unname(tax))
bim <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim)
# FALSE  TRUE 
#   87    17 
sum(st[,bim])/sum(st) # 0.02574814
sample.names <- gsub("./sfinderOutput.", "", fns)
sample.names <- gsub(".fastq.gz", "", sample.names)
rownames(st) <- sample.names
saveRDS(st, file.path(path.rds, "st.rds"))
saveRDS(tax, file.path(path.rds, "tax_Silva128.rds"))
# Reload processed data objects (can run code below from here in absence of input sequence data):
# st <- readRDS(file.path(path.rds, "st.rds"))
# tax <- readRDS(file.path(path.rds, "tax_Silva128.rds"))
# rsync -av RDS/* ak20x769@submit.unibe.ch:/storage/workspaces/vetsuisse_fiwi_mico4sys/fiwi_mico4sys001/PacBio/analysis/
ft <- sweep(st, 1, rowSums(st), "/")
df <- read.table("Metadata.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(df) <- df$SampleID
ps <- phyloseq(otu_table(ft, taxa_are_rows=FALSE), sample_data(df))
ord <- ordinate(ps, method="MDS", distance="bray")
dford <- cbind(df, ord$vectors)
ggplot(data=dford, aes(x=Axis.1, y=Axis.2, color=Location)) + geom_text(aes(label=SampleID))
dev.copy2pdf(file=paste("mds_bray2.", sep="", out.type="pdf"))
# Rarefy samples from each replicate down to a consistent read depth of 800 to remove systematic library size effects
# Set up file paths and define function
path.rare <- file.path(dirname(filts[[1]]), "Rare800")
if(!dir.exists(path.rare)) dir.create(path.rare)
# Define function that rarefies fastq files
rarefyFastq <- function(filt, fout, n) {
  require(ShortRead)
  if(file.exists(fout)) file.remove(fout)
  f <- FastqSampler(filt, n=n)
  fq.sample <- yield(f)
  close(f)
  nout <- writeFastq(fq.sample, fout)
  nout
}
# Seed RN generator for full replication
set.seed(100) 
rares <- file.path(path.rare, basename(filts))
nrec <- mapply(rarefyFastq, filts, rares, n=800)
names(nrec) <- sample.names
nrec
keep <- nrec == 800
# Now we go ahead and process the rarefied samples
drp.rare <- derepFastq(rares[keep])
dd.rare <- dada(drp.rare, err=err, multithread=TRUE, BAND_SIZE=32)
st.rare <- makeSequenceTable(dd.rare)
# We now merge the sequence tables into a data.frame with entries for each ASV found in each sample (in either replicate) along with the abundances detected in both replicates.k
sq.rare <- colnames(st.rare)
rownames(st.rare) <- sample.names[keep]
df <- melt(st.rare, varnames=c("Sample", "Sequence"), value.name="Abundance") # 385
# Remove the sequence/samples in which that sequence wasn't in that sample
df.rep <- df[df$Abundance > 0,] # 59
write.table(df.rep, file="table_seqs_after_rare2.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

stats = read.table("All_stats_long_format.txt", header=TRUE, stringsAsFactors=FALSE)
ggplot(stats, aes(x = Sample, y = Reads)) + geom_bar(aes(color = "#EFC000FF", fill = "#EFC000FF"), stat = "identity", position = position_stack())
dev.copy2pdf(file=paste("reads_raw.", sep="", out.type="pdf"))

ggplot(stats, aes(x = Sample, y = Reads)) + geom_bar(aes(color = Type, fill = Type), stat = "identity") +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF")) + scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) 
dev.copy2pdf(file=paste("reads_primers.", sep="", out.type="pdf"))

ggplot(stats, aes(x = Sample, y = Reads)) + geom_bar(aes(color = Type, fill = Type), stat = "identity") +
  scale_color_manual(values = c("#EFC000FF", "#FC4E07")) + scale_fill_manual(values = c("#EFC000FF", "#FC4E07")) 
dev.copy2pdf(file=paste("reads_filtered2.", sep="", out.type="pdf"))







