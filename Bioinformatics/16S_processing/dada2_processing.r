#this R script is largely following the dada2 tutorial available here: https://benjjneb.github.io/dada2/tutorial.html
#!/usr/bin/Rscript
#set_dependencies
library(dada2)
packageVersion("dada2") 
list.files() 


#r set create data and check QC 
#set paths for forward and reverse reads
path <- "trimmed"
fnFs <- sort(list.files(path, pattern="_1_trimmed.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_trimmed.fq.gz", full.names = TRUE))

# get the sample name from file pattern
sample.names <- sapply(strsplit(basename(fnFs), "_1_trimmed"), `[`, 1)

#some plots for sanity checks
pdf("QC_profile.pdf", height = 15, width = 15)
plotQualityProfile(fnFs) #problems to show all together
plotQualityProfile(fnRs) #more dramatic drop in reverse reads
dev.off()

#Filtering
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filtered_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, truncLen=c(270,200), rm.phix=TRUE, compress=TRUE, multithread=TRUE)
###keep filtered reads
saveRDS(filtered_out, "filtered_out.RData")
###filtering based on error rates
err_forward_reads <- learnErrors(filtFs, multithread=TRUE)
err_reverse_reads <- learnErrors(filtRs, multithread=TRUE)

#plot the error rates for sanity check
pdf("plotErrors.pdf", height = 15, width = 15)
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()

###Dereplication
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

###Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)

###Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, verbose=TRUE)

###keep filtered amplicons
saveRDS(merged_amplicons, "merged_amplicons.RData")

###Generating a count table
seqtab <- makeSequenceTable(merged_amplicons)
table(nchar(getSequences(seqtab)))

#remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.RData")

###get a little summary
# this is one quick way to look at sequences that have been lost, to know whether they held a lot in terms of abundance
sum(seqtab.nochim)/sum(seqtab)

#Overview of counts throughout:quick way to pull out how many reads were dropped at various points of the pipeline
# set a little function
getN <- function(x) sum(getUniques(x))

# make the summary table
summary_tab <- data.frame(row.names=sample.names, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
write.table(summary_tab, "summary_reads_table_microzebra_R220.txt", sep="\t", quote=F)


#Assign Taxanomy to ASVs
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", tryRC=T)
saveRDS(taxa, "taxa.RData")


###Write out tables for further processing}
#giving to seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
  }

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
