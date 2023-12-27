#read in packages
library(tidyverse)
library(decontam)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(iNEXT); packageVersion("iNEXT")
library(ape); packageVersion("ape")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(cowplot); packageVersion("cowplot")
library(plyr); packageVersion("plyr")
library(dplyr)
library(sjmisc); packageVersion("sjmisc")
library(data.table); packageVersion("data.table")
library(metacoder);packageVersion("metacoder")
library(hilldiv);packageVersion("hilldiv")
library(car);packageVersion("car")
library(lme4);packageVersion("lme4")
library(RColorBrewer);packageVersion("RColorBrewer")
library("ape")
library(patchwork)
library(ggpubr)
library(rstatix)

#Load in data
meta <- read.csv("Sampling_meta.csv" , sep=",",stringsAsFactors = TRUE)
asv.tab <- read.delim("Curated_Table.txt",sep=",")
asv.taxa <- read.delim("Curated_Tax.csv",sep=",")

#order the data a bit
rownames(asv.taxa) <- rownames(asv.tab)
rownames(meta) <- meta$Sample_ID
meta <- meta %>% mutate(sample_orf_control=if_else(Sample_type %in% c("Shield blank", "PCR control"),"CONTROL_SAMPLE","SAMPLE"))
meta$Sample_type <- factor(meta$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water","Shield blank"))

#prep for phyloseq
asv.taxa.df <- as.data.frame(asv.taxa)
asv.taxa.df$ASV.ID <- row.names(asv.taxa.df)

#make phyloseq object
ps <- phyloseq(otu_table(asv.tab, taxa_are_rows=TRUE),
                  tax_table(as.matrix(asv.taxa)),
                  sample_data(meta)
               )

#plot_library_size for decontam
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
fullsize <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_orf_control)) + geom_point() 
#looks good, likely no contamination based on this

#start decontamination
sample_data(ps)$is.neg <- sample_data(ps)$sample_orf_control == "CONTROL_SAMPLE"
#Then exclude the positive control, we don´t want that to interfere with the process, since its not a blank and not a sample. 
#Just go with the default settings since we are not troubled with alot of contamination
contam.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg")
#remove the contaminants
ps.filt.prev05 <- prune_taxa(!contam.prev05$contaminant,ps)
#Some sanity check of where the contaminants are present based on abundance yay.
ps.pa <- transform_sample_counts(ps, function(x) x/sum(x))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_orf_control == "CONTROL_SAMPLE", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_orf_control == "SAMPLE", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
contaminant=contam.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#Write out some tables
write.table(ps.filt.prev05@otu_table,"ASV_counts_deconprev.tsv",sep="\t", quote=F, col.names=NA)
write.table(ps.filt.prev05@tax_table, "ASVs_taxonomy_deconprev.tsv",sep="\t", quote=F, col.names=NA)

#Now we´re done with decontamination some analyses can start 
#read in those tables again so we can estimate based on taxonomy "NA"´s what taxonomic level is appropriate to merge on
asv_tab <- read.delim("ASV_counts_deconprev.tsv",sep="\t")
taxa_tab <- read.delim("ASVs_taxonomy_deconprev.tsv",sep="\t")

#Write simple functions for showing the % of the various taxonomic ranks with "NA"´s assigned
prop.na <- function(tax_tab,Rank,asv_tab){
  filt <- tax_tab  %>% filter(is.na(tax_tab[,Rank])) %>% rownames(.)
  prop <- asv_tab %>% filter(rownames(.) %in% filt) 
  sum <- sum(rowSums(prop))/sum(rowSums(asv.tab))
  print(sum)
  paste(round(sum*100,digits = 3),"%")}
#make it write you a table so it shows 
table.prop.na <- function(Rank,asv_tab,taxa_tab){
  df <- data.frame()
  for(tax_rank in Rank){
    a<- as.data.frame(cbind(prop.na(taxa_tab,tax_rank,asv_tab),tax_rank,no_ASV=sum(is.na(taxa_tab[,tax_rank])))) %>% 
      mutate(percent_ASV_NA = V1) %>% 
      select(-V1)
    df <-  rbind(a,df)
  }
  print(df)
}
#run those
rank <- c(colnames(tax_tab))
table.prop.na(rank,asv_tab,taxa_tab)
#We decide to merge on the genus level
phylobby <- ps.filt.prev05
phylobby <- tax_glom(phylobby,NArm = FALSE,taxrank = "Genus")
random_tree = rtree(ntaxa(phylobby), rooted=TRUE, tip.label=taxa_names(phylobby))
phylobby = merge_phyloseq(phylobby, random_tree)
#Generate sample based relative abundace (and filter low occuring ASVs and samples with low counts)
physeq_filt = prune_samples(sample_sums(phylobby) >= 2500, phylobby)
physeq_filt = subset_samples(physeq_filt,Sample_type!="Shield blank")
# relative abundance
physeq_norm = transform_sample_counts(physeq_filt, function(x) x/sum(x))
physeq_norm = phyloseq::filter_taxa(physeq_norm, function(x) var(x) > 0, TRUE)                                 
#Finally we hav the phyloseq object we want to use we can make some basic overview plots
TopNOTUs = names(sort(taxa_sums(physeq_norm),TRUE)[1:30])
physeq_30= prune_taxa(TopNOTUs, physeq_norm) 
#make barplot
bar_genus <- plot_bar(physeq_30,fill="Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(colour = "black",face = "bold"),
        strip.background = element_blank()) +
  labs(x = "Samples") +
  labs(y = "Relative Abundance") + 
  scale_fill_manual(values=harmony_joy_palette) +
  scale_color_manual(values=harmony_joy_palette) +
  ggforce::facet_row(~Sample_type,scales = "free_x",space = "free")

#Print plot                                    
#pdf("barplot_16S_genus.pdf", width = 15, height = 10)
#bar_genus
#dev.off()                                  
                                    
#Run gUnifrac ordinations and PERMANOVAS
a <- physeq_norm %>% 
  dist_calc("gunifrac",gunifrac_alpha = 0) %>%  
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 2), color = "Sample_type", size = 5) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size = 15)) + 
  labs(subtitle = "GUnifrac (alpha=0)") +
  scale_color_manual(values=childish)

#gunifrac 0.5
b <- physeq_norm %>% 
  dist_calc("gunifrac",gunifrac_alpha = 0.5) %>%  
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 2), color = "Sample_type", size = 5) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size = 15)) + 
  labs(subtitle = "GUnifrac (alpha=0.5)") +
  scale_color_manual(values=childish) 
#gunifrac 1
c <- physeq_norm %>% 
  dist_calc("gunifrac",gunifrac_alpha = 1) %>%  
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 2), color = "Sample_type", size = 5) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size = 15)) + 
  labs(subtitle = "GUnifrac (alpha=1)") +
  scale_color_manual(values=childish)

d <- physeq_norm %>% 
  dist_calc("bray") %>%  
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 2), color = "Sample_type", size = 5) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size = 15)) + 
  labs(subtitle = "GUnifrac (alpha=1)") +
  scale_color_manual(values=childish)
#permanovas
physeq_norm %>% 
  dist_calc("gunifrac",gunifrac_alpha = 0) %>%
  dist_permanova(variables = "Sample_type", n_perms = 99, seed = 123) %>%
  perm_get()
physeq_norm %>% 
  dist_calc("gunifrac",gunifrac_alpha = 0.5) %>%
  dist_permanova(variables = "Sample_type", n_perms = 99, seed = 123) %>%
  perm_get()
physeq_norm %>% 
  dist_calc("gunifrac",gunifrac_alpha = 1) %>%
  dist_permanova(variables = "Sample_type", n_perms = 99, seed = 123) %>%
  perm_get()
#print plots
#pdf("ordinations.pdf", width = 25, height = 10)
#a + b +c + plot_layout(guides="collect")
#dev.off()

#print bray ordination
#pdf("ordination_bray.pdf", width = 10, height = 7)
#d 
#dev.off()

#plot alpha diversities                                    
phylobby_filt <- subset_samples(phylobby,Sample_type!="Shield blank")

ASV_al <-  as.data.frame(phylobby_filt@otu_table)


sample_groups <- as.data.frame(as.matrix(phylobby_filt@sam_data))
sample_groups <- as.data.frame(sample_groups[,c(1,2)])

f_richness <- div_test(ASV_al,hierarchy = sample_groups,q=0,posthoc = TRUE)
f_Shannon <- div_test(ASV_al,hierarchy = sample_groups,q=1,posthoc = TRUE)
f_Simpson <- div_test(ASV_al,hierarchy = sample_groups,q=2,posthoc = TRUE)

div_test_plot(f_richness,posthoc=TRUE,threshold=0.05,col=childish)
div_test_plot(f_Shannon,posthoc=TRUE,threshold=0.05,col=childish)
div_test_plot(f_Simpson,col=childish)
