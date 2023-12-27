#Load libraries
library(readxl);packageVersion("readxl")
library(patchwork)
library(tidyverse);packageVersion("tidyverse")
library(ggplot2);packageVersion("ggplot2")
library(ggpubr);packageVersion("ggpubr")
library(ggfortify);packageVersion("ggfortify")
library(rstatix);packageVersion("rstatix")
library(phyloseq);packageVersion("phyloseq")
library(ape);packageVersion("ape")
library(metacoder);packageVersion("metacoder")
library(cowplot);packageVersion("cowplot")
library(boral);packageVersion("boral")
library(hilldiv);packageVersion("hilldiv")
library(vegan);packageVersion("vegan")
library(ADImpute)
library(MicEco)

#set some nice palettes
#load the lovely palettes
#old ugly one
ugly<- c("#C15444","#C1B244","#26736A","#C17A44","#2D8652","#2E428A","#A83861","#BEE0E9","#D4C57D","#4D3619","#883091","#B3B4E6","#CDCF6E","#6ECFA8","#1F1F1F","#5E5E5E","#0E220B","#CD936A","#62CB8C","#C65353","#C6BA53","#4FB3C4","#8FC757","#6B4424","#EC453C","#F4A590","#D6F7AB","#78F297", "#310E8B","#EA462A")


#old nice one
childish <- c("#B71C1C","#FFC107","#A7FFEB","#FFB74D","#757575","#FF5252","#BA68C8","#8BC34A","#F57F17","#009688","#FF80AB","#304FFE","#455A64","#69F0AE","#FF5722","#9C640C","#F1C40F","#A93226","#76448A","#17A589","#C6BA53","#4FB3C4","#8FC757","#6B4424","#EC453C","#F4A590","#D6F7AB","#78F297", "#310E8B","#EA462A")


#same colors as in anvio
anvi_palette <- c("#E8BC0E","#F04343","#B357DE","#C77F12","#0E12E8")

#chatgpt generated palette
chat <- c("#FF0000", "#FF1700", "#FF2F00", "#FF4700", "#FF5E00", "#FF7600", "#FF8D00", "#FFA500", "#FFBC00", "#FFD400",
          "#FFEB00", "#F4FF00", "#D7FF00", "#BBFF00", "#9EFF00", "#81FF00", "#65FF00", "#48FF00", "#2BFF00", "#0FFF00",
          "#00FF19", "#00FF36", "#00FF52", "#00FF6F", "#00FF8B", "#00FFA8", "#00FFC4", "#00FFE1", "#00FFFF", "#00E1FF",
          "#00C4FF", "#00A8FF", "#008BFF", "#006FFF", "#0052FF", "#0036FF", "#0019FF", "#0F00FF", "#2B00FF", "#4800FF",
          "#6500FF", "#8100FF")
#Lets go and extract some basic goods from the data
meta <- read.csv("meta_sheet_final.csv")
#filter out the samples that were not sequenced or failed sequencing
meta <- meta %>% 
  filter(!Sample_type%in%c("Shield blank","gut")) %>% 
  filter(Sample_ID!="W_002") 
#re-order for the rest of the analyses and ordering on plots
meta$Sample_type <- factor(meta$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))

#average DNA concentration
av_conc <- meta %>% 
  group_by(Sample_type) %>% summarize(mean=mean(Qubit.novogene),
                                      sd=sd(Qubit.novogene),
                                      mean_cov=mean(Depth_zeb_map))
#extraction of basic overivew numbers

#total raw reads
nr_raw <- sum(meta$raw_reads)
#total unmapped reads
nr_unmapped <- sum(meta$unmapped)
#total_quality_pass
nr_quality_pass <- sum(meta$rmad_rmdup_repair)

av_numbs <- meta %>% 
  group_by(Sample_type) %>% 
  summarize(mean_qual=mean(Qubit.novogene),
            mean_host_fil=mean(raw_reads),
            mean_perc_qual=mean(rmad_rmdup_repair),
            mean_perc_host_qual=mean(unmapped))

#Make function for plotting the basic stats
plot_stuff <- function(col_name,y_name,just){
  meta %>% 
    group_by(Sample_type) %>% 
    summarize(mean=mean({{col_name}}),sd=sd({{col_name}})) %>%     
    ggplot(aes(x=Sample_type,y=mean,fill=Sample_type)) +
    geom_bar(stat="identity",alpha=0.8) + 
    geom_text(aes(label=round(mean,2)),hjust=just,vjust=1.5,check_overlap = FALSE) +
    scale_fill_manual(values=anvi_palette) + 
    geom_errorbar(aes(ymin=mean-sd,
                      ymax=mean+sd),linewidth=0.4) + 
    labs(y=y_name,fill="Sample type") +
    theme_cowplot()} 

#make them basic stat plots
a <- plot_stuff(Qubit.novogene,"DNA concentration(ng/ul)",-1) 
b <- plot_stuff(perc_raw_repair,"% reads removed in Quality filtering",3)
c <- plot_stuff(perc_clean_unmapped,"% host reads",-0.5)
d <- plot_stuff(Depth_zeb_map,"Zebrafish depth of coverage(X)")

# add seome basic assembly stats
assembly_stat <- read.delim("assembly_stat.txt")
assembly_stat$Sample_type <- factor(assembly_stat$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))
#plot some assembly stats to add to the plot
e <- ggplot(assembly_stat,aes(x=Sample_type,fill=Sample_type)) + geom_bar(aes(y=nr_contigs/1000),stat="identity") +   scale_fill_manual(values=anvi_palette) + labs(y= "Number contigs(Kb)") +
  theme_cowplot()
 
f <- ggplot(assembly_stat,aes(x=Sample_type,fill=Sample_type)) + geom_bar(aes(y=N50/1000),stat="identity",position = "dodge") + labs(y="N50(Kb)") +
   scale_fill_manual(values=anvi_palette)+
  theme_cowplot()

g <- ggplot(assembly_stat,aes(x=Sample_type,fill=Sample_type)) +
geom_bar(aes(y=GC_content),stat="identity") + labs(y="% GC-content") +
  scale_fill_manual(values=anvi_palette)+
  theme_cowplot()

#make a lovely large plot
#pdf("sample_processing.pdf",width=20,height=7)
#a +b +c + plot_spacer() + e +f +g + plot_layout(nrow = 1)  & 
#coord_flip() & 
#scale_x_discrete(limits=rev(levels(meta$Sample_type))) & 
#theme(axis.text.y = element_blank(),
#axis.title.y = element_blank(),
#       legend.position = "none", 
#        axis.ticks.y = element_blank(),
#       axis.title.x = element_text(size=12, face="bold")) 
#Basic stats done
#Now compare quality of MAGs among sample types
MAG_qual_stat <- read.csv("MAG_qual_stat.csv")
#set order
MAG_qual_stat$Sample_type <- factor(MAG_qual_stat$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))

#tidy time
plot_stats <- pivot_longer(MAG_qual_stat,3:7,values_to = "Value",names_to = "type")
plot_stats <- plot_stats %>% mutate(two=ifelse(type %in% c("nr_actual","nr_expected"),"Actual vs expected nr MAGs","high quality vs lower quality MAGS")) %>% filter(type!="perc_exp")
#g

nr_mag_compare <- ggplot(plot_stats) + geom_point(aes(x=Sample_type,y=Value,color=Sample_type,shape=type),size=7,alpha=0.8) + ggforce::facet_row(~two,scales="free_x",space="free") + 
theme_cowplot() + 
theme(strip.background = element_blank(),
      strip.text = element_text(size=12,face="bold")) +
scale_color_manual(values = anvi_palette) + 
scale_shape_manual(name="",labels=c("Nr actual MAGs","Nr expected MAGs","Nr high quality MAGs","Nr lower quality MAGs"),values = c(19,15,17,18)) +  
  scale_y_continuous(breaks=scales::pretty_breaks(n=12)) + 
  labs(y="Number of MAGs")
#Make rarefaction curves of ORFs among the samples
#Read in the ORF coverage data and tidy it a bit
gene_cov <- read.csv("MAG_REF-GENE-COVERAGES.txt",sep="")
rownames(gene_cov) <- gene_cov$key
gene_cov <- gene_cov %>% select(-key)
gene_cov <- as.matrix(gene_cov)
class(gene_cov) <- "integer"
#Plot rarefaction curves
rarecurve(t(gene_cov), step = 50000, sample = gene.raremax, col = "black", cex = 0.6, label = FALSE, xlab = "Contig Coverage", ylab = "Gene calls (ORFs)")
#Lovely looks like we expected now lets get the data into a tidy frame so we can plot
rare_curve <- rarecurve(t(gene_cov), step = 50000, col = "black", cex = 0.6, label = FALSE, xlab = "Contig Coverage", ylab = "Gene calls (ORFs)",tidy = TRUE)
#Now we add a grouping variable to the rare_curve data frame based on the sample and plot
rare_mod <- rare_curve %>% 
  mutate(sample_type=str_split(Site,"_",n=2,simplify = TRUE)[,1]) %>% 
  mutate(sample_type= ifelse(Site=="merged_SZ_002","SZ",sample_type)) %>%
  mutate(SAMPLE_TYPE=ifelse(sample_type=="W","Water",ifelse(sample_type=="SZ","Squeezed gut",ifelse(sample_type=="WG","Whole gut",ifelse(sample_type=="DG","Intestinal content","Feces"))))) %>%
  mutate(SAMPLE_TYPE=factor(SAMPLE_TYPE,levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))) %>%
  ggplot() + geom_line(aes(x=Sample,y=Species,group=Site,color=SAMPLE_TYPE),size=2) 

#plot abd add some aesthetics to it
#and print if wanted
#pdf("ORF_rarefy.pdf",width=9,height=6)
rare_mod + 
  theme_cowplot()+
  scale_color_manual(values=anvi_palette,name="Sample Type") +
  labs(x="Read depth",y="Gene calls(ORFs)") +
  theme(axis.line = element_line(size=1,color="black"),
        axis.title = element_text(size = 16,face="bold"),
        axis.text = element_text(size = 12,face="bold"),
        legend.position = c(0.80,0.20),
        legend.text = element_text(face="bold",size=14),
        legend.title = element_text(face="bold",size=16))
#dev.off()
#rarefaction analyses done
#Now its time to normalise the MAG abundance based on MAG coverage and MAG size
#Load in MAG coverage data and summary
ag_cov <- read.csv("mean_coverage.txt",sep="")
rownames(mag_cov) <- mag_cov$bins
mag_cov <- mag_cov %>% select(-bins)
mag_sum <- read.delim("SUMMARY/bins_summary.txt")
#make data frame for the size
size <- mag_sum %>% select(bins,total_length) %>% rename(hgnc_symbol=bins,transcript_length=total_length)
#tpm_normalisation
mag.tmp <- NormalizeTPM(mag_cov,tr_length = size,scale = 1e+06)
#Normalising done
#Create tax table
tax_table <- mag_sum[,c(1,8:14)]
rownames(tax_table) <- tax_table$bins 
taxa_table <- tax_table %>% 
  select(-bins) %>% mutate_all(na_if,"")
#read into phyloseq the normalised mag abundances
rownames(meta) <- meta$exp_ID 
meta <- meta %>% select(-exp_ID)
#Make phyloseq object
py_mag <- phyloseq(otu_table(as.matrix(mag.tmp),taxa_are_rows=TRUE),
                   tax_table(as.matrix(taxa_table)),
                   sample_data(meta))
#Lovely
#Make a quick barplot to get the genus level overview
#This is however much better viewed using anvi-interactive as per figure 3,
#but we just do this to get an overview and an idea of whats up and to see the most abundant MAGs among the sample types
barplot_genus <- plot_bar(py_mag,"Sample_ID",fill="t_genus") + 
  geom_col(aes(color=t_genus),position = "stack",show.legend = FALSE) + 
  scale_fill_manual(values=childish,name="Genus") + 
  scale_color_manual(values=childish)+ 
  ggforce::facet_row(~Sample_type, scales="free_x",space = "free") + 
  guides(fill=guide_legend(ncol = 1)) + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90,face = "bold"),
        axis.text.y = element_text(face="bold",size = 12),
        axis.title = element_text(face="bold",size = 12),
        strip.text = element_text(face="bold",size=12),
        strip.background.x = element_blank(),
        legend.title = element_text(face="bold",size = 12),
        legend.text = element_text(face="bold",size=12),
        panel.border = element_blank())
Now have a look a the top taxa across the dataset and how their abundances are distributed
N <- 30
barplot(sort(taxa_sums(py_mag), TRUE)[1:N]/nsamples(py_mag), las=2,col="#9C640C")+ title(main="Cumulative fractional abundance of top 30 in thesamples",ylab = "Fractional abundance")

#the fololowing the MAGs that make up the most of the dataset
#Here are the ones which make up more than 5 % 
"WG_ZS_MAG_00001" "SZ_ZS_MAG_00001" "DG_ZS_MAG_00008"
#here are the ones which make up more than 1%
Top10 <- c("WG_ZS_MAG_00001",
           "SZ_ZS_MAG_00001",
           "DG_ZS_MAG_00008",
           "PW_ZS_MAG_00004", 
           "PW_ZS_MAG_00008",
           "DG_ZS_MAG_00003",
           "SZ_ZS_MAG_00011", 
           "SZ_ZS_MAG_00009",
           "DG_ZS_MAG_00001",
           "W_ZS_MAG_00013")

#function to get an the individual abundances
mean_ab <- function(MAG) {
  ab <- abundances %>% 
    group_by(Sample_type) %>%
    filter(OTU==MAG) %>% 
    summarise(mean=mean(Abundance)*100)
  ab$mean <- paste(round(ab$mean,digits=3),"%")
  print(ab)
  }
#for loop to get all the mean abundances 
result_list <- data.frame()

for (name in Top10) {result <- mean_ab(name)
print(result)
result_list <- c(result_list,result,name)}
merged_result <- do.call(rbind,result_list)
print(merged_result)
#yay now lets make an oridation plot
GP.ord <- ordinate(py_mag,"PCoA","bray")
#plot and color by sample types (we also make a plot where we color a gradient of host filtered reads
ordi.mag <- plot_ordination(py_mag,GP.ord,type="samples",color="Sample_type") + 
  geom_point(size=7) +
   theme_cowplot() +
  scale_color_manual(values=anvi_palette) +
  theme(axis.text.x = element_text(angle = 90,face = "bold"),
        axis.text.y = element_text(face="bold",size = 12),
        axis.title = element_text(face="bold",size = 12),
        strip.text = element_text(face="bold",size=12),
        strip.background.x = element_blank(),
        legend.title = element_text(face="bold",size = 12),
        legend.text = element_text(face="bold",size=12),
        panel.border = element_blank(),
        axis.line = element_line(size=1,color="black"))
#PERMANOVA for sample types
adonis2(t(as.data.frame(py_mag@otu_table))~Sample_type,data=data.frame(as.vector(py_mag@sam_data)),method="bray")

#Now lets compare the MAG abundances to the 16S abundances on the genus level
#If you follow the 16S code you will see an object called "physeq_norm"
#Set that to the object physeq_16 so we can make a little comparison
physeq_16S <- physeq_norm
#merge taxa on the genus level in the mag dataset
py_mag_genus <- tax_glom(py_mag,taxrank = "t_genus",NArm = FALSE)
#get the top 10 genera from the MAG dataset
TopNOTUs = names(sort(taxa_sums(py_mag_genus),TRUE)[1:10])
physeq_10_mag= prune_taxa(TopNOTUs, py_mag_genus)
#get the top 10 genera from the 16S dataset
TopNOTUs_16S = names(sort(taxa_sums(physeq_16S),TRUE)[1:10])
physeq_10_16S= prune_taxa(TopNOTUs_16S, physeq_16S) 
#remove the two samples from the 16S dataset that failed shotgun sequencing
physeq_10_16S <- physeq_10_16S %>% subset_samples(!Sample_ID %in% c("W_002","SB_002"))
top10_mag <- data.frame(physeq_10_mag@tax_table)
top10_16S <- data.frame(physeq_10_16S@tax_table)
common_mag <- top10_mag %>% filter(t_genus %in% top10_16S$Genus)
common_16S <- top10_16S %>% filter(Genus %in% top10_mag$t_genus)
filt_mag <- rownames(common_mag)
filt_16S <- rownames(common_16S)  
physeq_common_mag <- prune_taxa(filt_mag,physeq_10_mag)
physeq_common_16S <- prune_taxa(filt_16S,physeq_10_16S)
common_mag_ab <- data.frame(t(data.frame(physeq_common_mag@otu_table)))
common_16S_ab <- data.frame(t(data.frame(physeq_common_16S@otu_table)))
common_16S_ab$Sample_ID <- rownames(common_16S_ab)
#get them all into a metadata thingy
meta_new <- data.frame(py_mag_genus@sam_data)
together <- cbind(common_mag_ab,meta_new) 
together <- merge(together,common_16S_ab,by="Sample_ID")
#Now sort the dataframe in the way that we want and first define the which belong to the same genus
Cetobacterium <- c("ASV_1","WG_ZS_MAG_00001")
Aeromonas <- c("ASV_3","SZ_ZS_MAG_00001")
Pseudomonas <- c("ASV_32","PW_ZS_MAG_00020")
Ideonella <- c("ASV_7","PW_ZS_MAG_00004")
#plot and compare the abundances
compare_ab <- together %>% 
  select(c(filt_16S,filt_mag,Sample_type,Sample_ID)) %>%
  pivot_longer(cols = 1:10,values_to = "abundances",names_to = "taxa_ID") %>% 
  mutate(origin = ifelse(taxa_ID %in% filt_16S,"16S","Shotgun")) %>%
  mutate(genus=ifelse(taxa_ID %in% Cetobacterium,"Cetobacterium",
                      ifelse(taxa_ID %in% Aeromonas,"Aeromonas",
                             ifelse(taxa_ID %in% Ideonella, "Ideonella",
                                    ifelse(taxa_ID %in% Pseudomonas,"Pseudomonas","ZOR0006"))))) %>% 
  group_by(genus,origin,Sample_type) %>%
  filter(genus!="Ideonella") %>%
  summarise(mean=mean(abundances)) %>% 
  pivot_wider(names_from = "origin",values_from = "mean") %>%
  ungroup() %>%
  ggplot(aes(x=log2(`16S`),y=log2(Shotgun),color=genus,group=1)) + 
  geom_point(size=7) + 
  #geom_smooth(se=FALSE,method = "lm",linetype="dashed") +
  geom_abline(slope=1,intercept=0) +
  #coord_equal() +
  #xlim(-15,0) +
  #ylim(-15,0) +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = childish) +
    theme_cowplot() +
  ggforce::facet_row(~Sample_type,space="free")  +
  theme(axis.text.x = element_text(angle = 90,face = "bold"),
        axis.text.y = element_text(face="bold",size = 12),
        strip.text = element_text(face="bold",size=12),
        axis.title = element_text(face="bold",size=12),
        strip.background.x = element_blank(),
        legend.title = element_text(face="bold",size = 12),
        legend.text = element_text(face="bold",size=12),
        panel.border = element_blank()) +
  labs(x="log2 relavtive abundance  16S",y="log2 relavtive abundance shotgun")
#Awesome a great relationship. We chose log2 transformation of the relative abundances to give the lower abundance taxa
#a chance
#Ok time to compare the functional data among sample types
MAG_REF.GENE.CALLS <- read.delim("MAG_REF-GENE-CALLS.txt")
MAG_REF.GENE.COVERAGES <- read.delim("MAG_REF-GENE-COVERAGES.txt")
gene_coverages <- MAG_REF.GENE.COVERAGES
#read in Pfam functions
Pfam_functions <-read.delim("Pfam_functions.txt")
#read in COG data
COG_CAT <- read.delim("COG_CAT.txt")
#read in the mapping data for normalisation of gene coverage data
default <- read.delim("default.txt")

#lets have a look at the COG categries
#Keep only the genes that have an assigned COG_CATEGORY
gene_coverages <- gene_coverages %>% filter(key %in% COG_CAT$gene_callers_id)
#clear up for phyloseq and normalisation
rownames(gene_coverages) <- gene_coverages$key
gene_coverages <- gene_coverages %>% select(-key)
#gene cov ready
#now get the COG data ready
rownames(COG_CAT) <- COG_CAT$gene_callers_id
COG_CAT <- COG_CAT %>% select(-gene_callers_id)
COG_CAT$FUNCTION <- COG_CAT$function.
COG_CAT <- COG_CAT %>% select(-function.)
#Here we keep only the first function assigned to the gene
#This seems to be a 'user' choice as far as I can ascertain... 
COG_CAT <- COG_CAT %>% mutate(FUNCTION=str_split(FUNCTION,"!!!",simplify = TRUE)[,1],
                   accession=str_split(accession,"!!!",simplify=TRUE)[,1])
#COG data ready
#normalise gene coverages
#make sure that the deviding is correct
rep_total <- rep(default$total_reads_mapped,each=nrow(gene_coverages))
gene_cov_norm <- gene_coverages/rep_total
gene_cov_norm <- gene_cov_norm * 1000000
#Make phyloseq object
phy_cog <- phyloseq(otu_table(as.matrix(gene_cov_norm),taxa_are_rows=TRUE),
                   tax_table(as.matrix(COG_CAT)),
                   sample_data(meta))
#Make venn diagram
venn_cog <- ps_venn(phy_cog,group="Sample_type",fill=anvi_palette,alpha=0.8)
#make venn diagram a data_frame to extract and plot further
lol <- ps_venn(phy_cog,group="Sample_type",plot=FALSE)
#Focus in on the three different fractions, just for fun create separe data.frames for those
shared_all <- data.frame(lol$`Whole gut__Intestinal content__Squeezed gut__Feces__Water`)
shared_wo_wg <- data.frame(lol$`Intestinal content__Squeezed gut__Feces__Water`) 
pw_w_share <- data.frame(lol$`Feces__Water`)
#focus in on the different fractions create seperate phyloseq objects for them
pw_w_func <- prune_taxa(pw_w_share$lol.Feces__Water,phy_cog)
wo_wg_func <- prune_taxa(shared_wo_wg$lol..Intestinal.content__Squeezed.gut__Feces__Water.,phy_cog)
shared_all_func <- prune_taxa(shared_all$lol..Whole.gut__Intestinal.content__Squeezed.gut__Feces__Water.,phy_cog)
#shared poop water fraction
w_frac <- sample_sums(pw_w_func)/sample_sums(phy_cog)*100
dw_frac <- data.frame(w_frac)
#shared among all except whole gut samples
wo_wg_frac <- sample_sums(wo_wg_func)/sample_sums(phy_cog)*100
dwo_wg_frac <- data.frame(wo_wg_frac)
#core fraction shared among all samples
shared_frac <- sample_sums(shared_all_func)/sample_sums(phy_cog)*100
dshared_frac <- data.frame(shared_frac)
#get sample_types from meta for plotting
meta_comp_frac <- meta %>% select(Sample_type) 
meta_comp_frac$names <- rownames(meta_comp_frac)
meta_comp_frac <- meta_comp_frac %>% arrange(names)
#make new data frame for plotting the three fractions and make tidy
comp_frac <- cbind(dw_frac,dwo_wg_frac,dshared_frac,meta_comp_frac)
#lets first extract the means 
comp_frac <- pivot_longer(comp_frac,1:3,names_to="type",values_to = "frac")
comp_frac$Sample_type <- factor(comp_frac$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))
#lets first extract the means 
comp_frac %>% group_by(Sample_type,type) %>% summarise(mean=mean(frac),sd=sd(frac))
#Then make plot
frac_cog <- ggplot(comp_frac) +
  geom_point(aes(x=type,y=frac/100,color=Sample_type),size=5,position = position_dodge(width = 2),alpha=0.5) +
  theme_cowplot() +
  theme(strip.text = element_blank()) +
  ggforce::facet_row(~type,scales = "free",space="free") + 
  scale_color_manual(values = anvi_palette) +
  scale_y_continuous(labels = scales::percent,breaks = scales::pretty_breaks(n=8)) +
  labs(y="% of total normailised gene coverage",x="")
#happy
#print to pdf if wished
#pdf("frac_cog.pdf",width=10,height=4)
frac_cog
#dev.off()
#Lets see if removing functions that are higher in coverage in the water than the fecal samples
#does something

#get the mean coverage of genes in Water and Feces
mean_abundance_water <- rowMeans(otu_table(phy_cog)[, sample_data(phy_cog)$Sample_type == "Water"])
mean_abundance_feces <- rowMeans(otu_table(phy_cog)[, sample_data(phy_cog)$Sample_type == "Feces"])
#make list of gene calls to keep in the phyloseq object
#that is genes with higher coverage in the fecal samples compared ot the water samples
taxa_to_remove <- names(which(mean_abundance_water < mean_abundance_feces))
#make phyloseq object with gene calls we want to keep
phy_cog_filt <- prune_taxa(taxa_to_remove,phy_cog)
#make_venn_diagram of these functions
venn_cog_post <- ps_venn(phy_cog_filt,group="Sample_type",fill=anvi_palette,alpha=0.8)
#make the venn diagrem a non plot object so we can see the new shared functions
ps_venn_hcw_fil <- ps_venn(phy_cog_filt,group="Sample_type",plot = FALSE)

#now repeat the same analysis of the different fractions after removing all the shit.
shared_all_post <- data.frame(ps_venn_hcw_fil$`Whole gut__Intestinal content__Squeezed gut__Feces__Water`)
shared_wo_wg_post <- data.frame(ps_venn_hcw_fil$`Intestinal content__Squeezed gut__Feces__Water`) 
pw_w_share_post <- data.frame(ps_venn_hcw_fil$`Feces__Water`)

pw_w_func_post <- prune_taxa(pw_w_share_post$ps_venn_hcw_fil.Feces__Water,phy_cog_filt)
wo_wg_func_post <- prune_taxa(shared_wo_wg_post$ps_venn_hcw_fil..Intestinal.content__Squeezed.gut__Feces__Water.,phy_cog_filt)
shared_all_func_post <- prune_taxa(shared_all_post$ps_venn_hcw_fil..Whole.gut__Intestinal.content__Squeezed.gut__Feces__Water.,phy_cog_filt)

#shared poop water fraction
w_frac_post <- sample_sums(pw_w_func_post)/sample_sums(phy_cog)*100
dw_frac_post <- data.frame(w_frac_post)
#shared among all except whole gut samples
wo_wg_frac_post <- sample_sums(wo_wg_func_post)/sample_sums(phy_cog)*100
dwo_wg_frac_post <- data.frame(wo_wg_frac_post)
#core fraction shared among all samples
shared_frac_post <- sample_sums(shared_all_func_post)/sample_sums(phy_cog)*100
dshared_frac_post <- data.frame(shared_frac_post)

#no need to do anything to the metadata its already up there and called meta_comp_frac

#make new data frame for plotting the three fractions and make tidy
comp_frac <- cbind(dw_frac_post,dwo_wg_frac_post,dshared_frac_post,meta_comp_frac)
#lets first extract the means 
comp_frac <- pivot_longer(comp_frac,1:3,names_to="type",values_to = "frac")
comp_frac$Sample_type <- factor(comp_frac$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))
#lets first extract the means 
comp_frac %>% group_by(Sample_type,type) %>% summarise(mean=mean(frac),sd=sd(frac))

#Then make plot
frac_cog_post <- ggplot(comp_frac) +
  geom_point(aes(x=type,y=frac/100,color=Sample_type),size=5,position = position_dodge(width = 2),alpha=0.5) +
  theme_cowplot() +
  theme(strip.text = element_blank()) +
  ggforce::facet_row(~type,scales = "free",space="free") + 
  scale_color_manual(values = anvi_palette) +
  scale_y_continuous(labels = scales::percent,breaks = scales::pretty_breaks(n=8)) +
  labs(y="% of total normailised gene coverage",x="")
###Now we make a venn diagram from the gene calls that are present only in the Top10_mags
#first read in all 10 gene_calls files for the MAGS
file_names <- c("WG_ZS_MAG_00001-gene_calls.txt", "SZ_ZS_MAG_00001-gene_calls.txt", "DG_ZS_MAG_00008-gene_calls.txt", 
                "PW_ZS_MAG_00004-gene_calls.txt", "PW_ZS_MAG_00008-gene_calls.txt", "DG_ZS_MAG_00003-gene_calls.txt", 
                "SZ_ZS_MAG_00011-gene_calls.txt", "SZ_ZS_MAG_00009-gene_calls.txt", "DG_ZS_MAG_00001-gene_calls.txt", "W_ZS_MAG_00013-gene_calls.txt")
object_names <- c("WG_ZS_MAG_00001.gene_calls", "SZ_ZS_MAG_00001.gene_calls", "DG_ZS_MAG_00008.gene_calls", 
                  "PW_ZS_MAG_00004.gene_calls", "PW_ZS_MAG_00008.gene_calls", "DG_ZS_MAG_00003.gene_calls", 
                  "SZ_ZS_MAG_00011.gene_calls", "SZ_ZS_MAG_00009.gene_calls", "DG_ZS_MAG_00001.gene_calls","W_ZS_MAG_00013.gene_calls")
data_list <- list()
for (i in seq_along(file_names)) {
  directory_name <- gsub("-.*", "", file_names[i])
  full_path <- file.path("~/Documents/PhD/experiments/sampling_trial_for_host-vs-bacterial/Anvio/DREP_ALL/SUMMARY/bin_by_bin", directory_name, file_names[i])
  data_list[[object_names[i]]] <- read.delim(full_path, header = TRUE)
}

#get the files in into the environment
list2env(data_list, envir = .GlobalEnv)

#lets merge them into a single data.frame of tons of functions
top_mag_gene_calls <- 
  rbind(WG_ZS_MAG_00001.gene_calls,SZ_ZS_MAG_00001.gene_calls,DG_ZS_MAG_00008.gene_calls, PW_ZS_MAG_00004.gene_calls, PW_ZS_MAG_00008.gene_calls, DG_ZS_MAG_00003.gene_calls, SZ_ZS_MAG_00011.gene_calls, SZ_ZS_MAG_00009.gene_calls, DG_ZS_MAG_00001.gene_calls,W_ZS_MAG_00013.gene_calls)
#read back in the gene calls data frame for the entire dataset
COG_CAT <- read.delim("COG_CAT.txt")
COG_CAT_merge <- COG_CAT %>% filter(gene_callers_id %in% top_mag_gene_calls$gene_callers_id)
cov_gene_top10 <- MAG_REF.GENE.COVERAGES %>% filter(key %in% COG_CAT_merge$gene_callers_id)
#clear up for phyloseq and normalisation
rownames(cov_gene_top10) <- cov_gene_top10$key
cov_gene_top10 <- cov_gene_top10 %>% select(-key)
rownames(COG_CAT_merge) <- COG_CAT_merge$gene_callers_id
COG_CAT_merge <- COG_CAT_merge %>% select(-gene_callers_id)
COG_CAT_merge$FUNCTION <-COG_CAT_merge$function.
COG_CAT_merge <- COG_CAT_merge %>% select(-function.)
#Here we keep only the first function assigned to the gene
COG_CAT_merge <- COG_CAT_merge %>% mutate(FUNCTION=str_split(FUNCTION,"!!!",simplify = TRUE)[,1],
                   accession=str_split(accession,"!!!",simplify=TRUE)[,1])
#COG data ready
#normalise gene coverages
#make sure that the deviding is correct
rep_total <- rep(default$total_reads_mapped,each=nrow(cov_gene_top10))

gene_cov_norm_top10 <- cov_gene_top10/rep_total
gene_cov_norm_top10 <- gene_cov_norm_top10 * 1000000

phy_cog_top10 <- phyloseq(otu_table(as.matrix(gene_cov_norm_top10),taxa_are_rows=TRUE),
                   tax_table(as.matrix(COG_CAT_merge)),
                   sample_data(meta))

top10_venn <- ps_venn(phy_cog_top10,group="Sample_type",fill=anvi_palette)
top10_noplot <- ps_venn(phy_cog_top10,group="Sample_type",fill=anvi_palette,plot = FALSE)

#Lets have a look at the functional proportions of these guys from the total

shared_all_top10 <- data.frame(top10_noplot$`Whole gut__Intestinal content__Squeezed gut__Feces__Water`)
shared_wo_wg_top10 <- data.frame(top10_noplot$`Intestinal content__Squeezed gut__Feces__Water`) 
pw_w_share_top10 <- data.frame(top10_noplot$`Feces__Water`)

pw_w_func_top10 <- prune_taxa(pw_w_share_top10$top10_noplot.Feces__Water,phy_cog_top10)
wo_wg_func_top10 <- prune_taxa(shared_wo_wg_top10$top10_noplot..Intestinal.content__Squeezed.gut__Feces__Water.,phy_cog_top10)
shared_all_func_top10 <- prune_taxa(shared_all_top10$top10_noplot..Whole.gut__Intestinal.content__Squeezed.gut__Feces__Water.,phy_cog_top10)

#shared poop water fraction
w_frac_top10 <- sample_sums(pw_w_func_top10)/sample_sums(phy_cog)*100
dw_frac_top10 <- data.frame(w_frac_top10)
#shared among all except whole gut samples
wo_wg_frac_top10 <- sample_sums(wo_wg_func_top10)/sample_sums(phy_cog)*100
dwo_wg_frac_top10 <- data.frame(wo_wg_frac_top10)
#core fraction shared among all samples
shared_frac_top10 <- sample_sums(shared_all_func_top10)/sample_sums(phy_cog)*100
dshared_frac_top10 <- data.frame(shared_frac_top10)


comp_frac <- cbind(dw_frac_top10,dwo_wg_frac_top10,dshared_frac_top10,meta_comp_frac)
#lets first extract the means 
comp_frac <- pivot_longer(comp_frac,1:3,names_to="type",values_to = "frac")
comp_frac$Sample_type <- factor(comp_frac$Sample_type, levels=c("Whole gut", "Intestinal content", "Squeezed gut","Feces","Water"))
#lets first extract the means 
comp_frac %>% group_by(Sample_type,type) %>% summarise(mean=mean(frac),sd=sd(frac))

#Then make plot
frac_cog_top10 <- ggplot(comp_frac) +
  geom_point(aes(x=type,y=frac/100,color=Sample_type),size=5,position = position_dodge(width = 2),alpha=0.5) +
  theme_cowplot() +
  theme(strip.text = element_blank()) +
  ggforce::facet_row(~type,scales = "free",space="free") + 
  scale_color_manual(values = anvi_palette) +
  scale_y_continuous(labels = scales::percent,breaks = scales::pretty_breaks(n=8)) +
  labs(y="% of total normailised gene coverage",x="")

#lovely now we filter out the ideonella whihc is the one that is in higest abundance in the water samples
#And barely present in the others and see what effect that has. 



#find_table_to_prune away the ideonella
taxa_temove <- PW_ZS_MAG_00004.gene_calls$gene_callers_id
taxa_keep <- top_mag_gene_calls %>% filter(!gene_callers_id %in% taxa_temove) %>% select(gene_callers_id)
taxa_keep <- as.character(taxa_keep$gene_callers_id)
#Lets see how it looks for the top10
phy_cog_top10_wo <- prune_taxa(taxa_keep,phy_cog_top10)
top10_venn_ideo <- ps_venn(phy_cog_top10_wo,group="Sample_type",fill=anvi_palette)
#and then overall

taxa_keep <- COG_CAT %>% filter(!gene_callers_id %in% taxa_temove) %>% select(gene_callers_id)
taxa_keep <- as.character(taxa_keep$gene_callers_id)

phy_cog_full_ideo <- prune_taxa(taxa_keep,phy_cog)
full_venn_ideo <- ps_venn(phy_cog_full_ideo,group="Sample_type",fill=anvi_palette)

#I will write some functions for this whole mess at some point
#but anyways we do the same for Pfam functions as we have done with the above with slightly modified code 

#OK venn stuff done
#Now lets compare to EBI
#get total % mapped unmapped to this and ebi´s reference
mean_perc_map <- meta %>% 
  summarise(mean_drep=mean(perc_map_DREP),
                                    mean_ebi=mean(perc_map_EBI),
                                    mean_ori=mean(perc_ori))
per_type_mean_perc_map <- meta %>% 
  group_by(Sample_type) %>%  
  summarise(mean_drep=mean(perc_map_DREP),
            mean_ebi=mean(perc_map_EBI),
            sd_drep=sd(perc_map_DREP),
            sd_ebi=sd(perc_map_EBI)
            )
#Lovely lets put this into a lovely biorender figure

#finally lets make a correlation of number of variable positions reported per sample and correlate that to 
#the coverage of the Cetobacterium MAG.

read.delim("meta_sheet_final.txt")

#read in departure from consensus data
read.delim("view.txt")

#Wrangle a bit
var <- pivot_longer(view,cols=2:16) %>% 
  filter(value!=0) %>% 
  select(-contig) %>% 
  group_by(name) %>% 
  summarize(n()) %>%
  mutate(n_pos=`n()`) %>%
  select(-`n()`)
#make final frame
final_var <- meta_sheet_final %>% 
  arrange(samples) %>%  
  cbind(.,var)

#cor_test
cor.test(final_var$n_pos, final_var$Ceto_coverage)

#plot
ggplot(final_var,aes(x=n_pos,y=Ceto_coverage)) + geom_point() + stat_cor(method = "pearson") 
#yay
