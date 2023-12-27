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

#Per group average 

av_numbs <- meta %>% 
  group_by(Sample_type) %>% 
  summarize(mean_qual=mean(Qubit.novogene),
            mean_host_fil=mean(raw_reads),
            mean_perc_qual=mean(rmad_rmdup_repair),
            mean_perc_host_qual=mean(unmapped))


plot_stuff <- function(col_name,y_name,just){
  meta %>% 
    group_by(Sample_type) %>% 
    summarize(mean=mean({{col_name}}),sd=sd({{col_name}})) %>%     ggplot(aes(x=Sample_type,y=mean,fill=Sample_type)) +
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
assembly_stat <- read.delim("~/Documents/PhD/experiments/sampling_trial_for_host-vs-bacterial/Anvio/DREP_ALL/test_man/assembly_stat.txt")
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
a +b +c + plot_spacer() + e +f +g + plot_layout(nrow = 1)  & 
  coord_flip() & 
  scale_x_discrete(limits=rev(levels(meta$Sample_type))) & 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none", 
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12, face="bold")) 
