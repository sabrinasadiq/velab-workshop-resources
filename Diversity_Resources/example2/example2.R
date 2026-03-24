# Load packages
library(tibble)
library(doBy)
library(ggplot2)
library(ggpattern)
library(gtable)
library(grid)
library(gridExtra)
library(lattice)
library(dplyr)
library(Hmisc)
library(reshape2)
library(cowplot)
library(MASS)
library(viridis)
library(RColorBrewer)
library(treemap)
library(vegan)
library(phyloseq)
library(tidyr)
library(ape)
library(phytools)
library(phangorn)
library(multcomp)
library(gridExtra)
library(ggpubr)
library(ggh4x)
library(jtools)
library(pals)
library(car)
library(ggstance)
library(bestglm)
library(forcats)

# Set working directory
setwd("/path/to/file/Diversity_Resources/example2")
# Import data for quick glimpse figures
virome2 <-  read.csv("ViromeComp2.csv",na.strings=c("", "NA"), header=TRUE,sep=",")

#######################################################################################################################

rel2<-melt(virome2, id.vars=c("Library", "Temperature", "Group", "Time", "Collected", 
                              "Richness", "Shannon", "Shannon.effective", "TotalReads", "SumReads"),
           variable.name="Family", value.name="Count")

rel2$Count[is.na(rel2$Count)] <- 0
rel2$Pct_rel2<-(rel2$Count/rel2$SumReads)*100

# Reordering for facets later on
rel2$Temperature <- factor(rel2$Temperature, levels = c("2-8°C", "-30°C", "-80°C"))

rel2$FacetGroup <- factor(
  paste(rel2$Collected, rel2$Temperature, sep = "\n"),
  levels = c("September\n2-8°C", "September\n-30°C", "September\n-80°C", "April\n-80°C"))

# Cleaner library and taxonomy order, better colour palette
v1<-rel2 %>%
  mutate(Library = fct_relevel(Library,
                               "C0A", "C0B", "C5A", "C5B", "C10A", "C10B", "C15A", "C15B", 
                               "T0A", "T0B", "T5A", "T5B", "T10A", "T10B", "T15A", "T15B", 
                               "E0A", "E0B", "E5A", "E5B", "E10A", "E10B", "E15A", "E15B", 
                               "L2A", "L2B", "L4A", "L4B", "L8A", "L8B", "L12A", "L12B")) %>%
  mutate(Family = fct_relevel(Family,
                              "Amalgaviridae","Astroviridae","Birnaviridae","Bunyavirales","Alsuviricetes",
                              "Haploviricotina","Amabiliviricetes","Botourmiaviridae","Leviviricetes","Mitoviridae",
                              "Nodamuvirales","Partitiviridae","Picobirnaviridae","Durnavirales","Permutotetraviridae",
                              "Picornavirales","Reovirales","Tolivirales","Ghabrivirales","Other",)) %>%
  ggplot(aes(x = Library, y = Pct_rel2, fill = Family)) +
  ylab("Percentage of RNA virus reads (%)") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,100.1), expand = c(0.005,0.005)) +
  scale_fill_manual(values = c("#BAED57","#00B0F0","#FF56B1","#FF7E79","#FFFE99",
                               "#5E5E5E","#EDE0D4","#D4BCA7","#C89B7C","#9C6644",
                               "#91D2FC","#4DC7C9","#00A1A4","#006466","#7A81FF",
                               "#6B4785","#D883FF","#337FC3","#941651","#8CD791")) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=64, color="#000000", margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=16, color="#000000", margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=16, color="#000000", margin = margin(0, 0, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=8, color="#000000"),
    axis.text.y = element_text(size=14, color="#000000"),
    legend.position = "right"
  ) 
ggsave("v1_RelAbundFam_clean.pdf", width = 9, height = 6, units = "in")

# Collected:Temperature
z1<-v1 + facet_grid(cols = vars(FacetGroup), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("z1_RelAbundFam_facet.pdf", width = 10, height = 7, units = "in")

#######################################################################################################################

dirt2<-melt(virome2, id.vars=c("Library", "Temperature", "Group", "Time", "Collected", 
                               "Richness", "Shannon", "Shannon.effective", "TotalReads", "SumReads"),
            variable.name="Family", value.name="Count")

dirt2$Count[is.na(dirt2$Count)] <- 0
dirt2$Abundance<-(dirt2$Count/dirt2$TotalReads)*100
dirt2$Count <- as.numeric(dirt2$Count)

Abund<-aggregate(Abundance ~ Library + Temperature + Group + Time + Collected +  
                   Richness +	Shannon +	Shannon.effective + TotalReads + SumReads, data=dirt2, FUN=sum)

Abund$Group <- factor(Abund$Group, levels = c("2-8°C", "-30°C", "-80°C", "-80°C (long)"))
Abund$Time <- factor(Abund$Time, levels = c("1 day", "5 days", "10 days", "15 days", "4 weeks", "8 weeks", "12 weeks"))

v2<- Abund %>%
  mutate(Library = fct_relevel(Library,
                               "C0A", "C0B", "C5A", "C5B", "C10A", "C10B", "C15A", "C15B", 
                               "T0A", "T0B", "T5A", "T5B", "T10A", "T10B", "T15A", "T15B", 
                               "E0A", "E0B", "E5A", "E5B", "E10A", "E10B", "E15A", "E15B", 
                               "L2A", "L2B", "L4A", "L4B", "L8A", "L8B", "L12A", "L12B")) %>%
  ggplot(aes(x = Library, y = Abundance, fill = Group)) +
  ylab("Abundance of RNA viral reads (%)") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,2), expand = c(0.01,0.01)) +
  scale_fill_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=32, color="#000000", margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=14, color="#000000", margin = margin(0, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=14, color="#000000", margin = margin(0, 20, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=12, color="#000000"),
    axis.text.y = element_text(size=12, color="#000000"),
    legend.position = "right"
  ) 
ggsave("v2_TotalAbund.pdf", width = 13, height = 6, units = "in")

f2<-v2 + facet_grid(cols = vars(Time), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 13))
ggsave("f2_TotalAbund_facet.pdf", width = 13, height = 6, units = "in")

v3<- Abund %>%
  mutate(Library = fct_relevel(Library,
                               "C0A", "C0B", "C5A", "C5B", "C10A", "C10B", "C15A", "C15B", 
                               "T0A", "T0B", "T5A", "T5B", "T10A", "T10B", "T15A", "T15B", 
                               "E0A", "E0B", "E5A", "E5B", "E10A", "E10B", "E15A", "E15B", 
                               "L2A", "L2B", "L4A", "L4B", "L8A", "L8B", "L12A", "L12B")) %>%
  ggplot(aes(x = Library, y = Shannon, fill = Group)) +
  ylab("Shannon diversity") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,2), expand = c(0.01,0.01)) +
  scale_fill_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=32, color="#000000", margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=12, color="#000000", margin = margin(0, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=14, color="#000000", margin = margin(0, 20, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=12, color="#000000"),
    axis.text.y = element_text(size=12, color="#000000"),
    legend.position = "right"
  ) 
ggsave("v3_TotalShannon.pdf", width = 13, height = 6, units = "in")

f3<-v3 + facet_grid(cols = vars(Time), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 13))
ggsave("f3_TotalShannon_facet.pdf", width = 13, height = 6, units = "in")

z2<-ggarrange(f2 + theme(axis.text.x  = element_text(size = 14),
                         axis.text.y  = element_text(size = 16, hjust = 1),
                         axis.title.x  = element_text(size = 18, color = "white"),
                         axis.title.y  = element_text(size = 18),
                         axis.ticks.x  = element_blank(),
                         plot.margin = unit(c(1, 1, 0, 1), "cm"),
                         legend.position = "none"),
              f3 + theme(axis.text.x  = element_text(size = 14),
                         axis.text.y  = element_text(size = 16, hjust = 1),
                         axis.title.x  = element_text(size = 18, margin = unit(c(0.5, 0, 0, 0), "cm")),
                         axis.title.y  = element_text(size = 18),
                         axis.ticks.x  = element_blank(),
                         plot.margin = unit(c(0, 1, 1, 1), "cm"),
                         legend.position = "none"),
              heights = c(6,6),
              nrow=2, ncol=1, labels = "AUTO",
              font.label = list(size = 24, color = "black", face = "bold"),
              hjust = c(-1.5, -1.5), vjust = c(2.7, -0.6), align="v")
ggsave("z2_CombinedPlots.pdf", width = 14, height = 12, units = "in")

#######################################################################################################################

data2 <- read.csv("metadata2.csv",na.strings=c("", "NA"), header=TRUE,sep=",")
summary(data2)

soil2<-melt(data2, id.vars=c("Library", "Temperature", "Group", "Time_grouped", "Time_days", "Collected", 
                             "Concentration", "Mass", "Norm_Conc", "RIN", "Yield", "TotalReads", "SumReads"),
            variable.name="Family", value.name="Count")

soil2$Count[is.na(soil2$Count)] <- 0
soil2$Abundance<-(soil2$Count/soil2$TotalReads)*100
soil2$Abundance2<-(soil2$Abundance)*1000000 #This should be above 0, I have modified slightly
soil2$Abunance2
soil2$Count <- as.numeric(soil2$Count)

Abund<-aggregate(Abundance ~ Library + Temperature + Group + Time_grouped + Time_days + Collected + 
                   Concentration + Mass + Norm_Conc + RIN + Yield + TotalReads + SumReads, data=soil2, FUN=sum)

###################
### RHEA SCRIPT ###
###################

# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count=sum(x[x>1]^0)
  return(count)
}

# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total=sum(x)
  se=-sum(x[x>0]/total*log(x[x>0]/total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total=sum(x)
  se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
  return(se)
}

####################
### /RHEA SCRIPT ###
####################

# DCAST - Remake otu table wide formq with Abundance2
otu_table<- dcast(soil2, Family~Library, fun.aggregate = mean, value.var = "Abundance2")
otu_table$Family<-NULL

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]
otu_table[is.na(otu_table)] <- 0

# Order and transpose OTU-table
my_otu_table <- otu_table[,order(names(otu_table))] 
my_otu_table <-data.frame(t(my_otu_table))

# Apply diversity functions to table#
otus_div_stats<-data.frame(my_otu_table[,0]) 
otus_div_stats$Richness<-apply(my_otu_table,1,Species.richness)
otus_div_stats$Shannon<-apply(my_otu_table,1,Shannon.entropy)
otus_div_stats$Shannon.effective<-apply(my_otu_table,1,Shannon.effective)

otus_div_stats
write.table(otus_div_stats, file="AlphaDiversity2.txt")

# Add in MetaData
df<-tibble::rownames_to_column(otus_div_stats, "Library")
df2<-merge(x = data2, y = df, by = "Library", all = TRUE)

dflong<-melt(df2, id.vars=c("Library", "Temperature", "Group", "Time_grouped", "Time_days", "Collected", 
                            "Concentration", "Mass", "Norm_Conc", "RIN", "Yield", "TotalReads", "SumReads"))

#######################################################################################################################

# Visualising viral read abundance
p1<-Abund%>%
  mutate(Group = fct_relevel(Group,
                             "2-8°C", "-30°C", "-80°C", "-80°C (long)")) %>%
  ggplot(aes(Group, Abundance, fill=Group)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,2), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Group), size=1) +
  xlab("Temperature") + 
  ylab("Abundance") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("p1_AbundanceByTemp.pdf", width = 5, height = 6, units = "in")

p2<-ggplot(Abund, aes(factor(Time_days), Abundance, fill=factor(Time_days))) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,2), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Abundance") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("p2_AbundanceByTime.pdf", width = 5, height = 6, units = "in")

p3<-ggplot(Abund, aes(Time_days, Abundance)) +
  scale_y_continuous(limits = c(0,2), expand = c(0.025,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0.025,0), breaks = seq(0,100, by = 20)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Abundance") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("p3_AbundanceByTime2.pdf", width = 5, height = 6, units = "in")

p15<-ggarrange(p1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
               p3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
               nrow=1, ncol=2, align="h")
ggsave("p15_AbundanceByVariables.pdf", width = 10, height = 6, units = "in")

# Visualising viral richness (traditionally species, in this case, "taxa")
r1<-dflong[which(dflong$variable=="Richness"),]%>%
  mutate(Group = fct_relevel(Group,
                             "2-8°C", "-30°C", "-80°C", "-80°C (long)")) %>%
  ggplot(aes(Group, value, fill=Group)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Group), size=1) +
  xlab("Temperature") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("r1_RichnessByTemp.pdf", width = 5, height = 6, units = "in")

r2<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(factor(Time_days), value, fill=factor(Time_days))) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("r2_RichnessByTime.pdf", width = 5, height = 6, units = "in")

r3<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(Time_days, value)) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0.025,0), breaks = seq(0,100, by = 20)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("r3_RichnessByTime2.pdf", width = 5, height = 6, units = "in")

r15<-ggarrange(r1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")),
               r3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
               nrow=1, ncol=2, align="h")
ggsave("r15_RichnessByVariables.pdf", width = 10, height = 6, units = "in")

# Visualising Shannon diversity
s1<-dflong[which(dflong$variable=="Shannon"),]%>%
  mutate(Group = fct_relevel(Group,
                             "2-8°C", "-30°C", "-80°C", "-80°C (long)")) %>%
  ggplot(aes(Group, value, fill=Group)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,2), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Group), size=1) +
  xlab("Temperature") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("s1_ShannonByTemp.pdf", width = 5, height = 6, units = "in")

s2<-ggplot(dflong[which(dflong$variable=="Shannon"),], aes(factor(Time_days), value, fill=factor(Time_days))) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,2), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("s2_ShannonByTime.pdf", width = 5, height = 6, units = "in")

s3<-ggplot(dflong[which(dflong$variable=="Shannon"),], aes(Time_days, value)) +
  scale_y_continuous(limits = c(0,2), expand = c(0.025,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0.025,0), breaks = seq(0,100, by = 20)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("s3_ShannonByTime2.pdf", width = 5, height = 6, units = "in")

s15<-ggarrange(s1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")),
               s3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
               nrow=1, ncol=2, align="h")
ggsave("s15_ShannonByVariables.pdf", width = 10, height = 6, units = "in")

# Visualising effective Shannon diversity
# "True" Shannon diversity - exponent of the Shannon diversity index
# ie, number of equally common species needed to produce given Shannon diversity index
t1<-dflong[which(dflong$variable=="Shannon.effective"),]%>%
  mutate(Group = fct_relevel(Group,
                             "2-8°C", "-30°C", "-80°C", "-80°C (long)")) %>%
  ggplot(aes(Group, value, fill=Group)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,6), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Group), size=1) +
  xlab("Temperature") + 
  ylab("Effective Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("t1_TrueDivByTemp.pdf", width = 5, height = 6, units = "in")

t2<-ggplot(dflong[which(dflong$variable=="Shannon.effective"),], aes(factor(Time_days), value, fill=factor(Time_days))) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,6), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Effective Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("t2_TrueDivByTime.pdf", width = 5, height = 6, units = "in")

t3<-ggplot(dflong[which(dflong$variable=="Shannon.effective"),], aes(Time_days, value)) +
  scale_y_continuous(limits = c(0,6), expand = c(0.025,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0.025,0), breaks = seq(0,100, by = 20)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Temperature), size=1) +
  xlab("Time stored (days)") + 
  ylab("Effective Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DCE3E6","#C8CFD2","#9EA7AA","#889296","#747F82","#4A565A","#354146")) + 
  scale_color_manual(values = c("#F79800","#B41851","#7031AF","#00660E")) +
  theme(legend.position = "none")
ggsave("t3_TrueDivByTime2.pdf", width = 5, height = 6, units = "in")

t15<-ggarrange(t1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")),
               t3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(5, 0 ,0, 0)),
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
               nrow=1, ncol=2, align="h")
ggsave("t15_TrueDivByVariables.pdf", width = 10, height = 6, units = "in")

#######################################################################################################################

# Combine into one figure
z3<-ggarrange(p15, r15, s15, t15,
              nrow=2, ncol=2, align="hv")
ggsave("z3_DiversityByVariables.pdf", width = 20, height = 12, units = "in")

#######################################################################################################################

######################
### BETA DIVERSITY ###
######################

tmp<- dcast(soil2, Family~Library, value.var = "Abundance")

# Make a string of the virus families
tax<-tmp$Family
as.data.frame(tax)
TAX<-as.matrix(tax)
TAX2 = tax_table(TAX)

# Make OTU table
tmp$Family<-NULL
colnames(tmp)<-c("sa1", "sa2", "sa3", "sa4", "sa5", "sa6", "sa7", "sa8", 
                 "sa9", "sa10", "sa11", "sa12", "sa13", "sa14", "sa15", "sa16", 
                 "sa17", "sa18", "sa19", "sa20", "sa21", "sa22", "sa23", "sa24", 
                 "sa25", "sa26", "sa27", "sa28", "sa29", "sa30", "sa31", "sa32")
OTU = otu_table(tmp, taxa_are_rows = TRUE)

# Create a metadata table called "map"
map_tmp<-(df2[,c(1:9)])
map <- sample_data(map_tmp)
str(map)

# Make it into a phyloseq object
physeq = phyloseq(OTU, TAX2, map)
str(physeq)

#######################################################################################################################

# Ordinate
soil2_bray <- ordinate(
  physeq = physeq,
  method = "NMDS",
  distance = "bray")

# Check nMDS plots
b1<-plot_ordination(physeq, soil2_bray, color="Temperature") + geom_point(alpha = 1, size = 5, stroke = 2) + 
  scale_colour_manual(values=c("#F79800","#B41851","#354875"),
                      breaks=c("2-8°C","-30°C","-80°C")) + 
  stat_ellipse() +
  theme(legend.position = "right")
ggsave("b1_BetaDiversity_nMDS_Temp.pdf", width = 11, height = 7, units = "in")

b2<-plot_ordination(physeq, soil2_bray, color="Group") + geom_point(alpha = 1, size = 5, stroke = 2) + 
  scale_colour_manual(values=c("#F79800","#B41851","#7031AF","#00660E"),
                      breaks=c("2-8°C","-30°C","-80°C","-80°C (long)")) + 
  stat_ellipse() +
  theme(legend.position = "right")
ggsave("b2_BetaDiversity_nMDS_Group.pdf", width = 11, height = 7, units = "in")

# For this figure, you may get the following warning message:
# Removed 3 rows containing missing values or values outside the scale range (`geom_path()`).
# It's just warning you there are not enough points to draw an ellipse around 3 of the tested categories (in this case, Time_grouped)
# This is because the samples taken at 4, 8, and 12 weeks have only two data points each - limitation of the project design
# So don't worry about it for this one :)
b3<-plot_ordination(physeq, soil2_bray, color="Time_grouped") + geom_point(alpha = 1, size = 5, stroke = 2) + 
  scale_colour_manual(values=c("#7C0A26", "#B91B1E", "#F28225", "#FFC845", "#78B773", "#3AA5A3", "#4D5E80"),
                      breaks = c("1 day","5 days","10 days", "15 days", "4 weeks", "8 weeks", "12 weeks")) + 
  stat_ellipse() +
  theme(legend.position = "right")
ggsave("b3_BetaDiversity_nMDS_Time.pdf", width = 11, height = 7, units = "in")

z4<-plot_ordination(physeq, soil2_bray, color="Collected") + geom_point(alpha = 1, size = 5, stroke = 2) + 
  scale_colour_manual(values=c("#00660E", "#BB6ABA"),
                      breaks=c("April", "September")) + 
  stat_ellipse() +
  theme(legend.position = "right")
ggsave("z4_BetaDiversity_nMDS_Month.pdf", width = 11, height = 7, units = "in")

#######################################################################################################################

# Statistics on beta diversity using adonis

soil_bray <- phyloseq::distance(physeq, method = "bray")
sampledf <- data.frame(sample_data(physeq))

# Temperature
vegan::adonis2(soil_bray ~ Temperature, data = sampledf)

bd.temp <- betadisper(soil_bray, sampledf$Temperature)
anova(bd.temp)

TukeyHSD(bd.temp)

# Collection month
vegan::adonis2(soil_bray ~ Collected, data = sampledf)

bd.month <- betadisper(soil_bray, sampledf$Collected)
anova(bd.month)

TukeyHSD(bd.month)


#######################################################################################################################

# Recommendation, but not explored today: Run these statistics on subset of data excluding April libraries
# Do storage conditions significantly affect the abundance and alpha/beta diversity of only the September samples?

#######################################################################################################################





