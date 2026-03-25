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
setwd("/path/to/file/Diversity_Resources/example1")
# Import data for quick glimpse figures
virome1 <-  read.csv("ViromeComp1.csv",na.strings=c("", "NA"), header=TRUE,sep=",")

#######################################################################################################################

rel1<-melt(virome1, id.vars=c("Library", "Environment", "Soil_type", "Land_use_broad", "Land_use", "Depth", 
                              "Richness", "Shannon", "Shannon.effective", "TotalReads", "SumReads"),
           variable.name="Family", value.name="Count")

rel1$Count[is.na(rel1$Count)] <- 0
rel1$Pct_rel1<-(rel1$Count/rel1$SumReads)*100

v1<-ggplot(rel1, aes(x = Library, y = Pct_rel1, fill = Family)) +
  ylab("Percentage of RNA virus reads (%)") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,100.1), expand = c(0.005,0.005)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=64, margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=16, margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=16, margin = margin(0, 0, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=8),
    axis.text.y = element_text(size=14),
    legend.position = "right"
  ) 
ggsave("v1_RelAbundFam.pdf", width = 9, height = 6, units = "in")

# Re-ordering variables for aesthetics later on
rel1$Soil_type <- factor(rel1$Soil_type, levels = c("Chromosol", "Vertosol", "Sodosol"))
rel1$Land_use <- factor(rel1$Land_use, levels = c("Cropping", "Pasture", "Native"))
rel1$Depth <- factor(rel1$Depth, levels = c("5 cm", "15 cm"))

# Cleaner library and taxonomy order, better colour palette
v2<-rel1 %>%
  mutate(Family = fct_relevel(Family,
                              "Astroviridae","Birnaviridae","Bunyavirales","Alsuviricetes","Monjiviricetes",
                              "Amabiliviricetes","Botourmiaviridae","Leviviricetes","Mitoviridae",
                              "Nodamuvirales","Durnavirales", "Picornavirales","Reovirales",
                              "Sobelivirales","Tolivirales","Ghabrivirales","Other")) %>%
  ggplot(aes(x = Library, y = Pct_rel1, fill = Family)) +
  ylab("Percentage of RNA virus reads (%)") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,100.1), expand = c(0.005,0.005)) +
  scale_fill_manual(values = c("#00B0F0","#FF56B1","#FF7E79","#FFFE99","#5E5E5E",
                               "#EDE0D4","#D4BCA7","#C89B7C","#9C6644",
                               "#91D2FC","#00A1A4","#6B4785","#D883FF",
                               "#BC91C1","#337FC3","#941651","#8CD791")) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=64, margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=16, margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=16, margin = margin(0, 0, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=8),
    axis.text.y = element_text(size=14),
    legend.position = "right"
  ) 
ggsave("v2_RelAbundFam_clean.pdf", width = 9, height = 6, units = "in")

f1<-v2 + facet_grid(cols = vars(Soil_type), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("f1_RelAbundFam_SoilType.pdf", width = 10, height = 7, units = "in")

f2<-v2 + facet_grid(cols = vars(Land_use_broad), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("f2_RelAbundFam_LandUse.pdf", width = 10, height = 7, units = "in")

f3<-v2 + facet_grid(cols = vars(Depth), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("f3_RelAbundFam_Depth.pdf", width = 10, height = 7, units = "in")

z1<-v2 + facet_nested(cols= vars(Soil_type, Land_use), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("z1_RelAbundFam_facet.pdf", width = 10, height = 7, units = "in")

#######################################################################################################################

dirt1<-melt(virome1, id.vars=c("Library", "Environment", "Soil_type", "Land_use_broad", "Land_use", "Depth", 
                               "Richness", "Shannon", "Shannon.effective", "TotalReads", "SumReads"),
            variable.name="Family", value.name="Count")

dirt1$Count[is.na(dirt1$Count)] <- 0
dirt1$Abundance<-(dirt1$Count/dirt1$TotalReads)*100

Abund<-aggregate(Abundance ~ Library + Environment + Soil_type + Land_use_broad + Land_use + Depth + 
                   Richness + Shannon + Shannon.effective + TotalReads + SumReads, data=dirt1, FUN=sum)

# Re-ordering, same as before
Abund$Soil_type <- factor(Abund$Soil_type, levels = c("Chromosol", "Vertosol", "Sodosol"))
Abund$Land_use <- factor(Abund$Land_use, levels = c("Cropping", "Pasture", "Native"))
Abund$Depth <- factor(Abund$Depth, levels = c("5 cm", "15 cm"))

v3<- Abund %>%
  ggplot(aes(x = Library, y = Abundance, fill = Soil_type)) +
  ylab("Abundance of RNA viral reads (%)") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,0.5), expand = c(0,0.01)) +
  scale_fill_manual(values = c("#F1A800","#226600","#0F8DFB"),
                    breaks = c("Chromosol","Vertosol","Sodosol")) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=32, color="#000000", margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=14, color="#000000", margin = margin(0, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=14, color="#000000", margin = margin(0, 20, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=12, color="#000000"),
    axis.text.y = element_text(size=12, color="#000000"),
    legend.position = "right"
  ) 
ggsave("v3_TotalAbund.pdf", width = 13, height = 6, units = "in")

f6<-v3 + facet_nested(cols= vars(Soil_type, Land_use), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("f6_TotalAbund_facet.pdf", width = 13, height = 6, units = "in")

v4<- Abund %>%
  ggplot(aes(x = Library, y = Shannon, fill = Soil_type)) +
  ylab("Shannon diversity") +
  xlab("Library") +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,2.5), expand = c(0,0.01)) +
  scale_fill_manual(values = c("#F1A800","#226600","#0F8DFB"),
                    breaks = c("Chromosol","Vertosol","Sodosol")) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold",size=32, color="#000000", margin = margin(10, 0, 10, 0)),
    axis.title.x = element_text(face="bold", size=14, color="#000000", margin = margin(0, 0, 0, 0)),
    axis.title.y = element_text(face="bold", size=14, color="#000000", margin = margin(0, 20, 0, 0)),
    axis.text.x = element_text(hjust=1, angle=45, size=12, color="#000000"),
    axis.text.y = element_text(size=12, color="#000000"),
    legend.position = "right"
  ) 
ggsave("v4_TotalShannon.pdf", width = 13, height = 6, units = "in")

f7<-v4 + facet_nested(cols= vars(Soil_type, Land_use), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right")
ggsave("f7_TotalShannon_facet.pdf", width = 13, height = 6, units = "in")

z2<-ggarrange(f6 + theme(axis.text.x  = element_text(size = 14),
                         axis.text.y  = element_text(size = 16, hjust = 1),
                         axis.title.x  = element_text(size = 18, color = "white"),
                         axis.title.y  = element_text(size = 18),
                         axis.ticks.x  = element_blank(),
                         plot.margin = unit(c(1, 1, 0, 1), "cm"),
                         legend.position = "none"),
              f7 + theme(axis.text.x  = element_text(size = 14),
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

data1 <-  read.csv("metadata1.csv",na.strings=c("", "NA"), header=TRUE,sep=",")

soil1<-melt(data1, id.vars=c("Library", "Environment", "Soil_type", "Land_use_broad", "Land_use", "Depth", 
                             "TotalReads", "SumReads"),
            variable.name="Family", value.name="Count")

soil1$Count[is.na(soil1$Count)] <- 0 # Turns empty cells to zeroes
soil1$Abundance<-(soil1$Count/soil1$TotalReads)*100 # Abundance calculated same as  before
soil1$Abundance2<-(soil1$Abundance)*1000000 # For Rhea script - needs to be above 1

Abund<-aggregate(Abundance ~ Library + Environment + Soil_type + Land_use_broad + Land_use + Depth + 
                   TotalReads + SumReads, data=soil1, FUN=sum)

###################
### RHEA SCRIPT ###
###################

# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts
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

# Lots of objects in this next section will be renamed when new data set analysed
# Marked specific calls to input data for ease of seeing how to apply script to new data
# Can change all object names (otu_table, my_otu_table, otus_div_stats, etc s) if you want...
# but that's effort and in most cases you won't be processing two distinct data sets in one session/script

# DCAST - Remake otu table wide form with Abundance2
# CALLS DATA HERE - CHANGE FOR NEW DATA SET ("soil1")
otu_table<- dcast(soil1, Family~Library, fun.aggregate = mean, value.var = "Abundance2")
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
write.table(otus_div_stats, file="AlphaDiversity1.txt")

# Add in MetaData
# CALLS DATA HERE - CHANGE FOR NEW DATA SET
df<-tibble::rownames_to_column(otus_div_stats, "Library")
df2<-merge(x = data1, y = df, by = "Library", all = TRUE)

dflong<-melt(df2, id.vars=c("Library", "Environment", "Soil_type", "Land_use_broad", "Land_use", "Depth", "TotalReads", "SumReads"))

#######################################################################################################################

# Visualising viral read abundance
p1<-ggplot(Abund, aes(Soil_type, Abundance, fill=Soil_type)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,0.5), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Soil type") + 
  ylab("Viral read bundance (% of total reads)") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F97100","#02885D","#2269D6"),
                    breaks = c("Chromosol","Vertosol","Sodosol")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("p1_AbundanceBySoil.pdf", width = 5, height = 6, units = "in")

p2<-ggplot(Abund, aes(Land_use_broad, Abundance, fill=Land_use_broad)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,0.5), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("Viral read bundance (% of total reads)") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("p2_AbundanceByUseBr.pdf", width = 5, height = 6, units = "in")

p3<-Abund %>%
  mutate(Land_use = fct_relevel(Land_use,
                                "Cropping","Native","Pasture")) %>%
  ggplot(aes(Land_use, Abundance, fill=Land_use)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,0.5), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("Viral read bundance (% of total reads)") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7","#0446A9"),
                    breaks = c("Cropping","Native","Pasture")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("p3_AbundanceByUse.pdf", width = 5, height = 6, units = "in")

p4<-Abund %>%
  mutate(Depth = fct_relevel(Depth, 
                             "5 cm", "15 cm")) %>%
  ggplot(aes(Depth, Abundance, fill=Depth)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,0.5), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Depth (cm)") + 
  ylab("Viral read bundance (% of total reads)") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#BB9984","#866A58")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("p4_AbundanceByDepth.pdf", width = 5, height = 6, units = "in")

p5<-Abund %>%
  mutate(Environment = fct_relevel(Environment, 
                                   "CC", "CN", "SP", "SN")) %>%
  ggplot(aes(Environment, Abundance, fill=Environment)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,0.5), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Environment") + 
  ylab("Viral read bundance (% of total reads)") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                    breaks = c("CC","CN","VC","VN","SP","SN")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("p5_AbundanceByEnv.pdf", width = 5, height = 6, units = "in")

p15<-ggarrange(p5 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0 ,0, 0)),
                          axis.title.y = element_text(face = "bold", margin = margin(0, 5, 0, 0))),
               p1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               p2 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               p3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               p4 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0 ,0, 0)),
                          axis.title.y = element_blank()),
               nrow=1, ncol=5, align="hv") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave("p15_AbundanceByVariables.pdf", width = 20, height = 6, units = "in")

# Visualising viral richness (traditionally species, in this case, "taxa")
r1<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(Soil_type, value, fill=Soil_type)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Soil type") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F97100","#02885D","#2269D6"),
                    breaks = c("Chromosol","Vertosol","Sodosol")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("r1_RichnessBySoil.pdf", width = 5, height = 6, units = "in")

r2<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(Land_use_broad, value, fill=Land_use_broad)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +  theme(legend.position = "none")
ggsave("r2_RichnessByUseBr.pdf", width = 5, height = 6, units = "in")

r3<-dflong[which(dflong$variable=="Richness"),]  %>%
  mutate(Land_use = fct_relevel(Land_use,
                                "Cropping","Native","Pasture")) %>%
  ggplot(aes(Land_use, value, fill=Land_use)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7","#0446A9"),
                    breaks = c("Cropping","Native","Pasture")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("r3_RichnessByUse.pdf", width = 5, height = 6, units = "in")

r4<-dflong[which(dflong$variable=="Richness"),]  %>%
  mutate(Depth = fct_relevel(Depth, 
                             "5 cm", "15 cm")) %>%
  ggplot(aes(Depth, value, fill=Depth)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Depth (cm)") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#BB9984","#866A58")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("r4_RichnessByDepth.pdf", width = 5, height = 6, units = "in")

r5<-dflong[which(dflong$variable=="Richness"),]  %>%
  mutate(Environment = fct_relevel(Environment, 
                                   "CC", "CN", "SP", "SN")) %>%
  ggplot(aes(Environment, value, fill=Environment)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,20), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Environment") + 
  ylab("Richness") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                    breaks = c("CC","CN","VC","VN","SP","SN")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("r5_RichnessByEnv.pdf", width = 5, height = 6, units = "in")

r15<-ggarrange(r5 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0 ,0, 0)),
                          axis.title.y = element_text(face = "bold", margin = margin(0, 5, 0, 0))),
               r1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               r2 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               r3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               r4 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0 ,0, 0)),
                          axis.title.y = element_blank()),
               nrow=1, ncol=5, align="hv") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave("r15_RichnessByVariables.pdf", width = 20, height = 6, units = "in")

# Visualising Shannon diversity
s1<-ggplot(dflong[which(dflong$variable=="Shannon"),], aes(Soil_type, value, fill=Soil_type)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,3), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Soil type") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F97100","#02885D","#2269D6"),
                    breaks = c("Chromosol","Vertosol","Sodosol")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("s1_ShannonBySoil.pdf", width = 5, height = 6, units = "in")

s2<-ggplot(dflong[which(dflong$variable=="Shannon"),], aes(Land_use_broad, value, fill=Land_use_broad)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,3), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("s2_ShannonByUseBr.pdf", width = 5, height = 6, units = "in")

s3<-dflong[which(dflong$variable=="Shannon"),]  %>%
  mutate(Land_use = fct_relevel(Land_use,
                                "Cropping","Native","Pasture")) %>%
  ggplot(aes(Land_use, value, fill=Land_use)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,3), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7","#0446A9"),
                    breaks = c("Cropping","Native","Pasture")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("s3_ShannonByUse.pdf", width = 5, height = 6, units = "in")

s4<-dflong[which(dflong$variable=="Shannon"),]  %>%
  mutate(Depth = fct_relevel(Depth, 
                             "5 cm", "15 cm")) %>%
  ggplot(aes(Depth, value, fill=Depth)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,3), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Depth (cm)") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#BB9984","#866A58")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("s4_ShannonByDepth.pdf", width = 5, height = 6, units = "in")

s5<-dflong[which(dflong$variable=="Shannon"),]  %>%
  mutate(Environment = fct_relevel(Environment, 
                                   "CC", "CN", "SP", "SN")) %>%
  ggplot(aes(Environment, value, fill=Environment)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,3), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Environment") + 
  ylab("Shannon diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                    breaks = c("CC","CN","VC","VN","SP","SN")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("s5_ShannonByEnv.pdf", width = 5, height = 6, units = "in")

s15<-ggarrange(s5 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0, 0, 0)),
                          axis.title.y = element_text(face = "bold", margin = margin(0, 15, 0, 0))),
               s1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               s2 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               s3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               s4 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0, 0, 0)),
                          axis.title.y = element_blank()),
               nrow=1, ncol=5, align="hv") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave("s15_ShannonByVariables.pdf", width = 20, height = 6, units = "in")

# Visualising effective Shannon diversity
# "True" Shannon diversity - exponent of the Shannon diversity index
# ie, number of equally common species needed to produce given Shannon diversity index
t1<-ggplot(dflong[which(dflong$variable=="Shannon.effective"),], aes(Soil_type, value, fill=Soil_type)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,10), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Soil type") + 
  ylab("True diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#F97100","#02885D","#2269D6"),
                    breaks = c("Chromosol","Vertosol","Sodosol")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("t1_TrueDivBySoil.pdf", width = 5, height = 6, units = "in")

t2<-ggplot(dflong[which(dflong$variable=="Shannon.effective"),], aes(Land_use_broad, value, fill=Land_use_broad)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,10), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("True diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("t2_TrueDivByUseBr.pdf", width = 5, height = 6, units = "in")

t3<-dflong[which(dflong$variable=="Shannon.effective"),]  %>%
  mutate(Land_use = fct_relevel(Land_use,
                                "Cropping","Native","Pasture")) %>%
  ggplot(aes(Land_use, value, fill=Land_use)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,10), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Land use") + 
  ylab("True diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#777777","#B7B7B7","#0446A9"),
                    breaks = c("Cropping","Native","Pasture")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("t3_TrueDivByUse.pdf", width = 5, height = 6, units = "in")

t4<-dflong[which(dflong$variable=="Shannon.effective"),]  %>%
  mutate(Depth = fct_relevel(Depth, 
                             "5 cm", "15 cm")) %>%
  ggplot(aes(Depth, value, fill=Depth)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,10), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Depth (cm)") + 
  ylab("True diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#BB9984","#866A58")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("t4_TrueDivByDepth.pdf", width = 5, height = 6, units = "in")

t5<-dflong[which(dflong$variable=="Shannon.effective"),]  %>%
  mutate(Environment = fct_relevel(Environment, 
                                   "CC", "CN", "SP", "SN")) %>%
  ggplot(aes(Environment, value, fill=Environment)) + geom_boxplot(alpha=0.6) +
  scale_y_continuous(limits = c(0,10), expand = c(0.025,0)) +
  geom_point(colour = "black", size = 1.5)+ geom_point(aes(colour = Environment), size=1) +
  xlab("Environment") + 
  ylab("True diversity") +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                    breaks = c("CC","CN","VC","VN","SP","SN")) + 
  scale_color_manual(values = c("#DE0000","#FF9F00","#004D10","#02C390","#0446A9","#1E88E5"),
                     breaks = c("CC","CN","VC","VN","SP","SN")) +
  theme(legend.position = "none")
ggsave("t5_TrueDivByEnv.pdf", width = 5, height = 6, units = "in")

t15<-ggarrange(t5 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0 ,0, 0)),
                          axis.title.y = element_text(face = "bold", margin = margin(0, 0, 0, 0))),
               t1 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               t2 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               t3 + theme(axis.title.x  = element_text(face = "bold", margin = margin(0, 0, 0, 0)),
                          axis.title.y = element_blank()),
               t4 + theme(axis.title.x  = element_text(face = "bold", margin = margin(-25, 0 ,0, 0)),
                          axis.title.y = element_blank()),
               nrow=1, ncol=5, align="hv") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave("t15_TrueDivByVariables.pdf", width = 20, height = 6, units = "in")

#######################################################################################################################

# Combine all 4 indices by all 5 ecological variables
# I edit this further in Illustrator to annotate pairwise relationships
z3<-ggarrange(p15 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
              r15 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
              s15 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
              t15 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
              nrow=4, ncol=1, align="hv") + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.2), "cm"))
ggsave("z3_AllDivByVariables.pdf", width = 20, height = 24, units = "in")

#################################################################################################################################################

##############################
### EXPLORATORY STATISTICS ###
######  & MODEL FITTING ######
##############################

### ABUNDANCE ###

###
# Sometimes something goes weird with data frame formatting, but I haven't 100% figured out what causes it
# But converting each variable into a tibble cleans it up and allows Tukey tests to be run
# Only needs to be done once per variable for abundance (Abund) and once per variable for all alpha diversity indices (df2)
# That's why you'll see it throughout abundance and richness stats, but not in Shannon and effective Shannon
# (because doing it for richness first applied it to the other alpha diversity stats too)
DATA <-as_tibble (DATA)%>%
  mutate (VARIABLE = factor (VARIABLE))
###

Abund <-as_tibble (Abund)%>%
  mutate (Environment = factor (Environment))
m.env<-lm(Abundance ~ Environment, data=Abund)
print(summary(m.env))
print(anova(m.env,test="Chisq"))

summary(glht(m.env,linfct = mcp(Environment="Tukey")))


Abund <-as_tibble (Abund)%>%
  mutate (Soil_type = factor (Soil_type))
m.type<-lm(Abundance ~ Soil_type, data=Abund)
print(summary(m.type))
print(anova(m.type,test="Chisq"))

summary(glht(m.type,linfct = mcp(Soil_type="Tukey")))


Abund <-as_tibble (Abund)%>%
  mutate (Land_use_broad = factor (Land_use_broad))
m.broad<-lm(Abundance ~ Land_use_broad, data=Abund)
print(summary(m.broad))
print(anova(m.broad,test="Chisq"))

summary(glht(m.broad,linfct = mcp(Land_use_broad="Tukey")))


Abund <-as_tibble (Abund)%>%
  mutate (Land_use = factor (Land_use))
m.use<-lm(Abundance ~ Land_use, data=Abund)
print(summary(m.use))
print(anova(m.use,test="Chisq"))

summary(glht(m.use,linfct = mcp(Land_use="Tukey")))


Abund <-as_tibble (Abund)%>%
  mutate (Depth = factor (Depth))
m.depth<-lm(Abundance ~ Depth, data=Abund)
print(summary(m.depth))
print(anova(m.depth,test="Chisq"))

summary(glht(m.depth,linfct = mcp(Depth="Tukey")))


# Model testing
mmodel_1<-lm(Abundance ~ Environment+Depth, data=Abund)
Anova(mmodel_1, type="II")

mmodel_2<-lm(Abundance ~ Soil_type+Land_use+Depth, data=Abund)
Anova(mmodel_2, type="II")

mmodel_3<-lm(Abundance ~ Soil_type+Depth, data=Abund)
Anova(mmodel_3, type="II")

mmodel_4<-lm(Abundance ~ Land_use+Depth, data=Abund)
Anova(mmodel_4, type="II")

mmodel_5<-lm(Abundance ~ Soil_type+Land_use, data=Abund)
Anova(mmodel_5, type="II")


AIC(mmodel_1, mmodel_2, mmodel_3, mmodel_4, mmodel_5, m.env, m.type, m.broad, m.use, m.depth)


# Don't use the Analysis of Variance Table from this anova - it's just to quickly get the list of models
anova(mmodel_1, mmodel_2, mmodel_3, mmodel_4, mmodel_5, m.env, m.type, m.broad, m.use, m.depth, test="Chisq")


# Compare a complex model with a simpler, nested model, ie, a subset of variables from the complex model
# If p<0.05 (*), complex model significantly better than nested model - missing variable importantly contributes to fit!
# If p>0.05, complex model not significantly different to nester model - missing value isn't contributing anything to fit
anova(mmodel_1, m.env, test="Chisq")

anova(mmodel_1, m.depth, test="Chisq")

anova(mmodel_2, mmodel_3, test="Chisq")

anova(mmodel_2, mmodel_4, test="Chisq")

anova(mmodel_2, mmodel_5, test="Chisq")

anova(m.type, mmodel_2, test="Chisq")

anova(m.use, mmodel_2, test="Chisq")

anova(m.depth, mmodel_2, test="Chisq")


### RICHNESS ###

df2 <-as_tibble (df2)%>%
  mutate (Environment = factor (Environment))
r.env<-glm(Richness ~ Environment, data=df2)
print(summary(r.env))
print(anova(r.env,test="Chisq"))

summary(glht(r.env,linfct = mcp(Environment="Tukey")))


df2 <-as_tibble (df2)%>%
  mutate (Soil_type = factor (Soil_type))
r.type<-glm(Richness ~ Soil_type, data=df2)
print(summary(r.type))
print(anova(r.type,test="Chisq"))

summary(glht(r.type,linfct = mcp(Soil_type="Tukey")))


df2 <-as_tibble (df2)%>%
  mutate (Land_use_broad = factor (Land_use_broad))
r.broad<-glm(Richness ~ Land_use_broad, data=df2)
print(summary(r.broad))
print(anova(r.broad,test="Chisq"))

summary(glht(r.broad,linfct = mcp(Land_use_broad="Tukey")))


df2 <-as_tibble (df2)%>%
  mutate (Land_use = factor (Land_use))
r.use<-glm(Richness ~ Land_use, data=df2)
print(summary(r.use))
print(anova(r.use,test="Chisq"))

summary(glht(r.use,linfct = mcp(Land_use="Tukey")))


df2 <-as_tibble (df2)%>%
  mutate (Depth = factor (Depth))
r.depth<-glm(Richness ~ Depth, data=df2)
print(summary(r.depth))
print(anova(r.depth,test="Chisq"))

summary(glht(r.depth,linfct = mcp(Depth="Tukey")))


# Model testing
rmodel_1<-glm(Richness ~ Environment+Depth, data=df2)
Anova(rmodel_1, type="II")

rmodel_2<-glm(Richness ~ Soil_type+Land_use+Depth, data=df2)
Anova(rmodel_2, type="II")

rmodel_3<-glm(Richness ~ Soil_type+Depth, data=df2)
Anova(rmodel_3, type="II")

rmodel_4<-glm(Richness ~ Land_use+Depth, data=df2)
Anova(rmodel_4, type="II")

rmodel_5<-glm(Richness ~ Soil_type+Land_use, data=df2)
Anova(rmodel_5, type="II")

# Trying some extra models because Land_use_broad was also significant for richness
# Note that Land_use and Land_use_broad are not in any models together
# That's because one is nested in the other (cropping and pasture are by definition ONLY agricultural and native = native, obviously)
rmodel_6<-glm(Richness ~ Soil_type+Land_use_broad+Depth, data=df2)
Anova(rmodel_6, type="II")

rmodel_7<-glm(Richness ~ Soil_type+Depth, data=df2)
Anova(rmodel_7, type="II")

rmodel_8<-glm(Richness ~ Land_use_broad+Depth, data=df2)
Anova(rmodel_8, type="II")

rmodel_9<-glm(Richness ~ Soil_type+Land_use_broad, data=df2)
Anova(rmodel_9, type="II")

AIC(rmodel_1, rmodel_2, rmodel_3, rmodel_4, rmodel_5, rmodel_6, rmodel_7, rmodel_8, rmodel_9, r.env, r.type, r.broad, r.use, r.depth)


anova(rmodel_1, rmodel_2, rmodel_3, rmodel_4, rmodel_5, rmodel_6, rmodel_7, rmodel_8, rmodel_9, r.env, r.type, r.broad, r.use, r.depth, test="Chisq")
# Model  1: Richness ~ Environment + Depth
# Model  2: Richness ~ Soil_type + Land_use + Depth
# Model  3: Richness ~ Soil_type + Depth
# Model  4: Richness ~ Land_use + Depth
# Model  5: Richness ~ Soil_type + Land_use
# Model  6: Richness ~ Soil_type + Land_use_broad + Depth
# Model  7: Richness ~ Soil_type + Depth
# Model  8: Richness ~ Land_use_broad + Depth
# Model  9: Richness ~ Soil_type + Land_use_broad
# Model 10: Richness ~ Environment
# Model 11: Richness ~ Soil_type
# Model 12: Richness ~ Land_use_broad
# Model 13: Richness ~ Land_use
# Model 14: Richness ~ Depth

anova(rmodel_6, rmodel_7, test="Chisq")

anova(rmodel_6, rmodel_8, test="Chisq")

anova(rmodel_6, rmodel_9, test="Chisq")

anova(rmodel_6, r.type, test="Chisq")

anova(rmodel_6, r.broad, test="Chisq")

anova(rmodel_6, r.depth, test="Chisq")

anova(rmodel_2, rmodel_3, test="Chisq")

anova(rmodel_2, rmodel_4, test="Chisq")

anova(rmodel_2, rmodel_5, test="Chisq")

anova(rmodel_2, r.type, test="Chisq")

anova(rmodel_2, r.use, test="Chisq")

anova(rmodel_2, r.depth, test="Chisq")

anova(rmodel_1, r.env, test="Chisq")

anova(rmodel_1, r.depth, test="Chisq")


### SHANNON ###

s.env<-lm(Shannon ~ Environment, data=df2)
print(summary(s.env))
print(anova(s.env,test="Chisq"))

summary(glht(s.env,linfct = mcp(Environment="Tukey")))


s.type<-lm(Shannon ~ Soil_type, data=df2)
print(summary(s.type))
print(anova(s.type,test="Chisq"))

summary(glht(s.type,linfct = mcp(Soil_type="Tukey")))


s.broad<-lm(Shannon ~ Land_use_broad, data=df2)
print(summary(s.broad))
print(anova(s.broad,test="Chisq"))

summary(glht(s.broad,linfct = mcp(Land_use_broad="Tukey")))


s.use<-lm(Shannon ~ Land_use, data=df2)
print(summary(s.use))
print(anova(s.use,test="Chisq"))

summary(glht(s.use,linfct = mcp(Land_use="Tukey")))


s.depth<-lm(Shannon ~ Depth, data=df2)
print(summary(s.depth))
print(anova(s.depth,test="Chisq"))

summary(glht(s.depth,linfct = mcp(Depth="Tukey")))


# Model testing
smodel_1<-lm(Shannon ~ Environment+Depth, data=df2)
Anova(smodel_1, type="II")

smodel_2<-lm(Shannon ~ Soil_type+Land_use+Depth, data=df2)
Anova(smodel_2, type="II")

smodel_3<-lm(Shannon ~ Soil_type+Depth, data=df2)
Anova(smodel_3, type="II")

smodel_4<-lm(Shannon ~ Land_use+Depth, data=df2)
Anova(smodel_4, type="II")

smodel_5<-lm(Shannon ~ Soil_type+Land_use, data=df2)
Anova(smodel_5, type="II")

AIC(smodel_1, smodel_2, smodel_3, smodel_4, smodel_5, s.env, s.type, s.broad, s.use, s.depth)


anova(smodel_1, smodel_2, smodel_3, smodel_4, smodel_5, s.env, s.type, s.broad, s.use, s.depth, test="Chisq")
# Model  1: Shannon ~ Environment + Depth
# Model  2: Shannon ~ Soil_type + Land_use + Depth
# Model  3: Shannon ~ Soil_type + Depth
# Model  4: Shannon ~ Land_use + Depth
# Model  5: Shannon ~ Soil_type + Land_use
# Model  6: Shannon ~ Environment
# Model  7: Shannon ~ Soil_type
# Model  8: Shannon ~ Land_use_broad
# Model  9: Shannon ~ Land_use
# Model 10: Shannon ~ Depth

anova(smodel_1, s.env, test="Chisq")

anova(smodel_1, s.depth, test="Chisq")


### TRUE DIVERSITY ###

t.env<-lm(Shannon.effective ~ Environment, data=df2)
print(summary(t.env))
print(anova(t.env,test="Chisq"))

summary(glht(t.env,linfct = mcp(Environment="Tukey")))


t.type<-lm(Shannon.effective ~ Soil_type, data=df2)
print(summary(t.type))
print(anova(t.type,test="Chisq"))

summary(glht(t.type,linfct = mcp(Soil_type="Tukey")))


t.broad<-lm(Shannon.effective ~ Land_use_broad, data=df2)
print(summary(t.broad))
print(anova(t.broad,test="Chisq"))

summary(glht(t.broad,linfct = mcp(Land_use_broad="Tukey")))

t.use<-lm(Shannon.effective ~ Land_use, data=df2)
print(summary(t.use))
print(anova(t.use,test="Chisq"))

summary(glht(t.use,linfct = mcp(Land_use="Tukey")))


t.depth<-lm(Shannon.effective ~ Depth, data=df2)
print(summary(t.depth))
print(anova(t.depth,test="Chisq"))

summary(glht(t.depth,linfct = mcp(Depth="Tukey")))


# Model testing
tmodel_1<-lm(Shannon.effective ~ Environment+Depth, data=df2)
Anova(tmodel_1, type="II")

tmodel_2<-lm(Shannon.effective ~ Soil_type+Land_use+Depth, data=df2)
Anova(tmodel_2, type="II")

tmodel_3<-lm(Shannon.effective ~ Soil_type+Depth, data=df2)
Anova(tmodel_3, type="II")

tmodel_4<-lm(Shannon.effective ~ Land_use+Depth, data=df2)
Anova(tmodel_4, type="II")

tmodel_5<-lm(Shannon.effective ~ Soil_type+Land_use, data=df2)
Anova(tmodel_5, type="II")

AIC(tmodel_1, tmodel_2, tmodel_3, tmodel_4, tmodel_5, t.env, t.type, t.broad, t.use, t.depth)


anova(tmodel_1, tmodel_2, tmodel_3, tmodel_4, tmodel_5, t.env, t.type, t.broad, t.use, t.depth, test="Chisq")
# Model  1: Shannon.effective ~ Environment + Depth
# Model  2: Shannon.effective ~ Soil_type + Land_use + Depth
# Model  3: Shannon.effective ~ Soil_type + Depth
# Model  4: Shannon.effective ~ Land_use + Depth
# Model  5: Shannon.effective ~ Soil_type + Land_use
# Model  6: Shannon.effective ~ Environment
# Model  7: Shannon.effective ~ Soil_type
# Model  8: Shannon.effective ~ Land_use_broad
# Model  9: Shannon.effective ~ Land_use
# Model 10: Shannon.effective ~ Depth

anova(tmodel_1, t.env, test="Chisq")

anova(tmodel_1, t.depth, test="Chisq")

#######################################################################################################################







