##############################################################################################
################################## TREE AESTHETICS SCRIPTS ###################################
################## Refined by Sabrina Sadiq over 5 years of trial and error ##################
########################## Tree Visualisation Workshop - 18/03/2026 ##########################
##############################################################################################

# Load packages
library(ggplot2)
library(ape)
library(ggtree)
library(phytools)
library(ggnewscale)
library(extrafont)
library(draw)

##############################################################################################

# Set working directory and import data
setwd("/path/to/folder/Resources")

metadata.zeta <- read.csv("metadata_picobirna.csv")
tree.zeta <- ape::read.tree("zetapicobirnavirus.nwk")

# Completely unedited midpoint-rooted tree
tree.zeta.r <- midpoint.root(tree.zeta)
zeta <- ggtree(tree.zeta.r)
ggsave("zeta_tree1.pdf", width = 7, height = 7.5, units = "in", limitsize = FALSE)

# Setting up: place a circle on any node > 80% bootstrap support
d.zeta <- zeta$data
d.zeta <- d.zeta[!d.zeta$isTip,]
d.zeta$label <- as.numeric(d.zeta$label)
d.zeta <- d.zeta[d.zeta$label > 80,]

# Now we're cooking
zeta <- ggtree(tree.zeta.r) %<+% metadata.zeta +
  geom_tiplab(aes(label=full_name, subset=isTip, color=broad_source),
              size = 2,
              fontface= "bold",
              hjust = -0.01) +
  scale_color_manual(values = c("#E40000", "#FED301", "#006300", "#00A1A4", 
                                "#DD00FF", "#0063FF", "#011893")) +
  theme(legend.position = "none") +
  geom_nodepoint(
    size = 1,
    colour = "black",
    data=d.zeta, aes(label=label)) +
  geom_treescale(
    x = NULL,
    y = -1,
    width = 0.5,
    offset = -2.2,
    color = "black",
    linesize = 1,
    fontsize = 5,
    family = "sans") +
  labs(title = "Zetapicobirnavirus",
       subtitle = "Picobirnaviridae, Durnavirales") +
  theme(plot.title = element_text(size = 30, face = "italic", hjust = 0.0, vjust = -1),
        plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.0, vjust = -1.5)) +
  ggplot2::xlim(0,2)
ggsave("zeta_tree2.pdf", width = 7, height = 7.5, units = "in", limitsize = FALSE)

##############################################################################################

# Import new data and midpoint root tree
metadata.monji <- read.csv("metadata_monji.csv")
tree.monji <- ape::read.tree("monjiviricetes.nwk")

tree.monji.r <- midpoint.root(tree.monji)

# Take a look at the tree
monji <- ggtree(tree.monji.r)

d.monji <- monji$data
d.monji <- d.monji[!d.monji$isTip,]
d.monji$label <- as.numeric(d.monji$label)
d.monji <- d.monji[d.monji$label > 80,]

monji <- ggtree(tree.monji.r) %<+% metadata.monji +
  geom_tiplab(aes(label=full_name),
              size = 2,
              fontface= "bold",
              color="black",
              hjust = -0.01) +
  geom_nodepoint(
    size = 1,
    colour = "black",
    data=d.monji, aes(label=label)) +
  labs(title = "Monjiviricetes",
       subtitle = "Negarnaviricota") +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, face = "italic", hjust = 0.0, vjust = -1),
        plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.0, vjust = -1.5)) +
  ggplot2::xlim(0,3)
ggsave("monji_tree1.pdf", width = 9, height = 14, units = "in", limitsize = FALSE)

# Add node labels
monji_nodes <- monji +
  geom_text2(
    mapping = aes(subset = !isTip,
                  label = node),
    size = 5,
    color = "darkred",
    hjust = 1,
    vjust = 1)  +
  ggplot2::xlim(0,3)
ggsave("monji_tree2.pdf", width = 9, height = 14, units = "in", limitsize = FALSE)

tree.monji.grouped <- groupClade(tree.monji.r, .node = c(130,160,177,182,191,126))
monji_grouped <- ggtree(tree.monji.grouped, aes(color = group))
ggsave("monji_tree3.pdf", width = 9, height = 14, units = "in", limitsize = FALSE)

# No detail. Colours kinda suck. Let's cook.
monji_clades <- ggtree(tree.monji.grouped) %<+% metadata.monji +
  aes(colour = group) +
  scale_color_manual(
    name = "Family",
    values = c("#000000", "#00B0F0", "#FF7E79", "#DDA447", 
               "#FF56B1", "#5A3C6F", "#0FAD24"),
    labels = c("NA", "Rhabdoviridae", "Paramyxoviridae", "Filoviridae",
               "Nyamiviridae", "Chuviridae", "Mymonaviridae")) +
  geom_tiplab(
    aes(label = full_name),
    size = 2,
    color = "black",
    fontface= "bold",
    offset = 0.001) +
  labs(title = "Monjiviricetes",
       subtitle = "Negarnaviricota") +
  theme(legend.position = "right",
        plot.title = element_text(size = 30, face = "italic", hjust = 0.0, vjust = -1),
        plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.0, vjust = -1.5)
  ) +
  geom_nodepoint(
    size = 1,
    colour = "black",
    data=d.monji, aes(label=label)) +
  geom_treescale(
    x = NULL,
    y = -1,
    width = 0.5,
    offset = -2.5,
    color = "black",
    linesize = 1,
    fontsize = 5,
    family = "sans") +
  ggplot2::xlim(0,3)
ggsave("monji_tree4.pdf", width = 9, height = 14, units = "in", limitsize = FALSE)

##############################################################################################

# Blank circular tree
monji_circular <- ggtree(tree.monji.r, branch.length='none', layout='circular') %<+% metadata.monji
ggsave("monji_tree5.pdf", width = 8, height = 8, units = "in", limitsize = FALSE)

# Heatmap for novel vs reference sequence
novel.monji <- data.frame("novel" = metadata.monji[,c("novel")])
rownames(novel.monji) <- metadata.monji$seq_ID

monji_ns <- monji_circular + new_scale_fill() 
monji_novel <- gheatmap(monji_ns, novel.monji,
                        offset = -6,
                        width = 1,
                        color = NA,
                        colnames = FALSE) +
  scale_fill_manual(name = "Viral sequence source",
                    values = alpha(c("#A70000", "#000000")),
                    breaks = c("Y", "N"),
                    labels = c("Novel (this study)", "Reference (NCBI)")) +
  labs(title = "Monjiviricetes",
       subtitle = "Negarnaviricota") +
  theme(legend.position = "right",
        plot.title = element_text(size = 42, face = "italic", hjust = 0.0, vjust = -0),
        plot.subtitle = element_text(size = 30, face = "italic", hjust = 0.0, vjust = -1))
ggsave("monji_tree6.pdf", width = 8, height = 8, units = "in", limitsize = FALSE)

# Heatmap for taxonomy
monji.fam <- data.frame("family" = metadata.monji[,c("family")])
rownames(monji.fam) <- metadata.monji$seq_ID

monji_ns2 <- monji_novel + new_scale_fill()
monji_taxonomy <- gheatmap(monji_ns2, monji.fam,
                           offset = 14, 
                           width = 0.25,
                           color = NA,
                           colnames = FALSE) +
  scale_fill_manual(name = "Family",
                    values = alpha(c("#00B0F0", "#FF7E79", "#DDA447", 
                                     "#FF56B1", "#5A3C6F", "#0FAD24")),
                    breaks = c("rhabdoviridae", "paramyxoviridae", "filoviridae",
                               "nyamiviridae", "chuviridae", "mymonaviridae"),
                    labels = c("Rhabdoviridae", "Paramyxoviridae", "Filoviridae",
                               "Nyamiviridae", "Chuviridae", "Mymonaviridae"))
ggsave("monji_tree7.pdf", width = 8, height = 8, units = "in", limitsize = FALSE)

##############################################################################################




