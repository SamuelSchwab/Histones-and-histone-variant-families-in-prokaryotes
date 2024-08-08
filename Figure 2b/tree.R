library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
tree = read.newick("ar53.sp_labels.tree")

# Drop tips without a phylum
tree = drop.tip(tree, 3407)
x = as_tibble(tree)

# Get phyla nodes
x[grep("p__", x$label),]
archaea_phyla = x[grep("p__", x$label),]$node

# Rename
d = data.frame(label = c("100.0:p__Iainarchaeota; c__Iainarchaeia", "71.0:p__Hadarchaeota",
                         '100.0:p__Halobacteriota','100.0:p__Thermoplasmatota',
                         '78.0:p__Thermoproteota','100.0:p__Asgardarchaeota',
                         '100.0:p__Methanobacteriota; c__Methanobacteria; o__Methanobacteriales',
                         '100.0:p__Hydrothermarchaeota; c__Hydrothermarchaeia; o__Hydrothermarchaeales',
                         '26.0:p__Methanobacteriota_A', '100.0:p__Methanobacteriota_B; c__Thermococci',
                         '100.0:p__Micrarchaeota; c__Micrarchaeia', '100.0:p__Altiarchaeota; c__Altiarchaeia',
                         '93.0:p__Nanoarchaeota; c__Nanoarchaeia','98.0:p__Aenigmatarchaeota; c__Aenigmatarchaeia',
                         '100.0:p__Nanohaloarchaeota; c__Nanosalinia; o__Nanosalinales',
                         '100.0:p__Undinarchaeota; c__Undinarchaeia; o__Undinarchaeales',
                         '100.0:p__B1Sed10-29; c__B1Sed10-29','100.0:p__EX4484-52; c__EX4484-52; o__EX4484-52',
                         '100.0:p__SpSt-1190; c__SpSt-1190; o__SpSt-1190; f__SpSt-1190'),
               label2 = c("Iainarchaeota", "Hadarchaeota",
                          "Halobacteriota","Thermoplasmatota",
                          "Thermoproteota","Asgardarchaeota",
                          "Methanobacteriota",
                          "Hydrothermarchaeota",
                          "Methanobacteriota_A", "Methanobacteriota_B",
                          "Micrarchaeota", "Altiarchaeota",
                          "Nanoarchaeota","Aenigmatarchaeota",
                          "Nanohaloarchaeota",
                          "Undinarchaeota","B1Sed10-29","EX4484-52",
                          "SpSt-1190"))
for (row in 1:nrow(d))
{
  x[grep(d[row, "label"], x$label),]$label = d[row, "label2"]
}
tree = as.phylo(x)

# Collapse phyla nodes
p = ggtree(tree, right=TRUE, size=1.4, linetype=1) + geom_nodelab(aes(subset=(node %in% archaea_phyla)), hjust=-0.05, size=7, fontface=0)
for (phylum in archaea_phyla)
{
  p = collapse(p, node=phylum)
}

# Draw and save
p
df_histones = read.csv("tree.csv")
rownames(df_histones) = df_histones$organism
df_histones$organism = NULL
names(df_histones)[names(df_histones) == "Fold.to.fold"] = "FtF"
names(df_histones)[names(df_histones) == "Coiled.coil.bridging"] = "CC"
names(df_histones)[names(df_histones) == "Nucleosomal"] = "Nuc"
names(df_histones)[names(df_histones) == "Double.histone.DUF1931"] = "DUF1931"
names(df_histones)[names(df_histones) == "Halobacterium.double.histone"] = "Halo"
names(df_histones)[names(df_histones) == "Jannaschii.bridging"] = "Mc"
names(df_histones)[names(df_histones) == "Poseidoniaceae.double.histone"] = "Poseidon"
names(df_histones)[names(df_histones) == "Thermoplasmatota.large.histone"] = "Thermo"

p2 = p + new_scale_fill()
p3 = gheatmap(p2, df_histones, offset=0.55, width=0.4, font.size=5, high="black", low="white", color="grey", colnames_angle=-45, hjust=0) + theme(legend.position = "none") + vexpand(0.03, direction=-1) + hexpand(-0.6, direction=1) +
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "none")
p3
ggsave("tree.pdf", p3, width = 10, height = 8, limitsize = FALSE, bg = "transparent")