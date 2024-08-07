library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
library(dplyr)
library(ggtreeExtra)
tree = read.newick("bac120.sp_labels_2.tree")

# Remove branches without a phylum
nodes_to_drop = c(62990,62993,63022,63135,63136,64011,64019,64050,
					66544,66545,80738,83144,90766,90968,90973,90994,
					91001,91005,91012,91014,91015,91019,91020,92347,
					92541,92543,92908,92936,94196,94247,94249,94285,
					94289,94486,103079,103086,103202,103210,103264,
					103585,103588,104650,104668,105051,105053,105054,
					106004,106077,106145,106152,106219,106614,106858,
					124213,124249)

offspring_nodes = offspring(tree, nodes_to_drop)
offspring_nodes
for (drop_list in offspring_nodes)
{
 nodes_to_drop = append(nodes_to_drop, drop_list)
}
nodes_to_drop
tree = drop.tip(tree, nodes_to_drop)
tree = drop.tip(tree, c("s__UBA5359_sp002410925", "s__Muirbacterium_halophilum", "s__SLNR01_sp007132905",
													"s__RUG730_sp900321865", "s__JAHJDO01_sp018812485", "s__DRYD01_sp011388315",
													"s__UBA1439_sp002329605", "s__UBA6266_sp002440765", "s__JAGOBX01_sp017999075",
													"s__UBA4092_sp002428325", "s__ARS69_sp002686915", "s__J088_sp003695505",
													"s__JAFGND01_sp016927185", "s__JAGOEH01_sp017997825", "s__JADJOY01_sp016713535",
													"s__BS750m-G25_sp016783505", "s__JACQOV01_sp016202395", "s__JACQTN01_sp016210785",
													"s__T1SED10-198M_sp003554345", "s__FEN-1099_sp003170555", "s__SpSt-318_sp011047235",
													"s__JAFGOL01_sp016926495", "s__CG2-30-70-394_sp001873295", "s__JACPSX01_sp016193065",
													"s__GCA-2686955_sp002686955", "s__SM23-81_sp001304015", "s__GCA-001730085_sp001730085",
													"s__BS750m-G34_sp016783345", "s__DSNQ01_sp011047215", "s__JACPUC01_sp016192455",
													"s__SpSt-165_sp011057335","s__NP-7_sp005883865","s__NP-7_sp005882275",
													"s__NP-7_sp005882445","s__SpSt-315_sp011046845","s__JACPRY01_sp016193565",
													"s__13-1-40CM-64-14_sp001917855","s__NP-4_sp005882295","s__NP-4_sp005882335",
													"s__NP-4_sp005888595","s__VGFA01_sp016869025","s__JACQHU01_sp016198805",
													"s__SpSt-190_sp011055945","s__JACQGZ01_sp016199305","s__CSP1-3_sp001443485",
													"s__CSP1-3_sp011047475"))
x = as_tibble(tree)

# Get phyla nodes
phyla = x[grep("p__", x$label),]$node

new_phyla = read.csv("new_phyla.csv")
x = left_join(x, new_phyla, by = "label")
names(x)[names(x) == "label"] = "old_label"
names(x)[names(x) == "new_label"] = "label"

df_histones = read.csv("histone_types.csv")
rownames(df_histones) = df_histones$phylum
names(df_histones)[names(df_histones) == "phylum"] = "label"
x = left_join(x, df_histones, by = "label")

# Collapse phyla nodes
tree = as.treedata(x)
p = ggtree(tree, layout="circular", size=0.5, linetype=1) + geom_nodelab(aes(subset=(node %in% phyla)), hjust=-0.05, size=3.5, fontface=0, align=TRUE, family='mono', linetype = "dotted", linesize = 0.2)# + geom_nodelab(aes(label=node), hjust=-0.05, size=0.5, fontface=0) + xlim(0,1.3)
for (phylum in phyla)
{
  p = collapse(p, node=phylum)
}


df_histones = read.csv("histone_types.csv")
rownames(df_histones) = df_histones$phylum
df_histones$phylum = NULL

df_test = data.frame(df_histones[c("Fold.to.fold")],df_histones[c("Bacterial.dimer")],
											df_histones[c("ZZ.finger")],df_histones[c("Double.histone.DUF1931")],df_histones[c("IHF.histone")],
											df_histones[c("Rdgc.histone")],df_histones[c("TM.histone")],df_histones[c("Beta.propeller")],
											df_histones[c("Undefined")])
p = p + new_scale_fill()
p = gheatmap(p, df_test, offset=0.85, width=0.3, font.size=0.5, high="#1b9e77", low="white", color="grey") + theme(legend.position = "none")

ggsave("cladogram.pdf", p, width = 15, height = 15)

# Change colors and add legend in Adobe Illustrator