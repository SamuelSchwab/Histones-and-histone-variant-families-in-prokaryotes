library(ggplot2)
library(ggtree)
library(treeio)
library(ggbreak)
library(dplyr)
library(ggnewscale)
library(ggtreeExtra)
library(MetBrewer)

tree = read.newick("T2_TBE.tree",node.label='support')

x = as_tibble(tree)

df_new_labels = read.csv("a_new_labels.csv")
df_categories = read.csv("a_categories.csv")
x = left_join(x, df_new_labels, by = "label")
x = left_join(x, df_categories, by = "label")
names(x)[names(x) == "label"] = "protein"
names(x)[names(x) == "organism_protein"] = "label"


  d_phyla = data.frame(Phylum=c("Aenigmatarchaeota", "Asgardarchaeota", "Iainarchaeota", "Spirochaetota", "Nanoarchaeota", "Nanoarchaeota", "Nanoarchaeota",
 								"Bacteria (superkingdom)","Thermoplasmatota", "Thermoplasmatota", "Bacteria (superkingdom)", "Micrarchaeota", "Aenigmatarchaeota",
 								"Methanobacteriota B", "Thermoproteota", "Iainarchaeota", "Asgardarchaeota", "Halobacteriota", "Nanohaloarchaeota",
 								"Iainarchaeota", "Bacteria (superkingdom)", "Aenigmatarchaeota", "Asgardarchaeota", "Nanoarchaeota",
 								"Bacteria (superkingdom)", "EX4484-52","Hadarchaeota", "Micrarchaeota", "Nanoarchaeota",
 								"Altiarchaeota"),
 						node=c(1388, 1422, 1554, 1483, 1569, 1626, 1689,
 								1698,2129, 2152, 2161, 2160, 2728,
 								2229, 2298, 2296, 2267, 2319, 2306,
 								1618, 1553, 1560, 1370, 1378,
 								1418, 1563, 1400, 2724, 1620,
 								2304))


colors = c("#8dd3c7","#fbb4ae","#ffffb3","#bebada","#fed9a6","#b3cde3","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")

tree = as.treedata(x)

p = ggtree(tree, size=0.3, layout="circular") + xlim(-2, NA) + geom_tiplab(size=0.6, hjust=0) + geom_nodelab(aes(label=node), hjust=1.5, size=0.2)
			geom_hilight(data=d_phyla, aes(node = node, fill = Phylum), type = "rect", gradient.direction = "rt", alpha = 1) +
			scale_fill_manual(values=colors)
p$layers = rev(p$layers)

p = p + new_scale_fill()

p = p + geom_fruit(geom=geom_tile,mapping=aes(y=label, fill=Type), offset = 0.13, width=0.13) + 
			scale_fill_manual(values=c("#7fc97f","#fdc086","#ffff99","#386cb0","#f0027f"))
			


p = p + new_scale_fill()
p = p + geom_point2(aes(subset=!isTip, 
					fill=cut(support, c(0, 0.5, 0.75, 1))), 
					shape=21, size=1) +
			scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
					name='Transfer bootstrap expectation', 
					breaks=c('(0.75,1]', '(0.5,0.75]', '(0,0.5]'), 
					labels=expression(TBE>=75,50 <= TBE * " < 75", TBE < 50)) + guides(fill = guide_legend(override.aes = list(size=10))) +
			theme(legend.position=c(1.1, 0.5),
				legend.background=element_rect(fill=NA),
				panel.border = element_blank(),
				legend.key = element_blank(),
				axis.ticks = element_blank(),
				axis.text.y = element_blank(),
				axis.text.x = element_blank(),
				panel.grid = element_blank(),
				panel.grid.minor = element_blank(), 
				panel.grid.major = element_blank(),
				panel.background = element_blank(),
				plot.background = element_rect(fill = "transparent",colour = NA),
				legend.title=element_text(size=35),
				legend.text=element_text(size=28),
				legend.key.size = unit(2, "cm"),
				legend.spacing.y = unit(1, "cm"))


ggsave("a_tree_fill.pdf", p, width = 45, height = 30, limitsize = FALSE, bg = "transparent")