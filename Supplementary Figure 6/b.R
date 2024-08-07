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

df_new_labels = read.csv("b_new_labels.csv")
df_categories = read.csv("b_categories.csv")
x = left_join(x, df_new_labels, by = "label")
x = left_join(x, df_categories, by = "label")
names(x)[names(x) == "label"] = "protein"
names(x)[names(x) == "organism_protein"] = "label"

  d_phyla = data.frame(Phylum=c("Patescibacteria", "Actinobacteriota", "Cyanobacteria",
  							"Myxococcota", "Elusimicrobiota", "Myxococcota",
  							"Myxococcota", "Planctomycetota", "Acidobacteriota",
  							"Bdellovibrionota", "Myxococcota", "Planctomycetota",
  							"Planctomycetota", "Patescibacteria", "Margulisbacteria",
  							"Planctomycetota", "UBA10199", "Planctomycetota",
  							"UBA10199", "Spirochaetota", "UBA10199",
  							"Myxococcota_A", "Planctomycetota", "Myxococcota_A",
  							"Myxococcota", "Planctomycetota", "Bdellovibrionota",
  							"Bdellovibrionota", "Chlamydiota", "Myxococcota_A",
  							"Planctomycetota", "Archaea (superkingdom)", "Archaea (superkingdom)",
  							"Archaea (superkingdom)", "Archaea (superkingdom)", "Archaea (superkingdom)",
  							"Archaea (superkingdom)", "Archaea (superkingdom)", "Archaea (superkingdom)",
  							"Archaea (superkingdom)", "Archaea (superkingdom)", "Archaea (superkingdom)",
  							"Spirochaetota", "Archaea (superkingdom)", "Myxococcota",
  							"Archaea (superkingdom)", "Archaea (superkingdom)", "Archaea (superkingdom)",
  							"Archaea (superkingdom)", "Archaea (superkingdom)", "Archaea (superkingdom)",
  							"Elusimicrobiota", "Archaea (superkingdom)", "Planctomycetota",
  							"Myxococcota", "Spirochaetota", "Proteobacteria",
  							"Bdellovibrionota", "Spirochaetota",
  							"Bdellovibrionota"),
 						node=c(2001, 2012, 2051,
 							2059, 2022, 2068,
 							2078, 2087, 2119,
 							2095, 1989, 1718,
 							1717, 1714, 1980,
 							1981, 1977, 1941,
 							1728, 1730, 1934,
 							1923, 1869, 1910,
 							1888, 1844, 1807,
 							1758, 1784, 2175,
 							2161, 2305, 2226,
 							2152, 2160, 2724,
 							1689, 2728, 1626,
 							1618, 1569, 1620,
 							1483, 1554, 2056,
 							1422, 1563, 1560,
 							1388, 1379, 1370,
 							2018, 2129, 1722,
 							1708, 1713, 1864,
 							1838, 1753,
 							1836))


colors = c("#8dd3c7","#fbb4ae","#bebada","#ffffb3","#fed9a6","#b3cde3","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")

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


ggsave("b_tree_fill.pdf", p, width = 45, height = 30, limitsize = FALSE, bg = "transparent")