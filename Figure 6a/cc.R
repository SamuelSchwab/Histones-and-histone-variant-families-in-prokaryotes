library(ggplot2)
library(ggtree)
library(treeio)
library(ggbreak)
library(dplyr)
library(ggnewscale)

tree = read.newick("cc.tree",node.label='support')

x = as_tibble(tree)
x

df_new_labels = read.csv("new_labels.csv")
x = left_join(x, df_new_labels, by = "label")

names(x)[names(x) == "label"] = "protein"
names(x)[names(x) == "organism_protein"] = "label"

d_phyla = data.frame(Phylum=c("Methanobacteriota", "Undefined", "B1Sed10-29", "Methanobacteriota_B", "Micrarchaeota", "Nanohaloarchaeota", "Aenigmatarchaeota", "EX4484-52", "Nanoarchaeota", "Undefined", "Undinarchaeota", "Altiarchaeota", "Aenigmatarchaeota", "Iainarchaeota", "Iainarchaeota", "Iainarchaeota"),
						node=c(206, 291, 292, 293, 332, 316, 329, 311, 350, 308, 302, 358, 371, 201, 381, 386))

colors = c('#8dd3c7','#ffffb3',"#ffed6f",'#ccebc5','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd')

tree = as.treedata(x)

p = ggtree(tree, size=1, layout="circular", ladderize = TRUE) +
			geom_tiplab(size=1.5) +
			geom_hilight(data=d_phyla, aes(node = node, fill = Phylum), type = "rect", gradient.direction = "rt", alpha = 1) +
			scale_fill_manual(values=colors)

p$layers = rev(p$layers)
p2 = p + new_scale_fill()

poutput = p2 + geom_point2(aes(subset=!isTip, 
					fill=cut(support, c(0, 0.5, 0.75, 1))), 
					shape=21, size=1.5) +
			scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
					name='Transfer bootstrap expectation', 
					breaks=c('(0.75,1]', '(0.5,0.75]', '(0,0.5]'), 
					labels=expression(TBE>=75,50 <= TBE * " < 75", TBE < 50)) + guides(fill = guide_legend(override.aes = list(size=7))) +
			geom_treescale(x=2.6, y=0, fontsize=6, linesize=1) +
			theme(legend.position=c(0.9, 0.55),
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
				legend.title=element_text(size=20),
				legend.text=element_text(size=15),
				legend.key.size = unit(1, "cm"),
				legend.spacing.y = unit(0.5, "cm"))

ggsave("cc_tree.pdf", poutput, width = 20, height = 20, limitsize = FALSE, bg = "transparent")