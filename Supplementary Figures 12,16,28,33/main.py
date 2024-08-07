import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def volcanii_ammar():

	data_name = "volcanii_ammar"
	unit = "FPKM"
	data = pd.read_csv("ammar.txt", sep="\t")
	data = data.sort_values("FPKM", ascending=False, ignore_index=True)

	data = data[data["protein_id"].str.contains("CDS") == True]
	data = data.reset_index(drop=True)

	y = data["FPKM"].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["FtF", "FtF(pHV4)", "Double", "RdgC hist."]
	bars.append(data.index[data["protein_id"] == "CDS_HVO_0196"].tolist()[0])
	bars.append(data.index[data["protein_id"] == "CDS_HVO_A0023"].tolist()[0])
	bars.append(data.index[data["protein_id"] == "CDS_HVO_0520+"].tolist()[0])
	bars.append(data.index[data["protein_id"] == "CDS_HVO_2265++"].tolist()[0])

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def volcanii_finn():

	data_name = "volcanii_finn"
	unit = "TPM"
	data = pd.read_csv("finn.txt", sep="\t")
	data = data.sort_values("mean_H26", ascending=False, ignore_index=True)
	y = data["mean_H26"].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["FtF", "FtF(pHV4)", "Double", "RdgC hist."]
	bars.append(data.index[data["transcript_id"] == "CDS_HVO_0196"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "CDS_HVO_A0023"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "CDS_HVO_0520+"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "CDS_HVO_2265++"].tolist()[0])

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def salinarum_sakrikar():

	data_name = "salinarum_sakrikar"
	unit = "reads"
	data = pd.read_csv("sakrikar.csv", sep="\t")
	data = data.assign(avg_WTCM=data.loc[:, ["WTCMA", "WTCMB", "WTCMC", "WTCM_1", "WTCM_2", "WTCM_3"]].mean(axis=1))
	data = data.assign(avg_WTM=data.loc[:, ["WTMA","WTMC","WTM_1","WTM_2","WTM_3"]].mean(axis=1))
	data = data.sort_values("avg_WTCM", ascending=False, ignore_index=True)
	y = data["avg_WTCM"].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["FtF", "FtF(2)", "FtF(3)", "Double", "Lrp"]
	bars.append(data.index[data["transcript_id"] == "VNG_RS08870"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "VNG_RS10710"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "VNG_RS12150"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "VNG_RS00550"].tolist()[0])
	bars.append(data.index[data["transcript_id"] == "VNG_RS08105"].tolist()[0])

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit


def salinarum_lopez():

	data_name = "salinarum_lopez"
	unit = "TPM"
	data = pd.read_csv("lopez.csv", sep="\t")
	data = data.assign(avg_timepoint1=data.loc[:, ["replicate1_timepoint1", "replicate2_timepoint1", "replicate3_timepoint1"]].mean(axis=1))
	data = data.assign(avg_timepoint2=data.loc[:, ["replicate1_timepoint2", "replicate2_timepoint2", "replicate3_timepoint2"]].mean(axis=1))
	data = data.assign(avg_timepoint3=data.loc[:, ["replicate1_timepoint3", "replicate2_timepoint3", "replicate3_timepoint3"]].mean(axis=1))
	data = data.assign(avg_timepoint4=data.loc[:, ["replicate1_timepoint4", "replicate2_timepoint4", "replicate3_timepoint4"]].mean(axis=1))

	data = data.sort_values("avg_timepoint1", ascending=False, ignore_index=True)
	y = data["avg_timepoint1"].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["FtF", "Double", "FtF(2)", "FtF(3)", "Lrp"]
	bars.append(data.index[data["Gene_name"] == "VNG2273H"].tolist()[0])
	bars.append(data.index[data["Gene_name"] == "VNG0134G"].tolist()[0])
	bars.append(data.index[data["Gene_name"] == "VNG6047H"].tolist()[0])
	bars.append(data.index[data["Gene_name"] == "VNG6479H"].tolist()[0])
	bars.append(data.index[data["Gene_name"] == "VNG2094G"].tolist()[0])

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def kodakarensis_jager():

	data_name = "kodakarensis_jager"
	unit = "reads"
	data = pd.read_csv("jager.csv", sep="\t")
	data = data.sort_values("baseMeanA_Pexp", ascending=False, ignore_index=True)
	y = data["baseMeanA_Pexp"].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["HTkA", "HTkB", "FtF", "Duf"]
	bars.append(data.index[data["ID"] == "TK2289"].tolist()[0]) # HTKB
	bars.append(data.index[data["ID"] == "TK1413"].tolist()[0]) # HTKA
	bars.append(data.index[data["ID"] == "TK1040"].tolist()[0]) # FTF
	bars.append(data.index[data["ID"] == "TK0750"].tolist()[0]) # DUF

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def onnurineus_cho(condition):

	unit = "reads"
	data_name = "onnurineus_cho_" + condition
	data = pd.read_csv("cho.csv", sep="\t")
	data = data.sort_values(condition, ascending=False, ignore_index=True) # YPS MMC MMF
	print(data)
	y = data[condition].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["HToA", "HToB", "FtF"]
	bars.append(data.index[data["ID"] == "TON_0185"].tolist()[0]) # A
	bars.append(data.index[data["ID"] == "TON_1235"].tolist()[0]) # B
	bars.append(data.index[data["ID"] == "TON_0744"].tolist()[0]) # FTF

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def interrogans_xue(condition):

	unit = "reads"
	data_name = "interrogans_xue_" + condition
	data = pd.read_csv("xue.csv", sep="\t")
	data = data.sort_values(condition, ascending=False, ignore_index=True)
	y = data[condition].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["FtF", "Dps", "SMC"]
	bars.append(data.index[data["GENE_NAME"] == "LA2458"].tolist()[0]) # FTF
	bars.append(data.index[data["GENE_NAME"] == "LA3598"].tolist()[0]) # dps
	#bars.append(data.index[data["GENE_NAME"] == "LA4332"].tolist()[0]) # YbaB
	bars.append(data.index[data["GENE_NAME"] == "LA1309"].tolist()[0]) # SMC
	#bars.append(data.index[data["GENE_NAME"] == "LA1447"].tolist()[0]) # LexA

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def bacteriovorus_karunker(condition):

	unit = "FPKM"
	data_name = "bacteriovorus_karunker_" + condition
	data = pd.read_csv("karunker.csv", sep="\t")
	data = data.sort_values(condition, ascending=False, ignore_index=True)
	y = data[condition].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["Bact. dimer", "ZZ", "IHFa", "IHFb", "HUa", "HUb", "SMC"]
	bars.append(data.index[data["locus"] == "Bd0055"].tolist()[0]) # Bacterial dimer
	bars.append(data.index[data["locus"] == "Bd3044"].tolist()[0]) # ZZ
	bars.append(data.index[data["locus"] == "Bd1639"].tolist()[0]) # IHFa
	bars.append(data.index[data["locus"] == "Bd0711"].tolist()[0]) # IHFa
	bars.append(data.index[data["locus"] == "Bd2104"].tolist()[0]) # HUa
	bars.append(data.index[data["locus"] == "Bd3382"].tolist()[0]) # HUb
	bars.append(data.index[data["locus"] == "Bd1158"].tolist()[0]) # HUb

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def smithii_hansen(condition):

	unit = "reads"
	data_name = "smithii_hansen_" + condition
	data = pd.read_csv("hansen.csv", sep="\t")
	data = data.sort_values(condition, ascending=False, ignore_index=True)
	y = data[condition].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["Nuc A", "Nuc B", "Nuc C", "CC", "Alba"]
	bars.append(data.index[data["ID"] == "Msm0213"].tolist()[0]) # Nuc A
	bars.append(data.index[data["ID"] == "Msm0844"].tolist()[0]) # Nuc B
	bars.append(data.index[data["ID"] == "Msm1260"].tolist()[0]) # Nuc C
	bars.append(data.index[data["ID"] == "Msm0410"].tolist()[0]) # CC
	bars.append(data.index[data["ID"] == "Msm1245"].tolist()[0]) # Alba

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit

def cereus_kristoffersen():

	unit = "RPKM"
	data_name = "cereus_kristoffersen"
	data = pd.read_csv("kristoffersen.csv", sep="\t")
	data = data.sort_values("RPKM", ascending=False, ignore_index=True)
	y = data["RPKM"].tolist()
	x = data.index.values.tolist()

	bars = []
	labels = ["HU 1", "HU 2", "HU 3", "RgdC hist."]
	bars.append(data.index[data["locus"] == "BCE_1637"].tolist()[0]) # HU 1
	bars.append(data.index[data["locus"] == "BCE_3756"].tolist()[0]) # HU 2
	bars.append(data.index[data["locus"] == "BCE_2407"].tolist()[0]) # HU 3
	bars.append(data.index[data["locus"] == "BCE_A0013"].tolist()[0]) # RgdC
	print(bars)

	legend_labels = []
	for i, bar in enumerate(bars):
		legend_labels.append([labels[i], bar])

	return x,y,legend_labels,data_name,unit



def main():

	color_paired_list = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99"]
	color_dark_list = ["#a6cee3", '#1b9e77','#d95f02','#7570b3','#e7298a','#e6ab02']
	color_set2_list = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']
	color_paired2_list = ['#a6cee3','#33a02c','#e31a1c','#ff7f00','#6a3d9a',"#b15928", "#ffff99", "#999999"]

	color_list = color_paired2_list

	for function in [volcanii_ammar(), volcanii_finn(), salinarum_sakrikar(), salinarum_lopez(), kodakarensis_jager(),
						onnurineus_cho("YPS"), onnurineus_cho("MMF"), onnurineus_cho("MMC"), interrogans_xue("E0"),
						interrogans_xue("T90"), bacteriovorus_karunker("AP"), bacteriovorus_karunker("GP"),
						smithii_hansen("HL1"), cereus_kristoffersen()]:
		x,y,legend_labels,name,unit = function

		fig1, ax1 = plt.subplots(figsize=(10, 7))
		barlist=plt.bar(x,y, color=color_list[0], width=1)
		handles = []
		labels = []
		for i, bar_items in enumerate(legend_labels):
			histone_type = bar_items[0]
			bar = bar_items[1]
			if bar > 2000 and name != "cereus_kristoffersen":
				plt.xlim((-8,2000))
				continue
			elif bar > 2500 and name == "cereus_kristoffersen":
				plt.xlim((-8,2500))
				continue
			if bar > 1:
				for k in range(bar-3,bar+3):
					barlist[k].set_color(color_list[i+1])
			else:
				for k in range(0,4):
					barlist[k].set_color(color_list[i+1])				
			handles.append(plt.Rectangle((0,0),1,1, color=color_list[i+1]))
			labels.append(histone_type)

		# The cereus_kristoffersen plot needs a larger x range
		if name == "cereus_kristoffersen":
			plt.xlim((-8,2500))
		else:
			plt.xlim((-8,2000))
		ax1.set_yscale('log')
		ax1.set_ylabel("Expression (" + unit + ")", fontsize=20)
		ax1.set_xlabel("Genes", fontsize=20)
		ax1.tick_params(axis='x', labelsize=16)
		ax1.tick_params(axis='y', labelsize=16)
		fig1.tight_layout()
		fig1.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.97,0.97), frameon=False, fontsize = 16)#, facecolor='white', framealpha=1)
		plt.savefig(name + ".png", dpi=300)
		plt.close()



if __name__ == '__main__':
	main()