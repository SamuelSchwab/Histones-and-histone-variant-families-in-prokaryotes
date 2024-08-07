import matplotlib.pyplot as plt
import numpy as np

def histogram(x, bins, density: bool = False, cumulative: bool = False,
	divider: bool = False, divider_value: float = 59.5, histtype: str = "bar",
	tickFontSize: int = 15, labelFontSize: int = 25, xlim: tuple = (0,100),
	filename: str = "histogram.png"):

	fig1, ax1 = plt.subplots(figsize=(10, 7))
	n = plt.hist(x, bins=bins, cumulative=cumulative, histtype=histtype, density=density)
	if cumulative:
		ymax = 1
		plt.yticks(np.linspace(0, ymax, num=11))
		if density:
			ax1.set_ylabel("Density (cumulative)", fontsize = labelFontSize)
		else:
			ax1.set_ylabel("Count (cumulative)", fontsize = labelFontSize)
	elif density:
		ymax = max(n[0]) + (1/10) * max(n[0])
		plt.yticks(np.linspace(0, ymax, num=11))
		ax1.set_ylabel("Density", fontsize = labelFontSize)
	else:
		ymax = max(n[0]) + 5
		ax1.set_ylabel("Count", fontsize = labelFontSize)
	if divider:
		plt.plot([divider_value, divider_value], [0, ymax], color="black", linestyle='solid', linewidth=4)
	ax1.margins(0)
	plt.ylim((0,ymax))
	plt.xlim(xlim)
	ax1.tick_params(axis='x', labelsize=tickFontSize)
	ax1.tick_params(axis='y', labelsize=tickFontSize)
	ax1.set_xlabel("Length (aa)", fontsize = labelFontSize)
	fig1.tight_layout()
	plt.savefig(filename, dpi=300)

def main():

	archaea = []
	with open("archaea_ftf_length.list", "r") as file:
		for item in file:
			archaea.append(int(item.strip()))

	bacteria = []
	with open("bacteria_ftf_length.list", "r") as file:
		for item in file:
			bacteria.append(int(item.strip()))

	bins = np.arange(0,121, 0.5)
	histogram(archaea, bins, density = True, cumulative = True, divider = True, divider_value = 59, histtype = "bar", xlim=(40,120), filename="archaea.png")
	histogram(bacteria, bins, density = True, cumulative = True, divider = True, divider_value = 59, histtype = "bar", xlim=(40,120), filename="bacteria.png")

if __name__ == "__main__":
	main()