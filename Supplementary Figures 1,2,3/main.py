import json
import matplotlib.pyplot as plt
import numpy as np

def load_db(path):
    with open(path, "r") as file:
        db = json.load(file)
    return db

def save_db(db, path):
    with open(path, "w") as file:
        json.dump(db, file, indent=4)

def plot_msa_hist_cat(name, category_map, categories, xlim, ylim, binsize):
    fontsize = 23
    labelsize = 13

    msa_lengths = load_db("msas.json")
    length_dict = {}
    colours = ["#80b1d3", "#8da0cb", "#b3de69", "#fb8072", "#bebada","#bc80bd",
                "#8dd3c7", "#ffed6f","#b15928", "#66c2a5", "#a6cee3", "#fdb462",
                "#fb9a99", "#fc8d62", "#a6d854", "#fccde5", "#ffffb3"]

    for item in msa_lengths:
        length = msa_lengths[item]
        category = category_map[item]
        if category in categories:
            if category not in length_dict:
                length_dict[category] = []
            length_dict[category].append(length)
    
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    i = 0
    for category in categories:
        plt.hist(length_dict[category], np.arange(0, xlim, binsize), color=colours[i], label=category)
        i += 1
    fig1.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes, frameon = False, fontsize = 15)
    plt.plot([100,100], [0,1000], color="black", lw=2.5, alpha=1)
    
    ax1.set_xlim(0, xlim)
    ax1.set_ylim(0, ylim)
    ax1.set_xlabel("MSA depth", fontsize=fontsize)
    ax1.set_ylabel("Counts", fontsize=fontsize)
    ax1.tick_params(axis="x", labelsize=labelsize)
    ax1.tick_params(axis="y", labelsize=labelsize)
    fig1.tight_layout()
    plt.savefig(name, dpi=300)
    #plt.show()

def plot_plddt_cat(name, category_map, multimer, categories, xlim, ylim, binsize):
    fontsize = 23
    labelsize = 13

    plddt_values = load_db("plddts.json")
    plddt_dict = {}
    colours = ["#80b1d3", "#8da0cb", "#b3de69", "#fb8072", "#bebada","#bc80bd",
                "#8dd3c7", "#ffed6f","#b15928", "#66c2a5", "#a6cee3", "#fdb462",
                "#fb9a99", "#fc8d62", "#a6d854", "#fccde5", "#ffffb3"]

    for item in plddt_values[multimer]:
        plddt = plddt_values[multimer][item]
        category = category_map[item]
        if category in categories:
            if category not in plddt_dict:
                plddt_dict[category] = []
            plddt_dict[category].append(plddt)
    
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    i = 0
    for category in categories:
        if category in plddt_dict:
            plt.hist(plddt_dict[category], np.arange(0, xlim, binsize), color=colours[i], label=category)
            i += 1
    fig1.legend(loc="upper left", bbox_to_anchor=(0,1), bbox_transform=ax1.transAxes, frameon = False, fontsize = 15)
    
    ax1.set_xlim(0, xlim)
    ax1.set_ylim(0, ylim)
    ax1.set_xlabel("pLDDT", fontsize=fontsize)
    ax1.set_ylabel("Counts", fontsize=fontsize)
    ax1.tick_params(axis="x", labelsize=labelsize)
    ax1.tick_params(axis="y", labelsize=labelsize)
    fig1.tight_layout()
    plt.savefig(name, dpi=300)
    #plt.show()

def plot_pae_cat(name, category_map, multimer, categories, xlim, ylim, binsize):
    fontsize = 23
    labelsize = 13

    pae_values = load_db("paes.json")
    pae_dict = {}
    colours = ["#80b1d3", "#8da0cb", "#b3de69", "#fb8072", "#bebada","#bc80bd",
                "#8dd3c7", "#ffed6f","#b15928", "#66c2a5", "#a6cee3", "#fdb462",
                "#fb9a99", "#fc8d62", "#a6d854", "#fccde5", "#ffffb3"]

    for item in pae_values[multimer]:
        pae = pae_values[multimer][item]
        category = category_map[item]
        if category in categories:
            if category not in pae_dict:
                pae_dict[category] = []
            pae_dict[category].append(pae)
    
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    i = 0
    for category in categories:
        if category in pae_dict:
            plt.hist(pae_dict[category], np.arange(0, xlim, binsize), color=colours[i], label=category)
            i += 1
    fig1.legend(loc="upper left", bbox_to_anchor=(0,1), bbox_transform=ax1.transAxes, frameon = False, fontsize = 15)
    
    ax1.set_xlim(0, xlim)
    ax1.set_ylim(0, ylim)
    ax1.set_xlabel("ipTM", fontsize=fontsize)
    ax1.set_ylabel("Counts", fontsize=fontsize)
    ax1.tick_params(axis="x", labelsize=labelsize)
    ax1.tick_params(axis="y", labelsize=labelsize)
    fig1.tight_layout()
    plt.savefig(name, dpi=300)
    #plt.show()

def main():
    category_map = load_db("categories.json")

    categories = ['Nucleosomal', 'Face-to-face', 'DUF1931', 'Halo double', 'Coiled-coil histone',
                  'Poseidoniia double', 'Undefined', 'Bacterial dimer', 'ZZ histone']    
    minor_categories = ['Methanococcales histone', 'Beta-propeller', 'RdgC histone',
                        'Transmembrane histone', 'IHF histone', 'Thermoplasmatota',
                        'Rab GTPase', 'Nanohalo coiled-coil', 'Phage histone']

    plot_msa_hist_cat("msa_major.png",category_map,categories,9500,360,100)
    plot_msa_hist_cat("msa_minor.png",category_map,minor_categories,3500,15,100)

    for oligomer in ["monomer", "dimer", "tetramer", "hexamer"]:
        plot_plddt_cat("plddt_" + oligomer + "_major.png",category_map,oligomer,categories,100,620,1)
        plot_plddt_cat("plddt_" + oligomer + "_minor.png",category_map,oligomer,minor_categories,100,15,2)

    for oligomer in ["dimer", "tetramer", "hexamer"]:
       plot_pae_cat("pae_" + oligomer + "_minor.png",category_map,oligomer,minor_categories,1,25,0.02)
    
    plot_pae_cat("pae_dimer_major.png",category_map,"dimer",categories,1,700,0.01)
    plot_pae_cat("pae_tetramer_major.png",category_map,"tetramer",categories,1,250,0.01)
    plot_pae_cat("pae_hexamer_major.png",category_map,"hexamer",categories,1,300,0.02)


if __name__ == "__main__":
    main()