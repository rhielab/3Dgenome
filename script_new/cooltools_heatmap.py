import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import cooltools
import cooler
from matplotlib.ticker import EngFormatter
import cooltools.lib.plotting
from matplotlib.colors import LogNorm

#Load all the packages beforehand; all the packages are probably required, so just copy and past above when plotting with cooltools.
bp_formatter = EngFormatter('b')

#Process argument, and converthem into resolution ,prefix name and file name.
arg_list = sys.argv[1:]
output_path = arg_list[1]
output_prefix = arg_list[2]
input_file = arg_list[0] #should be an .mcool file

#This is function is required to format the plot
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


#Loading cooler file.
#Let me explain how mcool files are desinged.
#With mcool files, there are essentially made up of multiple cool files with different resoltuions, usually starting from original cool file, whichi is 2.5kb here, and 2x the resolution of the original cool file.
#So in this case, because we used 2.5kb, it would be 5kb, 10kb, 20kb.... so on for the resolutions. In this script, I used 20kb, because it is the fastest.
clr_1kb_hic=  cooler.Cooler(input_file + '::resolutions/20000')


#Location of the plot to be at. Not sure why we need both region and extents, but both are required to plot a the spot 
region = 'chr7:31,000,000-43,000,000'
start, end = 31_000_000, 43_000_000
extents = (start, end, end, start)
#This sets the limit in which the scale is to be displayed; I set it to 100 here, but you can extend to 1000 if you have enough interactions in the area.
norm = LogNorm(vmin=1, vmax=100)
#Pick color; I picked white and red since it is easier to see
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)

#This is where you define how many plots you will have, along with the size of the plot. nrows = Number of rows, ncols = Number of columns; set sharex and sharey = True since we will use them at the same time 
f, ax = plt.subplots(
    figsize=(13, 10),
    nrows=1,
    ncols=1,
    sharex=True,
    sharey=True
)

#Where actual plotting happens. clr_1kb_hic.matrix takes cooler matrices we got from above, use cmap as color map, then scale with norm, then plot upto extent
im = ax.matshow(
    clr_1kb_hic.matrix(balance=False).fetch(region),
    cmap='fall',
    norm=norm,
    extent=extents
);
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='counts');
#Add the color bar on the right side

format_ticks(ax, rotate=False)
plt.tight_layout()
plt.savefig(output_path + "/" + output_prefix + ".png") #20kb_chr7_31M-43M_Micro_C_1_Bil.
#Format and save the figure
