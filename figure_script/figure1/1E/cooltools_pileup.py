import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cooltools
import cooler
import bioframe

import bbi
from scipy.stats import linregress
import cooltools.lib.plotting
from ipywidgets import interact
from matplotlib.gridspec import GridSpec

#Process argument, and converthem into resolution ,prefix name and file name.
arg_list = sys.argv[1:]
input_file = arg_list[0] #should be an .mcool file
input_loop = arg_list[1] #should be loop.bedpe file identified from other loop-recoginition tools
output_path = arg_list[2]
output_prefix = arg_list[3]

#Load cooler data. With mcool, the data has multiple resolution 
clr = cooler.Cooler(input_file + '::/resolutions/1000')

#Load hg38 view data 
hg38_arms = pd.read_csv("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure1/cooltools_pileup/input/hg38_view_df.bed", sep = "\t")

#Obtain resolution information from cooler data
resolution = clr.binsize

#Index the arm. For this purprose, I needed to use cooler file that filtered out chr that are not from chr1-22 and X/Y, since the program would later complain about index and chr names in cooler not matching. This is the reason I used 1 billlion Micro-C data instead of 3 billion data, since filtered cooler data took more than two weeks to be made with 3 billion data for some reason.
hg38_arms = hg38_arms.set_index("chrom").loc[clr.chromnames].reset_index()


#Road loop files in bedpe format
paired_sites = bioframe.read_table(input_loop, schema = 'bedpe')

flank = 250 # Length of flank to one side from the boundary, in basepairs
expected = cooltools.expected_cis(clr, view_df=hg38_arms, nproc=2, chunksize=1_000_000)

stack = cooltools.pileup(clr, paired_sites, view_df=hg38_arms, expected_df=expected, flank=100_000)

mtx = np.nanmean(stack, axis=2)

plt.imshow(
    np.log2(mtx),
    vmax = 1,
    vmin = -1,
    cmap='coolwarm')

plt.colorbar(label = 'log2 mean obs/exp')
ticks_pixels = np.linspace(0, flank*2//resolution,5)
ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
plt.xticks(ticks_pixels, ticks_kbp)
plt.yticks(ticks_pixels, ticks_kbp)
plt.xlabel('relative position, kbp')
plt.ylabel('relative position, kbp')

plt.savefig(output_path + '/' + output_prefix +".png")
