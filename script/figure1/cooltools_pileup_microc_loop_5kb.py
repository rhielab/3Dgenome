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

clr = cooler.Cooler('../uni1945_1kb_1_billion.filt.mcool::/resolutions/1000')


hg38_arms = pd.read_csv("../hg38_view_df.bed", sep = "\t")

resolution = clr.binsize


hg38_arms = hg38_arms.set_index("chrom").loc[clr.chromnames].reset_index()


paired_sites = bioframe.read_table('1_Billion_Mustache-5kb-all-loop.bedpe', schema = 'bedpe')

flank = 500 # Length of flank to one side from the boundary, in basepairs
expected = cooltools.expected_cis(clr, view_df=hg38_arms, nproc=2, chunksize=1_000_000)

stack = cooltools.pileup(clr, paired_sites, view_df=hg38_arms, expected_df=expected, flank=10_000)

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

plt.savefig("Micro_C_5kb_Mustache_Loop_Flank_10kb_500.png")




