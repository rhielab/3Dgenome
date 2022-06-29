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

bp_formatter = EngFormatter('b')


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
clr_1kb_hic=  cooler.Cooler('../uni1945_2.5kb_1_billion.filt.mcool::resolutions/5000')

region = 'chr7:36,000,000-38,000,000'
start, end = 36_000_000, 38_000_000
extents = (start, end, end, start)
norm = LogNorm(vmin=1, vmax=100)
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)


f, ax = plt.subplots(
    figsize=(13, 10),
    nrows=1,
    ncols=1,
    sharex=True,
    sharey=True
)


im = ax.matshow(
    clr_1kb_hic.matrix(balance=False).fetch(region),
    cmap='fall',
    norm=norm,
    extent=extents
);
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='counts');

format_ticks(ax, rotate=False)
plt.tight_layout()
plt.savefig("5kb_chr7_36M-38M_Micro_C_1_Bil_Micro_C.png")

