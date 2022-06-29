import sys
from neoloop.visualize.core import *
import cooler
import pandas as pd

#Pass the arguments
arg_list = sys.argv[1:]
input_cool = arg_list[0] #should be an .cool file
input_assembly = arg_list[1]
input_neoloop = arg_list[2]
output_path = arg_list[3]
output_prefix = arg_list[4]

clr = cooler.Cooler(input_cool)
assembly = pd.read_csv(input_assembly, sep = "\t", header = None)
x = assembly.loc[3:,'0':'3'].to_string(header = False, index =  False, index_names = False).split('\n')
assembly = x[2]
pdf_name = output_path + '/' + output_prefix + '.pdf' #C42B_Micro_C_3.5kb.no_RWPE1.C2.pdf
vis = Triangle(clr, assembly, n_rows=5, figsize=(7, 5.2), track_partition=[5,0.8,0.8,0.8, 0.5], correct = 'weight')
vis.matrix_plot(vmin=0)
vis.plot_chromosome_bounds(linewidth=2.5)
vis.plot_loops(input_neoloop, face_color='none', marker_size=40, cluster=True)
#The below command is used to add another tracks for RNA-seq data. Can be compatiable with other NGS data like ChIP-seq
#vis.plot_signal('RNA', 'C42B_RNA_seq_rep3.bigwig', label_size=10, data_range_size=9, max_value=120, color='#FF0000')

#specify the genes that you want to show within the loops
vis.plot_genes(filter_=['ARID1A','DISP3'], fontsize=9)
vis.plot_genes(fontsize=9)
vis.plot_chromosome_bar(name_size=11, coord_size=4.8)
vis.outfig(pdf_name)



