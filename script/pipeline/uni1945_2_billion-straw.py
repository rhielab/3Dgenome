
import straw
import numpy as np
import sys



chr_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
arg_list = sys.argv[1:]
resolution = int(arg_list[1])
prefix_name = arg_list[2]
file_name = arg_list[0]

#Generate sparse file format for TopDom

try:
	straw.straw('KR',file_name,'1','1','BP',resolution)
	test_chr = ""
except:
	test_chr= "chr"

for i in chr_list:
        try: 
                result = straw.straw('KR', file_name,(test_chr+i),(test_chr+i),'BP',resolution)
                chr_name = prefix_name + "_chr" + i + "_" + str(resolution) + 'bp.txt'
                with open(chr_name,"w") as f:
                        for i2 in range(len(result[0])):
                                print("{0}	{1}	{2}".format(result[0][i2], result[1][i2], result[2][i2]),file=f)
        except:
               if i == "Y":
                       print("Extracting chromosome Y failed, continuing")
               else:
                       raise # any other chromosome should exist
