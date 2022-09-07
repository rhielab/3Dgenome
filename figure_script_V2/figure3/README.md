## Figure3 Chromatin loops near structural variants.

### 3A. Numbers of inter and intra-chromosomal structural variants identified from Hi-C and Micro-C data are shown. 3B. Numbers of each category of structural variants identified from Hi-C and Micro-C data are shown. 3C. Numbers of loops identified around the structural variants from Hi-C and Micro-C data are shown at 5kb and 10kb resolutions.

We used NeoLoopFinder (Wang et al., 2021) to identify loops near the structural variants using Hi-C and Micro-C data. As the figures are showing the stats number from NeoLoopFinder, here we explain our code and relevant scripts to run NeoLoopFinder. To start with, we first runed hic_breakfinder (https://github.com/dixonlab/hic_breakfinder) to identify the structural variants. 

The necessary input files for hic_breakfinder are ```.bam``` as well as the intra and inter chromosomal expectation file. Here we used ```inter_expect_1Mb.hg38.txt``` and ```intra_expect_100kb.hg38.txt```, which are provided by hic_breakfinder and can be found at https://salkinstitute.box.com/s/m8oyv2ypf8o3kcdsybzcmrpg032xnrgx.

The example code to run hic_breakfinder:
```
hic_breakfinder \
--bam-file input.bam \ 
--exp-file-inter inter_expect_1Mb.hg38.txt \
--exp-file-intra intra_expect_100kb.hg38.txt \
--name OutputPrefix
```

We then used [prepare-SV-breakpoints.py](https://github.com/XiaoTaoWang/NeoLoopFinder/blob/master/scripts/prepare-SV-breakpoints.py) from [NeoLoopFinder GitHub](https://github.com/XiaoTaoWang/NeoLoopFinder) to convert the ```breaks.txt``` to ```sv.txt``` which is used as the input for NeoLoopFinder.
```
python prepare-SV-breakpoints.py \
OutputPrefix.breaks.txt \
OutputPrefix.sv.txt
```

And then we can run NeoLoopFinder. The necessary input file include ```.cool``` and the newly generated ```sv.txt```. 

Based on the instruction, need to first run ```cooler balance```.
```
cooler balance -p 1 EXAMPLE.cool
```
Noted that it is also required to run ```calculate-cnv```, ```segment-cnv``` and ```correct-cnv```. However, these functions didn't work with our Micro-C data so we skip these steps and just run ```assemble-complexSVs```:
```
assemble-complexSVs \
-O OutputPrefix \
-B OutputPrefix.sv.txt \
-H EXAMPLE.cool \
--balance-type ICE
```
This step will create a ```assemblies.txt``` which contains the structural variant events and also serves as the input for ```neoloop-caller```:
```
neoloop-caller \
-O OutputPrefix.neo-loops.txt \
-H EXAMPLE.cool \
--assembly OutputPrefix.assemblies.txt \
--no-clustering \
--prob 0.95 \
--balance-type ICE
```
The loops identified near structural variants are save in ```neo-loops.txt```. For more detailed instruction about NeoLoopFinder, please refer to https://github.com/XiaoTaoWang/NeoLoopFinder.

### 3D. Numbers of neoloops (loops newly gained due to the structural variants) that are shared (between any datasets) or unique among Hi-C 1 billion, Micro-C 1 billion, 2 billion and 3 billion data are shown.

Here is the loop-comparison analysis, we used the same script as desribed in [figure_script/figure2/2D](https://github.com/rhielab/3Dgenome/tree/main/figure_script/figure2/2D). Please refer for more details.

### 3E. An example heatmap of Micro-C data near the ARID1A gene that includes inversion structural variant is shown on the top. Under the heatmap, RNA-seq and RefSeq gene tracks are shown. Example neoloops newly gained due to the structural variants are circled in blue.

We used ```plotting_NEOheatmap.py``` to generate the heatmap together with multipe tracks of data. The necessary input to run this script include ```.cool``` with additional ```assemblies.txt``` and ```neo-loops.txt```, which are generate from NeoLoopFinder (Wang et al., 2021). Please refert to [
/figure_script/figure3/3A.B.C](https://github.com/rhielab/3Dgenome/tree/main/figure_script/figure3/3A.B.C) for more details on how to generate these two files.

Below is an example code:
```
python plotting_NEOheatmap.py \
EXAMPLE.cool \
FromNeoLoopFinder.assemblies.txt \
FromNeoLoopFinder.neo-loops.txt \
./output/ \
Outputprefix
```

PLEASE NOTED that if you want to add tracks based on your need, the additional ```.bigwig``` files are needed and in this case you need to modify ```plotting_NEOheatmap.py``` by yourself (see line 23 and 24). Please refer to https://github.com/XiaoTaoWang/NeoLoopFinder for more details.
