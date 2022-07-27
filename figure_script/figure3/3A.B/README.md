## Figure3 Chromatin loops near structural variants.

### 2(A) Numbers of inter and intra-chromosomal structural variants identified from Hi-C and Micro-C data are shown. 2(B) Numbers of each category of structural variants identified from Hi-C and Micro-C data are shown.

We used NeoLoopFinder (Wang et al., 2021) to identify loops near the structural variants using Hi-C and Micro-C data. As the figure is showing the stats number from NeoLoopFinder, here we show our code and relevant scripts to run this analysis. To start with, we first runed hic_breakfinder (https://github.com/dixonlab/hic_breakfinder) to identify the structural variants. 

The necessary input files for hic_breakfinder are ```.bam``` as well as the intra and inter chromosomal expectation file which can be found at https://salkinstitute.box.com/s/m8oyv2ypf8o3kcdsybzcmrpg032xnrgx. Here we used ```inter_expect_1Mb.hg38.txt``` and ```intra_expect_100kb.hg38.txt```.

The example code to run hic_breakfinder:
```
hic_breakfinder \
--bam-file input.bam \ 
--exp-file-inter inter_expect_1Mb.hg38.txt \
--exp-file-intra intra_expect_100kb.hg38.txt \
--name OutputPrefix
```

We then used ```convert_break_sv.py``` to convert the ```breaks.txt``` to ```sv.txt``` which is used as the input for NeoLoopFinder.
```
python convert_break_sv.py \
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
