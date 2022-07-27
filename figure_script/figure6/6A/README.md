## Figure6 Nucleosome phasing and DNA methylation levels around regulatory elements involved in loops.

### 6A. Average Micro-C signals around active promoters, enhancers, insulators, and NDRs without features that are in loop vs those that are not in loop.

The signal plots are generated using DeepTools (Ramírez et al, 2014). The necessary input file is ```.bigwig``` as well as the ```.bed``` files for different regulatory elements (e.g promoter, enhancer, insulator, NDR). Please refer to https://deeptools.readthedocs.io/en/develop/ for more information. Below is our example code ro run DeepTools:

```
computeMatrix reference-point --referencePoint center \
-b 1000 -a 1000 \
-S input.bigwig \
-R RegulatoryElement.bed \
-p 16 \
-o output.RPKM.gz
```
The output ```RPKM.gz``` is necessary for ```plotHeatmap``` function
```
plotHeatmap \
-m output.RPKM.gz \
-o output.pdf
```
and also used to perform t-test. Here we provide an Rscript ```t_test_4deeptools.R``` to perform t-test at the center of the signal between in-loop and not in-loop data. To use:
```
Rscript t_test_4deeptools.R /PATH/TO/inloop.RPKM.gz /PATH/TO/notinloop.RPKM.gz
```
It will print out the t-test result.
