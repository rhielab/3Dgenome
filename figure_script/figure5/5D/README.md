## Figure5 Promoter capture Micro-C data analysis.

### 5D. Fractions of active promoters that intersect with the loop anchors from Micro-C 1 billion, 2 billion, 3 billion data or promoter capture Micro-C data are shown (left). A fraction of active promoters that intersect with loop anchors from any datasets is shown in grey (In loop) while the one not in loop is shown in orange (Not in loop) (right).

Here we provide an Rscript ```checkCGI_promoter.R``` to run the analysis. Please noted that the panel itself is not showing any analysis relevant with CGI promoter but just promoter. You can still use this script to check how many promoters are involved in your loop data; specify any bed file for CGI reference and just ignore the CGI results.

However, if you are interested to check CGI promoter, then you need to specify the correct CpG Isaland reference. Please refer to [UCSC genome browser](https://genome.ucsc.edu) for it. 

To run this script:
```
source('checkCGI_promoter.R')
myloopdata <- read.table(file = '/PATH/TO/LoopData')
checkCGI_promoter(loop_data= myloopdata, 
output_path='./output', 
prefix = 'outputprefix')
```
In the output, you will see the promoter-relevant files ending with ```Promoter_inLoop.bed``` or ```Promoter_nonLoop.bed```. If you are using the correct CGI reference, then you can also check files containing with ```CGIPromoter```. You can check ```summary.txt``` for all the necessary stats.
