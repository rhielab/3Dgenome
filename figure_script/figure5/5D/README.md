## Figure5 Promoter capture Micro-C data analysis.

### 5D. Fractions of active promoters that intersect with the loop anchors from Micro-C 1 billion, 2 billion, 3 billion data or promoter capture Micro-C data are shown (left). A fraction of active promoters that intersect with loop anchors from any datasets is shown in grey (In loop) while the one not in loop is shown in orange (Not in loop) (right).

Here we provide an Rscript ```check_REGelement.R``` to find out the regulatory elements involved in loops. To run this script, the loop data under specific resolution as well as the regulatory element data (in ```.bed``` format) is necessary. Check below code for example usage to identify enhancers involved in loop under 5kb resolution.
```
source('check_REGelement.R')
myloopdata <- read.table(file = '/PATH/TO/LoopData')
#E.g you want to find the enhancers involved in loops.
myEnhdata <- read.table(file = '/PATH/TO/EnhancerData')

checkCGI_promoter(
loop_data = myloopdata, 
reg_data = myEnhdata,
reg_type = 'Enhancer', #To clarify, this option is just a character to help note the regulatory element you used.
resolution = 5000,
output_path='./output', 
prefix = 'outputprefix')
```
In the output, you will see two files ending with ```inLoop.bed``` or ```nonLoop.bed```. You can check ```summary.txt``` for the necessary stats.
