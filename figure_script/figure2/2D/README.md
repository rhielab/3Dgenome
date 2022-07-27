## Figure2 Comparison of Micro-C data in different read depth sequencing.

### 2B. Numbers of loops shared (between any datasets) or unique among Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data.

Here we provide an in-house Rscript ```Loop_comparison.R``` to run loop intersection analysis between two libraries, using a list of loop data as the input. 

The output files in the output folder include:
- All the original loop data added with comparison 0/1 matrix.
- For each loop data, output a subforder containning the shared/unique loop data
- Summary table for the loop comparison.

See below code for usage instruction:

Load the function. To run it, make sure you have ```fuzzyjoin``` and ```dplyr``` libraries in your R environment:
```
source('Loop_comparison_Tim.R')
```

Then load the loop data. You should have at least two loop data to run this analysis. But let's use four as example here:
```
loop_1 = read.table(file = 'loop1.tsv')
loop_2 = read.table(file = 'loop2.tsv')
loop_3 = read.table(file = 'loop3.tsv')
loop_4 = read.table(file = 'loop4.tsv')
```

Create a NAMED list in your prefered order. Here, let's say the priority to assign overlap is loop_4 > loop_3 > loop_2 > loop_1, then:
```
loop_list <- list(
LOOP4 = loop_4,
LOOP3 = loop_3,
LOOP2 = loop_2,
LOOP1 = loop_1
)
```

Run the function using the list just created. The overlap will be calculated based on the priority same as the order in the list in deafault.
```
loop_comparison(
loop_list = loop_list,
output_path= './output'
)
```

Or suddendly if you want a different priority. For example, loop_1 > loop_2 > loop_3 > loop_4, then:
```
loop_comparison(
loop_list = loop_list,
priority = c(1,2,3,4)
output_path= './output2'
)
```
