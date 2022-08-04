## Figure2 Comparison of Micro-C data in different read depth sequencing.

### 2C. Fractions of loops that have different lengths (distances) found from Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data are shown. 

Here we provide a script ```category_loop_distance.sh``` which is used to classify a specific loop data into different loop distance categoires. (<200kb, 200kb-400kb, 400kb-600kb, 600kb-800kb, 800kb-1Mb, > 1Mb) and summarize the fraction of each category. To use this script:

```
bash category_loop_distance.sh \
/PATH/TO/LOOPDATA
```

The deafault output path is your current working directory. You can also specify the output path:

```
bash category_loop_distance.sh \
/PATH/TO/LOOPDATA \
OUTPUT_PATH
```

Outputs include a bunch of ```.bedpe``` files for loops falling into specific category, along with a ```summary.txt``` where the number and fraction of loops detected in different categoires are shown.
