## Figure2 Comparison of Micro-C data in different read depth sequencing.

### 2B. Numbers of loops identified by Mustache at different resolutions from Hi-C 1 billion, Micro-C 1 billion, Micro-C 2 billion, and Micro-C 3 billion data are shown.

We used the loop data genrated by Mustache (Roayaei Ardakany et al., 2020) for different libraries (Hi-C 1 billion; Micro-C 1,2,3 billion data) at 1kb, 2kb, 5kb, 10kb, 25kb, 50kb respectively. The plot is showing the number of loops we get at different different resolutions among libraries.

If you are interested to use Mustache, please refer to https://github.com/ay-lab/mustache for more detailed instructions. The necessary input is ```.hic``` format. Here we provide our example code to run Mustache.

Need to specify the specific chromosome when running Mustache. If the interested resolution is less than 5kb (e.g. 2kb), we used specialized settings.
```
chr_prefix="chr"
resolution=2kb
for chr_number in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
mustache -f EXAMPLE.hic \
-ch $chr_prefix$chr_number \
-r $resolution \
-st 0.7 \
-pt 0.1 \
-o output.tsv
done
```

If the interested resolution is equal or more than 5kb, we used the deafault setting.
```
chr_prefix="chr"
resolution=5kb
for chr_number in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
mustache -f EXAMPLE.hic \
-ch $chr_prefix$chr_number \
-r $resolution \
-o output.tsv
done
```
