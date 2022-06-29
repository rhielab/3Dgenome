


for i in *-5kb*.tsv;
do
prefixing=$i
second="1kb_10kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 10000 && $5-$2 <= 100000) print}' $i > $real_prefix
done


for i in *-5kb*.tsv;
do
prefixing=$i
second="100kb_200kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 100000 && $5-$2 <= 200000) print}' $i > $real_prefix
done

for i in *-5kb*.tsv;
do
prefixing=$i
second="200kb_400kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 200000 && $5-$2 <= 400000) print}' $i > $real_prefix
done

for i in *-5kb*.tsv;
do
prefixing=$i
second="400kb_600kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 400000 && $5-$2 <= 600000) print}' $i > $real_prefix
done

for i in *-5kb*.tsv;
do
prefixing=$i
second="600kb_800kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 600000 && $5-$2 <= 800000) print}' $i > $real_prefix
done


for i in *-5kb*.tsv;
do
prefixing=$i
second="800kb_1MB.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 800000 && $5-$2 <= 1000000) print}' $i > $real_prefix
done

for i in *-5kb*.tsv;
do
prefixing=$i
second="greater_1MB.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 > 1000000) print}' $i > $real_prefix
done


