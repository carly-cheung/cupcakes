# bedtools for sequencing analysis

## Starting and logging in to a Spec Ops instance:

```
aegea launch --iam-role S3fromEC2 --ami-tags Name=czbiohub-specops -t t2.micro <your-image-name>
aegea ssh ubuntu@<your-image-name>
cd /mnt/data/

```
DON'T FORGET TO TERMINATE YOUR INSTANCE WHEN YOU'RE DONE!



## Installing bedtools:
```
sudo apt install bedtools
```

## Conveting a sam file to a bam file:

```
samtools view -bhS <input.sam> > <output.bam>

#or#

for i in *.sam; do samtools view -bhS $i > ${i:0:-3}bam; done;

mkdir sam/
mv *.sam sam/
```
## Converting a bam file to a bed file:
```
bedtools bamtobed -i <input.bam> > <output.bed>

#or#

for i in *.bam; do bedtools bamtobed -i $i > ${i:0:-2}ed; done;

mkdir bam/
mv *.bam bam/
```

## Counting reads that intersect with a feature:
```
bedtools intersect [options] -a <a_file.bed> -b <b_file.bed>

bedtools intersect -f 1 -c -a labSnpRanges.bed -b 100k_Q01_FLASH.bed > intersect_labSnpRanges_100k_Q01_FLASH.bed

#or#

bedtools intersect -f 0.1 -c -a labSnpRanges.bed -b 100k_Q01_FLASH.bed > intersect_labSnpRanges_100k_Q01_FLASH.bed

mkdir bed/
mv *.bed bed/
```

## Looking at your results:

```
more intersect_labSnpRanges_100k_Q01_FLASH.bed

grep DR1 intersect_labSnpRanges_100k_Q01_FLASH.bed

grep DR1 intersect_labSnpRanges_100k_Q0*

```
# Terminate your instance when you're done:

```
exit

aegea terminate <your-image-name>
```