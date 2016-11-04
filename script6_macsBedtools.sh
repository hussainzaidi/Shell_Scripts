
#!/bin/bash

#This is a script to automate macs14 and bedtools analysis on the ChIP data. It assumes that you have sorted bam files for all the fastq data in the same directory where you have BAM files. sorted bams are faster.
#Input format is: path/Chip1.bam path/Control1.bam path/Chip2.bam path/Control2.bam
#MACS is run only on the first Chip dataset. If path/Chip1.bam is followed by path/Chip2.bam as opposed to path/Control1.bam, MACS is run without any input control on Chip1

source /mnt/bekraid/haz4z/source_me.sh


args=("$@")

counter=0
MACSrun=0

echo " "

for var in "$@"
do

 counter=$((counter+1))
 if [[ $var != *.bam ]]
 then
  echo "BAM file not provided in $var"
 else
  pathname="${var%.bam}"
 fi
 samfile=$pathname".sam"
 bamfile=$pathname".bam"
 sortedbamprefix=$pathname"_sorted"
 sortedbam=$sortedbamprefix".bam"

if [ $MACSrun == 0 ]
then
 if [[ $pathname == *Chip* ]]
 then
 MACSrun=1
 summitspathname=$pathname
 if [ $counter == "$#" ]
 then
  echo Only one Chip file passed without control...
  echo Running MACS...
  macs14 -t $bamfile -f BAM --bw=223 --keep-dup=all --nolambda -m 2,30 -g 1.21e+07 -n $pathname -p 1.0e-02
 else
 varNext=${args[${counter}]}
  if [[ $varNext == *Input* ]]
  then
   pathnameInput=$varNext
   echo "Running MACS on $bamfile with $pathnameInput as control "
   macs14 -t $bamfile -c $pathnameInput -f BAM --bw=223 --keep-dup=all --nolambda -m 2,30 -g 1.21e+07 -n $pathname -p 1.0e-02 
  else
  echo Running MACS on sample $pathname without any input
  macs14 -t $bamfile -f BAM --bw=223 --keep-dup=all --nolambda -m 2,30 -g 1.21e+07 -n $pathname -p 1.0e-02 
  fi
 fi

#extend the summits to +/-300 bps:
echo "Running slop..."
bedtools slop -i $summitspathname"_summits.bed" -g /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3coords.genome -b 299 > $summitspathname"_summitsextended.bed"
 fi
fi

#samtools faidx sacCer3.fa
#above step creates an index for the reference genome

#copy the length of the chromosomes (column 2) from index file (.fai) into a .genome file for bedtools, and then go to the next step
#bedtools makewindows -g /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3coords.genome -w 600 > /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w600.bed

#create a bed file from sorted BAM with strand information
echo "Creating a sorted BED file for the sample..."
bedtools bamtobed -i $sortedbam > $sortedbamprefix".bed"


echo "Summing reads..."
bedtools intersect -a $sortedbamprefix".bed" -b $summitspathname"_summitsextended.bed" -f 0.5 -wa -wb > $pathname"ReadSummit.bed"
cut -f 7,8,9,10 $pathname"ReadSummit.bed" > $pathname"SummitOccurs.bed"
bedtools intersect -a $summitspathname"_summitsextended.bed" -b $pathname"SummitOccurs.bed" -c > $pathname"ReadSums.bed"

echo Creating binned data...
bedtools intersect -a /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w600.bed -b $sortedbam -c -sorted -bed > $pathname"BinnedData.bed"

#what we need is the background, i.e. we want to bin the data that does not come from peaks identified by MACS

echo "Creating BackgroundBins (coordinates without extended peaks)..."
bedtools intersect -a /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w600.bed -b $summitspathname"_peaks.bed" -sorted -v -bed > $pathname"BackgroundBins.bed"
 
echo Generating background file...

bedtools intersect -a $sortedbamprefix".bed" -b $pathname"BackgroundBins.bed" -f 0.5 -wa -wb > $pathname"ReadBackgroundBin.bed"
cut -f 7,8,9 $pathname"ReadBackgroundBin.bed" > $pathname"BackgroundBinOccurs.bed"
bedtools intersect -a $pathname"BackgroundBins.bed" -b $pathname"BackgroundBinOccurs.bed" -c > $pathname"Background.bed"

echo "--------"

done

echo Done!
