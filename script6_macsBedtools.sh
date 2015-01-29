
#!/bin/bash

#This is a script to automate macs14 and bedtools analysis on the ChIP data. It assumes that you have sorted bam files for all the fastq data in the same directory where you have BAM files. sorted bams are faster.
#Input format is: path/Chip1.bam path/Control1.bam path/Chip2.bam path/Control2.bam
#MACS is run only on the first Chip dataset. If path/Chip1.bam is followed by path/Chip2.bam as opposed to path/Control1.bam, MACS is run without any input control on Chip1

#new "algorithm" that does not rely on defining the overlap in terms of read length....for this, finding the number of reads that have 50percent overlap with a peak is done in three steps:
#intersect -f 0.5 -a reads -b summitsextended -wa -wb
#cut -f 7,8,9,10
#intersect -a summits -b cutfile -c

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

#sum the reads like for the background, given that the read length is 140 and bin size is 600:
#bedtools intersect -a $summitspathname"_summitsextended.bed" -b $sortedbamprefix".bed" -f 0.117 -sorted > $pathname"ReadSumsTemp.bed"
#bedtools intersect -a $summitspathname"_summitsextended.bed" -b $sortedbamprefix".bed" -f 0.117 -sorted -c > $pathname"ReadSums.bed"

echo "Summing reads..."
#new algorithm (as mentioned at the start of the script)
#intersect -f 0.5 -a reads -b summitsextended -wa -wb
#cut -f 7,8,9,10
#intersect -a summits -b cutfile -c
bedtools intersect -a $sortedbamprefix".bed" -b $summitspathname"_summitsextended.bed" -f 0.5 -wa -wb > $pathname"ReadSummit.bed"
cut -f 7,8,9,10 $pathname"ReadSummit.bed" > $pathname"SummitOccurs.bed"
bedtools intersect -a $summitspathname"_summitsextended.bed" -b $pathname"SummitOccurs.bed" -c > $pathname"ReadSums.bed"

#should we fix that fact that BinnedData.bed counts reads spanning two bins in both the bins? We don't use BinnedData.bed quantitatively, so it's OK for the time being.
#create binned data file for all the reads output in the SAM/BAM file, this is not what we need
#but it is a good way to visualize all data
echo Creating binned data...
bedtools intersect -a /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w600.bed -b $sortedbam -c -sorted -bed > $pathname"BinnedData.bed"

#what we need is the background, i.e. we want to bin the data that does not come from peaks identified by MACS

echo "Creating BackgroundBins (coordinates without extended peaks)..."
#bedtools intersect -a /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w600.bed -b $summitspathname"_summitsextended.bed" -sorted -v -bed > $pathname"BackgroundBins.bed"
bedtools intersect -a /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w600.bed -b $summitspathname"_peaks.bed" -sorted -v -bed > $pathname"BackgroundBins.bed"
#in the second method above, I am calculating background from area not identified as signal by MACS.


#maybe create coordinates without extended peaks from PeakSplitter...MACS peaks may miss paired peaks closed by
#also have to install peak splitter on CB48
#bedtools intersect -a /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.w150.bed -b FirstSample_peaksControlnolambdaBW150Wiggle.subpeaksC250V97.bed -v -bed > FirstSample.w150NoMACSSubPeaks.bed
 
echo Generating background file...

#bedtools intersect -a $pathname"BackgroundBins.bed" -b $sortedbam -f 0.117 -sorted > $pathname"BackgroundTemp.bed"
#above is without c so I can do wc -l and get the total number of reads in background.
#bedtools intersect -a $pathname"BackgroundBins.bed" -b $sortedbam -f 0.117 -sorted -c > $pathname"Background.bed"

#new algorithm (as mentioned at the start of the script)
#intersect -f 0.5 -a reads -b summitsextended -wa -wb
#cut -f 7,8,9,10
#intersect -a summits -b cutfile -c
bedtools intersect -a $sortedbamprefix".bed" -b $pathname"BackgroundBins.bed" -f 0.5 -wa -wb > $pathname"ReadBackgroundBin.bed"
cut -f 7,8,9 $pathname"ReadBackgroundBin.bed" > $pathname"BackgroundBinOccurs.bed"
bedtools intersect -a $pathname"BackgroundBins.bed" -b $pathname"BackgroundBinOccurs.bed" -c > $pathname"Background.bed"

#notice that we could have moved creating BackgroundBins.bed to the sample with summits, since we need only one file with Background bins, but it doesn't matter at this
#point if we have multiple copies of BackgroundBin.bed, one in each sample/input directory

echo "--------"

done

echo Done!

#is there a way to intersect Yeast20min.w150NoMACSPeaks.bed with Yeast20minBinnedData.bed and keep the column with count reads? That will be faster
#than intersecting with the BAM file again.

#macs14 -t FirstSample.bam -c /mnt/bekraid/haz4z/YeastCLK/WTInput60sec/WTInput60sec.bam -f BAM --bw=150 -m 2,30 --slocal=1000 --llocal=5000 -g 1.21e+07 -n FirstSample -p 1.0e-02
#tried running with mfold 3,30, but got no paired peaks!! So now running with 2,30
#above step produces bed files: peaks, summits; and xls file: peaks
#can I run MACS on a sorted BAm file to save time? find out
#macs14 -t FirstSample.bam -c /mnt/bekraid/haz4z/YeastCLK/OEInput60sec/WTInput60sec.bam -f BAM --bw=150 -m 2,30 --nolambda -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -c /mnt/bekraid/haz4z/YeastCLK/OEInput60sec/WTInput60sec.bam -f BAM --bw=300 -m 2,30 --slocal=1000 --llocal=1500 -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -c /mnt/bekraid/haz4z/YeastCLK/OEInput60sec/OEInput60sec.bam -f BAM --bw=300 -m 2,30 -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -f BAM --bw=300 -m 2,30 -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -f BAM --bw=300 -m 2,30 --slocal=1000 --llocal=2000 -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -c /mnt/bekraid/haz4z/YeastCLK/OEInput60sec/OEInput60sec.bam -f BAM --bw=150 --nolambda -m 2,30 -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -f BAM --bw=150 --nolambda -m 2,30 -g 1.21e+07 -n FirstSample -p 1.0e-02
#macs14 -t FirstSample.bam -c /mnt/bekraid/haz4z/YeastCLK/OEInput60sec/OEInput60sec.bam -f BAM --bw=150 --nolambda -m 2,30 --wig -g 1.21e+07 -n FirstSample -p 1.0e-02
#I'm using the --wig option to create a wig file which I can then use with the peakSplitter program

#sum the reads in the summitextended regions...have to think what to do with summits that are <600bps apart. This is done in two steps:
#1- identify the reads that fit 50% within an extended peak
#echo Building Reads within peaks...
#bedtools intersect -a FirstSample_sorted.bed -b FirstSampleDupNoLambdaWithModel_summitsextended.bed -f 0.5 -bed > FirstSampleDupNoLambdaWithModelReadsWithinPeaks.bed
#2- Intersect the peaks with the above identified reads to count peaks that have these reads:
#bedtools intersect -a FirstSampleDupNoLambdaWithModel_summitsextended.bed -b FirstSampleDupNoLambdaWithModelReadsWithinPeaks.bed -c -bed -sorted > FirstSampleDupNoLambdaWithModelReadSums.bed 
#the two step process is beause I could not get bedtools to intersect reads with peaks, count the reads, but display the peaks with the count. Hence, I broke it up...
#This counting needs to be corrected for summits that are within 600bps of each other.

#extend the reads in our sorted bed file...
#bedtools slop -i FirstSample_sorted.bed -g /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3coords.genome -l 0 -r 178 -s > FirstSample_sorted_extended.bed

