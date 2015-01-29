
#!/bin/bash

#This is initial script to automate bowtie2-build, bowtie2, samtools view, samtools sort, samtools index on the sacCer3 genome and ChIP data. You can feed in FASTQ, Sam, or BAM files, and get 
#sortedBAM files as output


source /mnt/bekraid/haz4z/source_me.sh

FastqToSam ()
{
 echo "Running bowtie with input $var and output $samfile..."
 bowtie2 -x /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3 -U $var -S $samfile -p 3 --time
#echo "hello"
}

SamToBam ()
{
#keep the header: -h
#skip alignments with MAPQ smaller than INT: -q INT
#skip alignments with INT in FLAG field: -F 4 (-F 4 => skip unmpapped read)

 echo "Running samtools view with input $samfile and output $bamfile..."
 samtools view -b -S -h -q 20 -F 4 -o $bamfile $samfile 
}

BamToSortedBam ()
{
 echo "Running samtools sort with input $bamfile and output prefix $sortedbamprefix..."
 samtools sort $bamfile $sortedbamprefix
 echo "Running samtools index with input $sortedbam..."
 samtools index $sortedbam $sortedbamprefix".bai"
}


#bowtie2-build /mnt/bekraid/genomes/S_cerevisiae_2011/sacCer3.fa sacCer3
 
#echo ${args[0]}
#for var in $args
#echo "$@"
#echo "$@"

for var in "$@"
do
 if [[ $var == *.fastq ]]
 then
 pathname="${var%.fastq}"
 samfile=$pathname".sam"
 bamfile=$pathname".bam"
 sortedbamprefix=$pathname"_sorted"
 sortedbam=$sortedbamprefix".bam"
 FastqToSam
 SamToBam
 BamToSortedBam

 elif [[ $var == *.sam ]]
 then
 pathname="${var%.sam}"
 samfile=$pathname".sam"
 bamfile=$pathname".bam"
 sortedbamprefix=$pathname"_sorted"
 sortedbam=$sortedbamprefix".bam" 
 SamToBam
 BamToSortedBam

 elif [[ $var == *.bam && $var != *_sorted* ]]
 then
 pathname="${var%.bam}"
 samfile=$pathname".sam" 
 bamfile=$pathname".bam" #same as $var
 sortedbamprefix=$pathname"_sorted"
 sortedbam=$sortedbamprefix".bam" 
 BamToSortedBam

# elif [[ $var == *_sorted.bam ]]
# then
# pathname="${var%_sorted.bam}"
# samfile=$pathname".sam" 
# bamfile=$pathname".bam" 
# sortedbamprefix=$pathname"_sorted"
# sortedbam=$sortedbamprefix".bam" 
fi
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

