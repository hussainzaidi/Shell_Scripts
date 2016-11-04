
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
