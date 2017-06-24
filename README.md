# Shell_Scripts
Scripts for automating running bowtie, samtools, MACS, and bedtools on ChIP-Seq fastq files. 
These scripts came about after a few rounds of iteration on what works and what does not work when identifying peak regions in ChIP-Seq
data. The scripts will take the data from fastq to all the way through to bed files with peak regions identified and extended 300 bps
in each direction (a different extension size can be set in the scripts).
