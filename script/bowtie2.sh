#!/usr/bin/sh

#usage:
#args1: read1.fq.gz
#args2: read2.fq.gz
#args3: output_prefix

bowtie2 --phred33 -q --end-to-end -p 20 \
	-x ../rDNA_refSeq/rDNA	\
	-1 $1	\
	-2 $2	\
	-S $3".sam"
# use mapping reads
samtools view -bS -h -q 30 $3".sam" > $3".bam"
samtools sort -m 1G -@ 5 $3".bam" $3".sort"
samtools index $3".bam" 
