#/bin/sh

#usage MACs.sh sample.bam control.bam name

macs14 -t $1	\
	-c $2	\
	-n $3	\
	-w --space=50

