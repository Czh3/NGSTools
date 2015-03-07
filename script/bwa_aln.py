#!/usr/bin/env python2

import sys
from string import Template
import os

ref = '/home/zhangc/MeDIP/zongle/rDNA_refSeq/rDNA.reform.fasta'



if len(sys.argv) < 3:
	sys.exit('usage: python bwa_aln.py fq1.gz fq2.gz')


command = '''
bwa aln $ref $fq1 > $fq1.sai
bwa aln $ref $fq2 > $fq2.sai
bwa sampe $ref $fq1.sai $fq2.sai $fq1 $fq2 | samtools view -h -q 1 -S -b - > $outbam
rm $fq1.sai $fq2.sai
'''


command = Template(command)
command = command.substitute(ref=ref, fq1=sys.argv[1], fq2=sys.argv[2], outbam=os.path.basename(sys.argv[1]).split('.')[0]+'.bam')

os.system(command)

