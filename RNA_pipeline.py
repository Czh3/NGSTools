#-*- coding: UTF-8 -*-
import sys
#sys.path.append('/home/zhangc/bin/git/TEST')

import NGSTools
import argparse
import os

#argparse arguments
parser = argparse.ArgumentParser(description='A pipeline of RNA_seq data analysis. <zhangchao3@hotmail.com>')
parser.add_argument('-s', '--sampleList', help='sample list for RNA samples information:\nOne line per sample:\nsampleName\tsampleCondition\tfastq1Path\tfastq2Path', required=True)
parser.add_argument('-o', '--outDir', help='output dir', default='out')
parser.add_argument('-d', '--dataType', help='fastq data type:raw data or clean data', choices=['raw', 'clean'], default='clean')
parser.add_argument('-t', '--tools', help='tools using to call DEGs: cufflinks or DESeq', choices=['cufflinks', 'DESeq'], default='cufflinks')
parser.add_argument('-c', '--config', help='the config file of NGSTools package.', default='~/.NGSTools.cfg')
parser.add_argument('--debug', help='debug mode', default=False)
args = parser.parse_args()


#init
__VERSION__ = 'V0.1'

if not os.path.isfile(args.sampleList):
	sys.exit('ArgumentError:\tSampleList must be a file.')

args.outDir = os.path.abspath(args.outDir)
if not os.path.exists(args.outDir):
	os.mkdir(args.outDir)

if args.debug:
	_run = False
else:
	_run = True

condition = {}
transcripts = []

##
gtf = ''
genome = ''

for line in open(args.sampleList):
	if line.startswith('#'):
		continue

	cols = line.strip().split('\t')
	sample = {
		'name' : cols[0],
		'condition' :	cols[1],
		'fq1' : cols[2],
		'fq2' : cols[3],
		'bam' : ''
	}


	########################## 0. init #########################
	#__init__(self, sampleName, outdir, fq1, fq2='', quanlityBase='32', cfgfile='~/.NGSTools.cfg'):
	mySample = NGSTools.NGSTools(sample['name'], args.outDir, sample['fq1'], sample['fq2'], cfgfile=os.path.abspath(args.config))
	gtf = mySample.gtf
	genome = mySample.genome
	
	#################### 1. Quality Control ####################

	###### 1.1 cut adapter ######
	if args.dataType == 'raw':
		mySample.cutadapter(run = _run)
	else:
		pass

	##### 1.2 fastqc #####
	mySample.fastqc(run = _run)

	######################## 2. Mapping ########################
	
	sample['bam'] = mySample.tophat2(run = _run)


	if condition.has_key(sample['condition']):
		condition[sample['condition']][sample['name']] = sample['bam']
	else:
		condition[sample['condition']] = {sample['name'] : sample['bam']}


	######################## 3. DEGs calling ########################

	##### 3.1 cufflinks #####
	if args.tools == 'cufflinks':
		cuffdir = os.path.join(args.outDir, 'cufflinks')
		if not os.path.exists(cuffdir):
			os.mkdir(cuffdir)

	# cufflinks #
	command = 'cufflinks -p 4 -g %s -o %s %s' % (mySample.gtf, os.path.join(cuffdir, sample['condition']+'_'+sample['name']), sample['bam'])
	NGSTools.writeCommands(command, cuffdir+'/cufflinks_%s.sh' % sample['name'], _run)

	transcripts.append(os.path.join(cuffdir, sample['condition']+'_'+sample['name'], 'transcripts.gtf'))


if args.tools == 'cufflinks':
	# cuffmerge #
	with open(cuffdir+'/assemblies.txt', 'w') as writer:
		writer.write('\n'.join(transcripts))

	command = 'cuffmerge -o %s -g %s -s %s -p 10 %s' % (cuffdir+'/merged_asm', gtf, genome, cuffdir+'/assemblies.txt')
	NGSTools.writeCommands(command, cuffdir+'/cuffmerge.sh', _run)

	# cuffdiff #
	cuffCondition = []
	cuffBam = {}
	for c in condition.keys():
		cuffCondition.append(c)
		cuffBam[c] = ','.join(condition[c].values())

	if len(cuffCondition) != 2:
		sys.exit('Error: condition')

	command = 'cuffdiff -o %s -b %s -p 10 -L %s -u %s %s %s' % (cuffdir+'/cuffdiff', genome, cuffCondition[0]+','+cuffCondition[1], cuffdir+'/merged_asm/merged.gtf', cuffBam[cuffCondition[0]], cuffBam[cuffCondition[1]])
	NGSTools.writeCommands(command, cuffdir+'/cuffdiff.sh', _run)






