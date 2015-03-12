#-*- coding: UTF-8 -*-
import sys
sys.path.append('/home/zhangc/bin/git/TEST')

import NGSTools
import argparse
import os

#argparse arguments
parser = argparse.ArgumentParser(description='A pipeline of RNA_seq data analysis. <zhangchao3@hotmail.com>')
parser.add_argument('-s', '--sampleList', help='sample list for RNA samples information:\nOne line per sample:\nsampleName\tsampleCondition\tfastq1Path\tfastq2Path', required=True)
parser.add_argument('-o', '--outDir', help='output dir', default='out')
parser.add_argument('-d', '--dataType', help='fastq data type:raw data or clean data. if (clean data): not run cutadapter', choices=['raw', 'clean'], default='clean')
parser.add_argument('-a', '--analysis', help='analysis of the pipeline to do.[1:QC, 2:Mapping, 3:Cufflinks, 4:DESeq2, 5:DEXSeq, 6:GATK]', default='1,2,4')
parser.add_argument('-c', '--config', help='the config file of NGSTools package.', default='~/.NGSTools.cfg')
parser.add_argument('--debug', help='debug mode', default=False)
args = parser.parse_args()


# init
__VERSION__ = 'V0.1'

if not os.path.isfile(args.sampleList):
	sys.exit('ArgumentError:\tSampleList must be a file.')

args.outDir = os.path.abspath(args.outDir)
if not os.path.exists(args.outDir):
	os.mkdir(args.outDir)

analy = args.analysis.split(',')
analy = [int(i) for i in analy]

if max(analy) > 2 and 2 not in analy:
	sys.exit('analysis error:\tMust mapping the reads in advance')

QC = Mapping = Cufflinks = DESeq2 = DEXSeq = GATK = False

if 1 in analy:
	QC = True
if 2 in analy:
	Mapping = True
if 3 in analy:
	Cufflinks = True
if 4 in analy:
	DESeq2 = True
if 5 in analy:
	DEXSeq = True
if 6 in analy:
	GATK = True

if args.debug:
	_run = False
else:
	_run = True

condition = {}
transcripts = []
countsFiles = {}
finalBam = {}

##
cfg = NGSTools.getConfig(os.path.abspath(args.config))


for line in open(args.sampleList):
	if line.startswith('#'):
		continue

	if line == '\n':
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
	

	if QC:
		#################### 1. Quality Control ####################

		###### 1.1 cut adapter ######
		if args.dataType == 'raw':
			mySample.cutadapter(run = _run)
		else:
			pass

		##### 1.2 fastqc #####
		mySample.fastqc(run = _run)
	
	
	if Mapping:
		######################## 2. Mapping ########################
		
		sample['bam'] = mySample.tophat2(run = _run)


		if condition.has_key(sample['condition']):
			condition[sample['condition']][sample['name']] = sample['bam']
		else:
			condition[sample['condition']] = {sample['name'] : sample['bam']}


		if GATK:

			# remove duplicates
			mySample.picard_rmdup(run = _run)

			# splitN
			mySample.splitN(run = _run)

			# realign
			mySample.realn(run = _run)

			# recal
			recalBam = mySample.recal(run = _run)

			finalBam[recalBam] = sample['condition']


	########################  DEGs calling preparation ########################

	if Cufflinks:
		##### 3. cufflinks #####
		cuffdir = os.path.join(args.outDir, 'cufflinks')
		if not os.path.exists(cuffdir):
			os.mkdir(cuffdir)

		# cufflinks #
		command = 'cufflinks -p 4 -o %s %s' % (os.path.join(cuffdir, sample['condition']+'_'+sample['name']), sample['bam'])
		NGSTools.writeCommands(command, cuffdir+'/cufflinks_%s.sh' % sample['name'], _run)

		transcripts.append(os.path.join(cuffdir, sample['condition']+'_'+sample['name'], 'transcripts.gtf'))

	if DESeq2:
		##### 4. DESeq2 #####
		count = mySample.HTSeq_count(run = _run)
		countsFiles[count] = sample['condition']


########################  DEGs calling ########################
if Cufflinks:
	# cuffmerge #
	with open(cuffdir+'/assemblies.txt', 'w') as writer:
		writer.write('\n'.join(transcripts))

	command = 'cuffmerge -o %s -g %s -s %s -p 10 %s' % (cuffdir+'/merged_asm', cfg.gtf, cfg.genome, cuffdir+'/assemblies.txt')
	NGSTools.writeCommands(command, cuffdir+'/cuffmerge.sh', _run)

	# cuffdiff #
	cuffCondition = []
	cuffBam = {}
	for c in condition.keys():
		cuffCondition.append(c)
		cuffBam[c] = ','.join(condition[c].values())

	if len(cuffCondition) != 2:
		sys.exit('Error: condition')

	command = 'cuffdiff -o %s -b %s -p 10 -L %s -u %s %s %s' % (cuffdir+'/cuffdiff', cfg.genome, cuffCondition[0]+','+cuffCondition[1], cuffdir+'/merged_asm/merged.gtf', cuffBam[cuffCondition[0]], cuffBam[cuffCondition[1]])
	NGSTools.writeCommands(command, cuffdir+'/cuffdiff.sh', _run)

if DESeq2:
	# deseq2 #
	deseqDir = os.path.join(args.outDir, 'DESeq')
	try:
		os.mkdir(deseqDir)
	except:
		pass
	
	Rcommand = NGSTools.deseq2(countsFiles, deseqDir)
	with open(deseqDir+'/deseq2.R', 'w') as shell:
		shell.write(Rcommand)
	if _run:
		os.system('R %s' % Rcommand)


########################  DEU calling ########################
if DEXSeq:
	# DEXSeq #
	dexseqDir = os.path.join(args.outDir, 'DEXSeq')
	try:
		os.mkdir(dexseqDir)
	except:
		pass


########################  SNP calling ########################
if GATK:
	# GATK HC (in RNA-seq mode) call variations
	outdir = os.path.join(args.outDir, 'variation')
	NGSTools._mkdir(outdir)

	for condition in set(finalBam.values()):
		bams = []

		for bam in finalBam.keys():
			if finalBam[bam] == condition:
				bams.append(bam)

		out = os.path.join(outdir, condition)
		NGSTools._mkdir(out)

		mergedBam = '%s/%s.bam' % (out, condition)
		script = NGSTools.picard_merge(bams, mergedBam, cfg)

		NGSTools.writeCommands(command, '%s/picard_merge_%s.sh' % (out, condition), _run)


		script = NGSTools.GATK_HC(mergedBam, cfg, out, condition)

		NGSTools.writeCommands(command, '%s/gatk_HC_%s.sh' % (out, condition), _run)	


		rawVcf = '%s/%s.raw.snps.indels.vcf' % (out, condition)
		fltVcf = '%s/%s.flt.snps.indels.vcf' % (out, condition)
		script = NGSTools.GATK_filter(mergedBam, rawVcf, fltVcf, cfg)

		NGSTools.writeCommands(command, '%s/gatk_filter_%s.sh' % (out, condition), _run) 














