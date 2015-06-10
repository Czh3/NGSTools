#-*- encoding=utf8 -*-
import sys
sys.path.append('/home/zhangc/bin/git/TEST')

import NGSTools
import argparse
import os

from multiprocessing import Process, Manager

#argparse arguments
parser = argparse.ArgumentParser(description='A pipeline of RNA_seq data analysis. <zhangchao3@hotmail.com>',
								formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-s', '--sampleList',
					help="sample list for RNA samples information.\n"
					" A file each line contains:\n"
					"sampleName\tsampleCondition\tfastq1Path\tfastq2Path",
					required=True)
parser.add_argument('-o', '--outDir',
					help='The pipeline output dir',
					default='out')
parser.add_argument('-d', '--dataType',
					help='fastq data type:\n'
					'raw data or clean data.\n'
					' if (clean data): not run cutadapter',
					choices=['raw', 'clean'],
					default='clean')
parser.add_argument('-a', '--analysis',
					help='analysis of the pipeline to do.\n'
					'Here is some software to choose to analy\n'
					'[1:QC, quality control\n'
					' 2:Mapping, align the reads to reference genome\n'
					' 3:Cufflinks, assemble with cufflinkes\n'
					' 4:DESeq2, call DEGs(different expression genes) using DESeq2 package\n'
					' 5:DEXSeq, call DEUs(different exon usages) using DEXSeq package\n'
					' 6:GATK, call SNP on mRNA using GATK.\n'
					' 7:GFold, call DEGs without biological replicates using GFold]',
					default='1,2')
parser.add_argument('-c', '--config',
					help='the config file of NGSTools package.',
					default='~/.NGSTools.cfg')
parser.add_argument('--debug',
					help='debug mode',
					default=False)
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

QC = Mapping = Cufflinks = DESeq2 = DEXSeq = GATK = GFold  = False

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
if 7 in analy:
	GFold = True

global _run
_run = ''
if args.debug == False:
	_run = True
else:
	_run = False

##
cfg = NGSTools.getConfig(os.path.abspath(args.config))



def processSample(line, condition, transcripts, countsFiles, finalBam):

	cols = line.strip().split('\t')

	if len(cols) == 3:
		# single end library
		fq2 = '-'
	else:
		# paired end
		fq2 = cols[3]

	sample = {
		'name' : cols[0],
		'condition' :	cols[1],
		'fq1' : cols[2],
		'fq2' : fq2,
		'bam' : ''
	}


	########################## 0. init #########################
	#__init__(self, sampleName, outdir, fq1, fq2='', quanlityBase='32', cfgfile='~/.NGSTools.cfg'):
	mySample = NGSTools.NGSTools(sample['name'], args.outDir, sample['fq1'], sample['fq2'], cfgfile=os.path.abspath(args.config))
	

	if QC:
		#################### 1. Quality Control ####################

		###### 1.1 cut adapter ######
		if args.dataType == 'raw':
			#mySample.cutadapter(adapter5='', adapter3='AATGATACGGCGACCACCGAGATCT', run = _run)
			mySample.cutadapter(run = _run)
			#mySample.rm_lowQual(run = _run)
		else:
			pass

		##### 1.2 fastqc #####
		mySample.QC_fastqc(run = _run)
	
	
	if Mapping:
		######################## 2. Mapping ########################
		
		sample['bam'] = mySample.tophat2(run = _run)


		if condition.has_key(sample['condition']):
			#condition[sample['condition']][sample['name']] = sample['bam']
			condition[sample['condition']] += ","+sample['bam']
		else:
			#condition[sample['condition']] = {sample['name'] : sample['bam']}
			condition[sample['condition']] = sample['bam']

		if GFold:

			# GFold count
			mySample.gfoldCount(run = _run)
		
		if DESeq2:
			# DESeq2
			count = mySample.HTSeq_count(run = _run)
			countsFiles[count] = sample['condition']

		if GATK:

			# remove duplicates
			mySample.rmdup(run = _run)

			# picard reorder
			mySample.picard_reorder(run = _run)

			# splitN
			mySample.splitN(run = _run)

			# realign
			realnBam = mySample.realn(run = _run)

			# recal need known SNP site

			# recal
			#recalBam = mySample.recal(run = _run)

			finalBam[realnBam] = sample['condition']

			# samtools call SNP/InDel
			mySample.samtools_call(run = _run)
			mySample.samtools_filter(run = _run)


	########################  DEGs calling preparation ########################

	if Cufflinks:
		##### 3. cufflinks #####
		cuffdir = os.path.join(args.outDir, 'cufflinks')
		if not os.path.exists(cuffdir):
			os.mkdir(cuffdir)

		# cufflinks #
		command = 'cufflinks -p 4 -g %s -o %s %s' % (cfg.gtf, os.path.join(cuffdir, sample['condition']+'_'+sample['name']), sample['bam'])
		NGSTools.writeCommands(command, cuffdir+'/cufflinks_%s.sh' % sample['name'], _run)

		transcripts.append(os.path.join(cuffdir, sample['condition']+'_'+sample['name'], 'transcripts.gtf'))






# process communication
mng = Manager()
condition = mng.dict()
transcripts = mng.list()
countsFiles = mng.dict()
finalBam = mng.dict()

# multi-process
record = []

for line in open(args.sampleList):
	if line.startswith('#') or line == '\n':
		continue
	
	sampleName = line.split('\t')[0]

	P = Process(name=sampleName, target=processSample, args=(line, condition, transcripts, countsFiles, finalBam))
	P.start()

	record.append(P)

for P in record:
	P.join()



########################  DEGs calling ########################
if Cufflinks:
	# cuffmerge #
	cuffdir = os.path.join(args.outDir, 'cufflinks')
	with open(cuffdir+'/assemblies.txt', 'w') as writer:
		writer.write('\n'.join(transcripts))

	command = 'cuffmerge -o %s -g %s -s %s -p 10 %s' % (cuffdir+'/merged_asm', cfg.gtf, cfg.genome, cuffdir+'/assemblies.txt')
	NGSTools.writeCommands(command, cuffdir+'/cuffmerge.sh', _run)

	# cuffdiff #

	if len(condition) != 2:
		print 'WARNING: condition'

	command = 'cuffdiff -o %s -b %s -p 10 -L %s -u %s %s %s' % (cuffdir+'/cuffdiff', cfg.genome, condition.keys()[0]+','+condition.keys()[1], cuffdir+'/merged_asm/merged.gtf', condition.values()[0], condition.values()[1])
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

if GFold:
	# GFold #
	gfoldDir = os.path.join(args.outDir, 'GFold')
	
	command = '\\\n\t'.join(['%s diff -s1 $1 -s2 $2 ' % cfg.gfold,
							'-suf .gfoldCount ',
							'-o $1"vs"$2'])

	with open(gfoldDir+'/gfold_diff.sh', 'w') as shell:
		shell.write(command)

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

		NGSTools.writeCommands(script, '%s/picard_merge_%s.sh' % (out, condition), False)


		script = NGSTools.GATK_HC(mergedBam, cfg, out, condition)

		NGSTools.writeCommands(script, '%s/gatk_HC_%s.sh' % (out, condition), False)	


		rawVcf = '%s/%s.raw.snps.indels.vcf' % (out, condition)
		fltVcf = '%s/%s.flt.snps.indels.vcf' % (out, condition)
		script = NGSTools.GATK_filter(mergedBam, rawVcf, fltVcf, cfg)

		NGSTools.writeCommands(script, '%s/gatk_filter_%s.sh' % (out, condition), False) 













