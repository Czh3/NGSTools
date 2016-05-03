#!/usr/bin/env python
# -*- encoding=utf8 -*-

import sys
#sys.path.append('/home/zhangc/bin/git/NGSTools')

import NGSTools
import argparse
import re
import os

from multiprocessing import Process, Manager

#argparse arguments
parser = argparse.ArgumentParser(description='A pipeline for WGS data analysis. <zhangchao3@hotmail.com>',
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-s', '--sampleList',
                    help="sample list for RNA samples information.\n"
                    " A file each line contains:\n"
                    "sampleID\tsampleName\tfastq1Path\tfastq2Path",
                    required=True)

parser.add_argument('-o', '--outDir',
                    help='The pipeline output dir',
                    default='out')

parser.add_argument('-c', '--config',
                    help='the config file of NGSTools package.',
                    default='~/.NGSTools.cfg')

parser.add_argument('-a', '--analysis',
                    help='analysis of the pipeline to do.\n'
                    'Here is some software to choose to analy\n'
                    '[1:QC, quality control\n'
                    ' 2:Mapping, align the reads to reference genome\n'
					' 3:picard_rmdup, remove PCR duplicates using picard\n'
                    ' 4:mpileup, call SNP using samtools mpileup]\n',
					default='1,2,3')

parser.add_argument('-b', '--qbase',
					help="quality base of base calling: 33(default) or 64\n",
					default='33')

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

QC = BWA  = RMDUP = MPILEUP = False

if 1 in analy:
	QC = True
if 2 in analy:
	BWA = True
if 3 in analy:
	RMDUP = True
if 4 in analy:
	MPILEUP = True


if args.debug:
	_run = False
else:
	_run = True



#read config file
cfg = NGSTools.getConfig(os.path.abspath(args.config))

# sample list parse
def sampleListParser():
	''' parse the sample list file.
	return a dictionary:
	{
		sampleName1:
		{
			sample ID:
			[
				fastq1 file path,
				fastq2 file path
			]
		},
		sampleName2:
		{
			...
		}
	}
	'''
	sampleListDict = {}

	for line in open(args.sampleList):
		if line.startswith('#') or line == "\n":
			continue
		# "sampleID\tsampleName\tfastq1Path\tfastq2Path"
		cols = line.strip().split()

		# fastq exists
		if not (os.path.exists(cols[2]) and os.path.exists(cols[3])):
			os.exit("fastq file don't exists!")

		if sampleListDict.has_key(cols[1]):
			if sampleListDict[cols[1]].has_key(cols[0]):
				sys.exit('duplication sample IDs are not allowed')
			else:
				sampleListDict[cols[1]][cols[0]] = [cols[2], cols[3]]
		else:
			sampleListDict[cols[1]] = {cols[0] : [cols[2], cols[3]]}

	return sampleListDict





def processSampleFromLibarary(sampleID, sampleIDList, libraryBamFileList):
	''''''

	########################## 0. init #########################
	#__init__(self, sampleName, outdir, fq1, fq2='', qualityBase='33', cfgfile='~/.NGSTools.cfg')

	mySample = NGSTools.NGSTools(sampleID, args.outDir, fq1=sampleIDList[0], fq2=sampleIDList[1], qualityBase=args.qbase, cfgfile=os.path.abspath(args.config))

	########################## 1. QC  #########################

	if QC:
		mySample.cutadapter(run=_run)
		mySample.QC_fastqc(run=_run)

	########################## 2. Mapping  #######################
	if BWA:
		finalBam = mySample.bwa(run=_run)
		#finalBam = mySample.bowtie2(mode='--end-to-end', run=_run)
		finalBam = mySample.samtools_sort(run=_run)

	libraryBamFileList.append(finalBam)

def processSample(sampleName, sampleNameDict):
	'''pipeline for sample '''

	# process communication
	mng = Manager()
	libraryBamFileList = mng.list()
	# This libraryBamFileList is a list contains seval bams from one sample.

	# init mult-processer
	record = []

	for sampleID in sampleNameDict:

		Processer = Process(name = sampleID, target = processSampleFromLibarary, args = (sampleID, sampleNameDict[sampleID], libraryBamFileList, ))

		Processer.start()

		record.append(Processer)

	# wait for processer
	for proc in record:
		proc.join()
		

	######################	2.1 post mapping  ####################
	bamFileDir = os.path.dirname(libraryBamFileList[0])
	mergedBamFilePath = os.path.join(bamFileDir, "%s_merged.bam" % sampleName)
	finalBamFilePath = mergedBamFilePath

	# merge the bam from each lane to one final bam
	command = NGSTools.picard_merge(libraryBamFileList, mergedBamFilePath, cfg)
	NGSTools.writeCommands(command, bamFileDir+'/picard_mergebam_'+sampleName+'.sh', run=_run)


	###################### 2.2. filter bam  ######################
	#command = 'samtools view -Sb -h -f 2 -q 10 %s > %s ' % (mergedBamFilePath, finalBamFilePath)
	#command = 'samtools view -Sb -h -q 10 %s > %s ' % (mergedBamFilePath, finalBamFilePath)
	#NGSTools.writeCommands(command, bamFileDir+'/filterBam_'+sampleName+'.sh', run=_run)

	###################### 3. romove duplicates ##################
	if RMDUP:
		command = NGSTools.picard_rmdup(mergedBamFilePath, True, cfg)
		NGSTools.writeCommands(command, bamFileDir+'/picard_rmdup_'+sampleName+'.sh', run=_run)

		finalBamFilePath = re.sub(r'.bam$', '.rmdup.bam', finalBamFilePath)


	####################### 4. call SNP ######################
	if MPILEUP:
		SNP_out = os.path.join(args.outDir, "SNP", sampleName)
		NGSTools._mkdir(SNP_out)

		command, rawVcf = NGSTools.bcftools_call(finalBamFilePath, cfg, outdir=SNP_out, sampleName=sampleName)
		NGSTools.writeCommands(command, SNP_out+'/bcftools_call_'+sampleName+'.sh', run=_run)
		
		outVcf = rawVcf.replace("vcf$", "flt.vcf")
		command = NGSTools.bcftools_filter(rawVcf, outVcf, cfg)
		NGSTools.writeCommands(command, SNP_out+'/bcftools_filter_'+sampleName+'.sh', run=_run)


########### main function

sampleListDict = sampleListParser()

# init mult-processer
records = []
for sampleName in sampleListDict:
	
	P = Process(name = sampleName, target = processSample, args = (sampleName, sampleListDict[sampleName], ))
	
	P.start()

	records.append(P)

for rec in records:
	rec.join()


	


