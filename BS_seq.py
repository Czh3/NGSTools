#!/usr/bin/env python
# -*- encoding=utf8 -*-

import sys
#sys.path.append('/home/zhangc/bin/git/NGSTools')

import NGSTools
import argparse
import os

from multiprocessing import Process, Manager

#argparse arguments
parser = argparse.ArgumentParser(description='A pipeline of Bisulfite-Seq(WGBS and RRBS) data analysis. <zhangchao3@hotmail.com>',
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

parser.add_argument('-r', '--rrbs',
					help='for RRBS library',
					action='store_true',
					default=False)

parser.add_argument('-a', '--analysis',
                    help='analysis of the pipeline to do.\n'
                    'Here is some software to choose to analy\n'
                    '[1:QC, quality control\n'
                    ' 2:bismark, align the reads to reference genome\n'
                    ' 3:bs_seeker2, align the reads to reference genome\n'
					' 4:picard_rmdup, remove PCR duplicates using picard\n'
                    ' 5:methylation_extractor, call methylation levels using the tools in bismark\n'
                    ' 6:swDMR, call methylation levels using swDMR]\n',
                    default='1,2')

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

QC = Bismark = Bs_seeker2 = Methylation_extractor = swDMR = False

if 1 in analy:
	QC = True
if 2 in analy:
	Bismark = True
if 3 in analy:
	Bs_seeker2 = True
if 4 in analy:
	RMDUP = True
if 5 in analy:
	Methylation_extractor = True
if 6 in analy:
	swDMR = True

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
		if line.startswith('#'):
			continue
		# "sampleID\tsampleName\tfastq1Path\tfastq2Path"
		cols = line.strip().split("\t")

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
	#__init__(self, sampleName, outdir, fq1, fq2='', quanlityBase='33', cfgfile='~/.NGSTools.cfg')

	mySample = NGSTools.NGSTools(sampleID, args.outDir, fq1=sampleIDList[0], fq2=sampleIDList[1], cfgfile=os.path.abspath(args.config))

	########################## 1. QC  #########################

	if QC:
		mySample.trim_Galore(args.rrbs, run=_run)
		mySample.QC_fastqc(run=_run)

	########################## 2. Mapping  #######################
	if Bismark:
		finalBam = mySample.Bismark(run=_run)
	elif Bs_seeker2:
		finalBam = mySample.Bs_seeker2(run=_run)
	else:
		print 'WARNING: not choose a mapping tools'

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
	finalBamFilePath = os.path.join(bamFileDir, "%s_final.bam" % sampleName)

	# merge the bam from each line to one final bam
	command = NGSTools.picard_merge(libraryBamFileList, finalBamFilePath, cfg)
	NGSTools.writeCommands(command, args.outDir+'/mapping/picard_mergebam_'+sampleName+'.sh', run=_run)

	###################### 3. romove duplicates ##################
	if RMDUP:
		command = NGSTools.picard_rmdup(finalBamFilePath, True, cfg)
		NGSTools.writeCommands(command, args.outDir+'/mapping/picard_rmdup_'+sampleName+'.sh', run=_run)

		finalBamFilePath = re.sub(r'.bam$', '.rmdup.bam', finalBamFilePath)

	
	##################### 4. DMR calling ##################
	if Methylation_extractor:
		NGSTools.methylation_extractor(finalBamFilePath, bamFileDir+'/'+sampleName, cfg)


# main function

sampleListDict = sampleListParser()

# init mult-processer
records = []
for sampleName in sampleListDict:

	P = Process(name = sampleName, target = processSample, args = (sampleName, sampleListDict[sampleName], ))
	
	P.start()

	records.append(P)

for rec in records:
	rec.join()


	













































