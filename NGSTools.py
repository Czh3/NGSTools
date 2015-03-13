#-*- coding: UTF-8 -*-
import os
import re
from time import ctime
from ConfigParser import ConfigParser
from string import Template


#############################################
#Name:	Next Generation Sequencing Tools
#Author:  Czh3 <zhangchao3@hotmail.com>
#Date:	2015-2-9
#Version: v0.1
#############################################


def isFile(file):

	if os.path.isfile(os.path.expanduser(file)):
		return True
	else:
		return False


def safeOpen(file, mode='r'):

	if file.endswith('.gz'):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)


def _isRun(shell, logic):

	if logic:
		shellLog = shell.replace('.sh', '.log')
		print 'SUBMIT [%s]: %s' % (ctime(), shell)
		os.system('sh %s 2> %s' % (shell, shellLog))
		print 'DONE [%s]: %s' % (ctime(), shell)
	else:pass


def _mkdir(dir):

	if not os.path.exists(dir):
		os.system('mkdir -p '+dir)


def writeCommands(command, file, run):

	'''write commands to the file'''
	with safeOpen(file, 'w') as shell:
		shell.write(command)
	_isRun(file, run)


def phred64to33(inputFq, outputFq):
	'''convert phred64 scaling fastq file format to phred33.'''

	from Bio import SeqIO
	
	if inputFq.endswith('.gz'):
		tmpIN = re.sub(r'.gz$', '', inputFq)
		assert not os.system('gzip -cdf %s > %s' % (inputFq, tmpIN))
	else:
		tmpIN = inputFq

	if outputFq.endswith('.gz'):
		tmpOUT = re.sub(r'.gz$', '', outputFq)
		SeqIO.covert(tmpIN, "fastq-illumina", tmpOUT, "fastq-sanger")
		assert not os.system('gzip -f %s' % tmpOUT)
	else:
		SeqIO.covert(tmpIN, "fastq-illumina", outputFq, "fastq-sanger")




def picard_merge(bamsList, outBam, cfg):
	'''merge bam with picard'''


	command = '\\\n\t'.join(['java -Xmx5g -jar %s/MergeSamFiles.jar' % cfg.picard,
								 'INPUT='+' INPUT='.join(bamsList),
								 'OUTPUT=' + outBam,
								 'TMP_DIR=' + os.path.dirname(outBam),
								 'USE_THREADING=true',
								 'SORT_ORDER=coordinate',
								 'VALIDATION_STRINGENCY=SILENT',
								 'MAX_RECORDS_IN_RAM=5000000'])
	command += '\nsamtools index ' + outBam

	return command

	

def GATK_HC(bam, cfg, outdir='', sampleName=''):
	''' GATK HaplotypeCaller '''

	command = '\\\n\t'.join(['java -Xmx5g -jar %s' % cfg.GATK,
									'-T HaplotypeCaller',
									'-R %s' % cfg.genome,
									'-I %s' % bam,
									'-dontUseSoftClippedBases',
									'-nct 10',
									#'--dbsnp %s',
									'-stand_call_conf 20',
									'-stand_emit_conf 20',
									'-o %s/%s.raw.snps.indels.vcf' % (outdir, sampleName)])

	return command


def GATK_filter(bam, inputVcf, outputVcf, cfg):
	# gatk filter

	tmp = re.sub(r'vcf$', 'ann.vcf', inputVcf)
	command = '\\\n\t'.join(['java -Xmx2g -jar %s' % cfg.GATK,
									'-T VariantAnnotator',
									'-R %s' % cfg.genome,
									'-I %s' % bam,
									'-o %s' % tmp,
									'-A coverage',
									'--variant %s' % inputVcf,
									'-L %s' % inputVcf,
									'-nt 10',
									'--dbsnp %s' % cfg.dbsnp]) 
	command += '\n'+'\\\n\t'.join(['java -Xmx2g -jar %s' % cfg.GATK,
									'-T VariantFiltration',
									'-R %s' % cfg.genome,
									'-o %s' % outputVcf,
									'--variant %s' % tmp,
									'--filterExpression "AB < 0.2 || MQ0 > 50"',
									'--filterName "Nov09filters"'])

	return command



class getConfig:
	'''get config file'''

	def __init__(self, cfgfile='~/.NGSTools.cfg'):

		#read config file
		config = ConfigParser()
		config.read(cfgfile)
		
		getConfig.genome = config.get('genome', 'fasta')
		getConfig.bowtie2Index = config.get('genome', 'bowtie2Index')
		getConfig.gtf = config.get('genome', 'gtf')

		getConfig.picard = config.get('tools', 'picard')
		getConfig.python = config.get('tools', 'python')
		getConfig.htseq = config.get('tools', 'htseq-count')
		getConfig.samtools = config.get('tools', 'samtools')
		getConfig.GATK = config.get('tools', 'GATK')
		getConfig.fastx = config.get('tools', 'fastx')

		getConfig.dbsnp = config.get('resource', 'dbsnp')
		getConfig.know_indel = config.get('resource', 'know_indel')


class NGSTools(getConfig):
	
	def __init__(self, sampleName, outdir, fq1, fq2='', quanlityBase='32', cfgfile='~/.NGSTools.cfg'):

		self.sampleName = sampleName
		self.outdir = os.path.abspath(outdir)
		_mkdir(self.outdir+'/data/'+self.sampleName)

		# link the raw fastq file to the work dir, and rename it
		assert isFile(fq1)

		if fq1.endswith('.gz'):
			self.fq1 = self.outdir+'/data/'+sampleName+'/'+sampleName+'.1.fq.gz'
		else:
			self.fq1 = self.outdir+'/data/'+sampleName+'/'+sampleName+'.1.fq'

		assert not os.system('ln -sf %s %s' % (fq1, self.fq1))

		if fq2 != '':
			assert isFile(fq2)
			if fq2.endswith('.gz'):
				self.fq2 = self.outdir+'/data/'+sampleName+'/'+sampleName+'.2.fq.gz'
			else:
				self.fq2 = self.outdir+'/data/'+sampleName+'/'+sampleName+'.2.fq'
			assert not os.system('ln -sf %s %s' % (fq2, self.fq2))
		else:
			self.fq2 = ''

		self.quanlityBase = quanlityBase
		
		self.bam = sampleName+'.bam'
		
		#read config file
		getConfig.__init__(self, cfgfile)
		
		
	

	def cutadapter(self, adapter5='AGATCGGAAGAGCGTCGTGTAGGGAAA', adapter3='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', run=True):
		'''cut illumina sequencing adapter'''

		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)
		
		if self.fq1.endswith('.gz'):
			cleanFq1 = re.sub(r'fq.gz$', 'rmAD.fq', os.path.abspath(self.fq1))
		else:
			cleanFq1 = re.sub(r'fq$', 'rmAD.fq', self.fq1)

		if self.fq2 != '':
			if self.fq2.endswith('.gz'):
				cleanFq2 = re.sub(r'fq.gz$', 'rmAD.fq', os.path.abspath(self.fq2))
			else:
				cleanFq2 = re.sub(r'fq$', 'rmAD.fq', self.fq2)

		if self.fq2 != '':
			command1 = 'cutadapt -a %s -e 0.01 -m 30 -O 5 -q 15 -o %s %s' % (adapter5, cleanFq2, self.fq2)
			self.fq2 = cleanFq2
		else:
			command1 = ''

		command2 = 'cutadapt -a %s -e 0.01 -m 30 -O 5 -q 10 -o %s %s' % (adapter3, cleanFq1, self.fq1)
		self.fq1 = cleanFq1
		
		writeCommands(command1+'\n'+command2, myOutdir+'/cutadapter_'+self.sampleName+'.sh', run)

	def fastx_lowQual(self, q=5, p=50, run=True):
		''' remove low quality reads by fastx '''

		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)
	
		cleanFq1 = re.sub(r'fq$', 'highQ.fq.gz', self.fq1)
		command = '%s/fastq_quality_filter -q %s -p %s -z -i %s -o %s\n' % (self.fastx, q, p, self.fq1, cleanFq1)
		
		cleanFq2 = re.sub(r'fq$', 'highQ.fq.gz', self.fq2)
		command += '%s/fastq_quality_filter -q %s -p %s -z -i %s -o %s' % (self.fastx, q, p, self.fq2, cleanFq2)

		writeCommands(command, myOutdir+'/rm_lowQ_'+self.sampleName+'.sh', run)
	
		self.fq1 = cleanFq1
		self.fq2 = cleanFq2


	def fastqc(self, run=True):
		'''quality control'''
		
		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)

		command = 'fastqc -o %s -t 6 -d %s %s %s' % (myOutdir, myOutdir, self.fq1, self.fq2)
		writeCommands(command, myOutdir+'/fastqc_'+self.sampleName+'.sh', run)


	def tophat2(self, run=True):
		'''mapping with tophat2'''

		myOutdir = self.outdir+'/mapping/'+self.sampleName
		_mkdir(myOutdir)

		_phredQual = ''
		if self.quanlityBase == '64':
			_phredQual = '--phred64-quals'

		command = 'tophat -p 6 -G %s %s -o %s %s %s %s' % (self.gtf, _phredQual, myOutdir, self.bowtie2Index, self.fq1, self.fq2)
		command += '\nmv %s/accepted_hits.bam %s' % (myOutdir, myOutdir+'/'+self.bam)
		
		writeCommands(command, myOutdir+'/tophat_'+self.sampleName+'.sh', run)
		self.bam = myOutdir+'/'+self.bam
		return self.bam


	def bowtie2(self, mode='--local', run=True):
		'''mapping with bowtie2 '''

		myOutdir = self.outdir+'/mapping/'+self.sampleName
		_mkdir(myOutdir)

		_phredQual = ''
		if self.quanlityBase == '64':
			_phredQual = '--phred64'

		command = 'bowtie2 %s -x %s -p 6 %s -1 %s ' % (_phredQual, self.bowtie2Index, mode, self.fq1)

		if self.fq2 != '':
			command += '-2 %s ' % self.fq2

		command += '| %s view -bS -h - > %s' % (self.samtools, myOutdir+'/'+self.bam)

		writeCommands(command, myOutdir+'/bowtie2_'+self.sampleName+'.sh', run)

		self.bam = myOutdir+'/'+self.bam
		return self.bam
		

	def bwa(self, run=True):
		'''mapping with bwa mem'''

		myOutdir = self.outdir+'/mapping'
		_mkdir(myOutdir)

		if self.quanlityBase == '64':
			print "The BWA-MEM does not support the phred64 scaling. Need to convert."
			if run:
				print "\tConverting to phred33."
				phred64to33(self.fq1, self.fq1.replace('.1.fq', '.phred33.1.fq'))
				self.fq1 = self.fq1.replace('.1.fq', '.phred33.1.fq')
				phred64to33(self.fq2, self.fq2.replace('.2.fq', '.phred33.2.fq'))
				self.fq2 = self.fq2.replace('.2.fq', '.phred33.2.fq')
				print "\tConvert done."

		command = 'bwa mem -t 6 -M -R "@RG\\tID:%s\\tSM:%s" %s %s ' % (self.sampleName, self.sampleName, self.genome, self.fq1)

		if self.fq2 != '':
			command += self.fq2

		command += ' | %s view -bS - > %s' % (self.samtools, myOutdir+'/'+self.bam)

		writeCommands(command, myOutdir+'/bwa_mem_'+self.sampleName+'.sh', run)
	
		self.bam = myOutdir+'/'+self.bam
		return self.bam


	def samtools_sort(self, run=True):
		'''sort bam with samtools'''

		sortedBam = re.sub(r'.bam$', '.sort.bam', self.bam)
		command = '%s sort -@3 -m 2G %s %s' % (self.samtools, self.bam, sortedBam)
		
		self.bam = sortedBam
		writeCommands(command, self.outdir+'/mapping/samtools_sort_'+self.sampleName+'.sh', run)

		return self.bam


	def picard_merge(self, bamsList, mergedSamplesName, run=True):
		'''merge bam with picard'''

		mergeBam = self.outdir+'/mapping/%s.bam' % mergedSamplesName

		command = '\\\n\t'.join(['java -Xmx5g -jar %s/MergeSamFiles.jar' % self.picard,
								 'INPUT='+' INPUT='.join(bamsList),
								 'OUTPUT='+mergeBam,
								 'TMP_DIR='+self.outdir,
								 'USE_THREADING=true',
								 'SORT_ORDER=coordinate',
								 'VALIDATION_STRINGENCY=SILENT',
								 'MAX_RECORDS_IN_RAM=5000000'])
		command += '\nsamtools index '+mergeBam

		self.bam = mergeBam
		writeCommands(command, self.outdir+'/maping/picard_mergebam_'+self.sampleName+'.sh', run)

		return self.bam


	def picard_rmdup(self, remove=True, run=True):
		'''mark or remove PCR duplicates in bam with picard'''

		rmdupBam = re.sub(r'.bam$', '.rmdup.bam', self.bam)

		if remove:
			remove = 'true'
		else:
			remove = 'false'

		METRICS_FILE = re.sub(r'.bam$', '.metrics', self.bam)

		command = '\\\n\t'.join(['java -Xmx6g -jar %s/MarkDuplicates.jar' % self.picard,
								 'TMP_DIR='+self.outdir,
								 'INPUT='+self.bam,
								 'OUTPUT='+rmdupBam,
								 'METRICS_FILE='+self.outdir+'/mapping/'+METRICS_FILE,
								 'VALIDATION_STRINGENCY=SILENT',
								 'ASSUME_SORTED=true',
								 'REMOVE_DUPLICATES='+remove,
								 'MAX_RECORDS_IN_RAM=5000000'])
		command += '\nsamtools index '+rmdupBam

		self.bam = rmdupBam
		writeCommands(command, self.outdir+'/mapping/'+self.sampleName+'/picard_rmdup_'+self.sampleName+'.sh', run)

		return self.bam


	def samtools_idxstats(self, run=True):
		'''stats mapping rate with samtools idxstats'''

		command = '%s idxstats %s > %s' % (self.samtools, self.bam, re.sub(r'.bam', '.stats', self.bam))

		writeCommands(command, self.outdir+'/mapping/samtools_idxstats_'+self.sampleName+'.sh', run)


	def HTSeq_count(self, run=True):
		'''count reads with HTSeq '''

		countFile = re.sub(r'.bam$', '.counts', self.bam)
		command = '%s %s -r pos -s no -f bam -a 10 %s > %s' % (self.python, self.htseq, self.bam, countFile)

		writeCommands(command, self.outdir+'/mapping/'+self.sampleName+'/htseq_count_'+self.sampleName+'.sh', run)

		return countFile


	def splitN(self, bam='', run=True):
		'''GATK SplitNCigarReads: Splits reads that contain Ns in their cigar string'''

		if bam == '':
			bam = self.bam

		outbam = re.sub(r'bam$', 'splitN.bam', bam)
		outdir = os.path.join(self.outdir, 'variation')
		_mkdir(outdir)

		command = '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T SplitNCigarReads',
									'-I %s' % self.bam,
									'--out %s' % outbam])
		self.bam = outbam

		writeCommands(command, outdir+'/gatk_splitN_'+self.sampleName+'.sh', run)

		return self.bam


	def realn(self, bam='', run=True):
		''' GATK realn '''

		if bam == '':
			bam = self.bam

		outbam = re.sub(r'bam$', 'realn.bam', bam)
		outdir = os.path.join(self.outdir, 'variation')
		_mkdir(outdir)

		# RealignerTargetCreator
		command = '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T RealignerTargetCreator',
									'-R %s' % self.genome,
									'-I %s' % self.bam,
									'-nt 10',
									#'-know %s'
									'-o %s/forIndelRealigner.intervals' % outdir ])
		
		# IndelRealigner
		command += '\n' + '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T IndelRealigner',
									'-R %s' % self.genome,
									'-I %s' % self.bam,
									'-targetIntervals %s' % outdir+'/forIndelRealigner.intervals',
									#'-know %s',
									'-o %s' % outbam ])

		self.bam = outbam
		writeCommands(command, outdir+'/gatk_realn_'+self.sampleName+'.sh', run)

		return self.bam


	def recal(self, bam='', run=True):
		''' GATK recal '''

		if bam == '':
			bam = self.bam

		outbam = re.sub(r'bam$', 'recal.bam', bam)
		outdir = os.path.join(self.outdir, 'variation')
		_mkdir(outdir)

		# BaseRecalibrator
		command = '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T BaseRecalibrator',
									'-I %s' % bam,
									'-R %s' % self.genome,
									'-nct 10',
									#'-knowSite %s',
									#'-knowSite %s',
									'-o %s/recal_data.table' % outdir ])

		# PrintReads
		command += '\n' + '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T PrintReads',
									'-R %s' % self.genome,
									'-I %s' % bam,
									'-nct 10',
									'-BQSR %s/recal_data.table' % outdir,
									'-o %s' % outbam ])

		self.bam = outbam
		writeCommands(command, outdir+'/gatk_recal_'+self.sampleName+'.sh', run)

		return self.bam


	def call(self,  run=True):
		''' GATK HaplotypeCaller '''

		outdir = os.path.join(self.outdir, 'variation', self.sampleName)
		_mkdir(outdir)

		command = GATK_HC(self.bam, self, outdir, self.sampleName)

		writeCommands(command, outdir+'/gatk_HC_'+self.sampleName+'.sh', run)

		self.rawVcf = '%s/%s.raw.snps.indels.vcf' % (outdir, self.sampleName)

		return self.rawVcf

	def filter(self, ):

		outdir = os.path.join(self.outdir, 'variation', self.sampleName)
		_mkdir(outdir)

		self.fltVcf = '%s/%s.flt.snps.indels.vcf' % (outdir, self.sampleName)
		GATK_filter(self.bam, self.rawVcf, self.fltVcf, self)

		return self.fltVcf



def deseq2(sampleCountPath_Condition, outdir):
	'''
	sampleCountPath_Condition is a hash:
	{sample1_count_path : condition,
	sample2_count_path : condition}
	return: the Rscript doing DEG calling using DESeq2
	'''

	Rscript = '''
library("DESeq2")

sampleFiles <- c($sampleF)

sampleCondition <- c($sampleC)
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("$C1", "$C2"))

dds <- DESeq(ddsHTSeq)
res <- results(dds)
res <- res[order(res$padj), ]

resSig <- subset(res, padj < 0.1)
write.table(as.data.frame(resSig), sep="\\t", file=paste($outdir, "DESeq.out.xls", collapse="/"))

plotMA(dds, ylim=c(-2,2), main="deseq2")
dev.copy(png, paste($outdir, "deseq2_MAplot.png", collapse="/")
dev.off()

'''

	sampleF = '", "'.join(sorted(sampleCountPath_Condition.keys()))
	sampleF = '"' + sampleF + '"'

	sampleC = '"'
	for i in sorted(sampleCountPath_Condition.keys()):
		sampleC += sampleCountPath_Condition[i] + '","'
	sampleC = sampleC[:-2] 

	try:
		(C1, C2) = set(sampleCountPath_Condition.values())
	except:
		print 'ERROR: condition must be 2.'

	Rscript = Template(Rscript)
	Rscript = Rscript.safe_substitute(sampleF=sampleF, sampleC=sampleC, C1=C1, C2=C2, outdir=outdir)

	return Rscript


def dexseq():
	'''under construction '''
	pass























