#-*- ceoding: UTF-8 -*-
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


	if len(bamsList) == 1:
		command = 'ln -sf %s %s' % (bamsList[0], outBam)

	else:
		command = '\\\n\t'.join(['java -Xmx5g -jar %s/MergeSamFiles.jar' % cfg.picard,
								 'INPUT='+' INPUT='.join(bamsList),
								 'OUTPUT=' + outBam,
								 'TMP_DIR=' + os.path.dirname(outBam),
								 'USE_THREADING=true',
								 'SORT_ORDER=coordinate',
								 'VALIDATION_STRINGENCY=SILENT'
								 ])
	command += '\nsamtools index ' + outBam

	return command

def picard_rmdup(INbam, remove, cfg):
	'''mark or remove PCR duplicates in bam with picard.
	output rmduped bam file named *.rmdup.bam by default'''

	rmdupBam = re.sub(r'.bam$', '.rmdup.bam', INbam)

	if remove:
		remove = 'true'
	else:
		remove = 'false'

	METRICS_FILE = re.sub(r'.bam$', '.metrics', rmdupBam)

	command = '\\\n\t'.join(['java -Xmx6g -jar %s/MarkDuplicates.jar' % cfg.picard,
								 'TMP_DIR='+os.path.dirname(rmdupBam),
								 'INPUT='+INbam,
								 'OUTPUT='+rmdupBam,
								 'METRICS_FILE='+METRICS_FILE,
								 'VALIDATION_STRINGENCY=SILENT',
								 'ASSUME_SORTED=true',
								 'REMOVE_DUPLICATES='+remove,
								 'MAX_RECORDS_IN_RAM=5000000'])
	command += '\n%s index %s' % (cfg.samtools, rmdupBam)


	return command


def GATK_HC(bam, cfg, outdir='', sampleName=''):
	''' GATK HaplotypeCaller .
	cfg: config file'''

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



def bcftools_call(bam, cfg, outdir='', sampleName=''):
	'''samtools variants calling '''

	command = '\\\n\t'.join(['%s  mpileup -q 1 -t DP,DV -ugf %s' % (cfg.samtools, cfg.genome),
								'%s' % bam,
								'| %s call -vmO z' % (cfg.bcftools),
								'-o %s/%s.vcf.gz' % (outdir, sampleName)])

	return (command, '%s/%s.vcf.gz' % (outdir, sampleName))


def bcftools_filter(inputVcf, outputVcf, cfg):
	''' snp/indel hard filter using samtools/bcftools '''

	command = '\\\n\t'.join(['%s filter -O z' % cfg.bcftools,
								'-o %s' % outputVcf,
								'-i \'%%QUAL>30 & DP > 10 & MQ > 40 \' %s' % (inputVcf)])

	command += '\ntabix -p vcf %s' % outputVcf

	return command


def methylation_extractor(bam, output_prefix, cfg):
	'''BS_seeker2 methylation extractor'''
	
	command = '\\\n\t'.join(['%s %s' % (cfg.python, cfg.bs_seeker),
							'--input=%s' % bam,
							'--output-prefix=%s' % output_prefix
							])

	return command


def bismark_methylation_extractor(bam, outdir, cfg):
	##bismark_methylation_extractor -p --cytosine_report --CX --no_overlap --multicore 5 --genome_folder ~/reference/human

	command = '\\\n\t'.join(['bismark_methylation_extractor -p --cytosine_report --CX --no_overlap',
							'--multicore 3',
							'-o %s' % outdir,
							'--genome_folder %s' % cfg.bismark_genome,
							'%s' % bam ])
	return command

class getConfig:
	'''get config file'''

	def __init__(self, cfgfile='~/.NGSTools.cfg'):

		#read config file
		config = ConfigParser()
		config.read(cfgfile)
		
		getConfig.genome = config.get('genome', 'fasta')
		getConfig.bismark_genome = config.get('genome', 'bismark_genome_dir')
		getConfig.genomeFai = config.get('genome', 'genomeFai')
		getConfig.bowtie2Index = config.get('genome', 'bowtie2Index')
		getConfig.gtf = config.get('genome', 'gtf')

		getConfig.picard = config.get('tools', 'picard')
		getConfig.python = config.get('tools', 'python')
		getConfig.htseq = config.get('tools', 'htseq-count')
		getConfig.samtools = config.get('tools', 'samtools')
		getConfig.bcftools = config.get('tools', 'bcftools')
		getConfig.GATK = config.get('tools', 'GATK')
		getConfig.fastx = config.get('tools', 'fastx')
		getConfig.cutadapt = config.get('tools', 'cutadapt')
		getConfig.fastqc = config.get('tools', 'fastqc')
		getConfig.bowtie = config.get('tools', 'bowtie')
		getConfig.tophat = config.get('tools', 'tophat')
		getConfig.gfold = config.get('tools', 'gfold')
		getConfig.trim_galore = config.get('tools', 'trim_galore')
		getConfig.bismark = config.get('tools', 'bismark')
		getConfig.bs_seeker = config.get('tools', 'bs_seeker')

		getConfig.dbsnp = config.get('resource', 'dbsnp')
		getConfig.know_indel = config.get('resource', 'know_indel')

		getConfig.NGSTools = config.get('NGSTools', 'NGSTools')

class NGSTools(getConfig):
	
	def __init__(self, sampleName, outdir, fq1, fq2='-', libType='fr-unstranded', qualityBase='33', cfgfile='~/.NGSTools.cfg'):

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

		# PE reads or SE reads
		if fq2 != '-':
			assert isFile(fq2)
			if fq2.endswith('.gz'):
				self.fq2 = self.outdir+'/data/'+sampleName+'/'+sampleName+'.2.fq.gz'
			else:
				self.fq2 = self.outdir+'/data/'+sampleName+'/'+sampleName+'.2.fq'
			assert not os.system('ln -sf %s %s' % (fq2, self.fq2))
		else:
			self.fq2 = ''

		self.libType = libType
		self.qualityBase = qualityBase
		
		self.bam = sampleName+'.bam'
		
		#read config file
		getConfig.__init__(self, cfgfile)
		
		
	

	def cutadapter(self, adapter5='AGATCGGAAGAGCGTCGTGTAGGGAAA', adapter3='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', run=True):
		'''cut illumina sequencing adapter'''

		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)
		
		if self.fq1.endswith('.gz'):
			cleanFq1 = re.sub(r'fq.gz$', 'rmAD.fq.gz', os.path.abspath(self.fq1))
		else:
			cleanFq1 = re.sub(r'fq$', 'rmAD.fq.gz', self.fq1)

		if self.fq2 != '':
			if self.fq2.endswith('.gz'):
				cleanFq2 = re.sub(r'fq.gz$', 'rmAD.fq.gz', os.path.abspath(self.fq2))
			else:
				cleanFq2 = re.sub(r'fq$', 'rmAD.fq.gz', self.fq2)

		if self.fq2 != '':
			# pair-end reads:
			tmpFq1 = re.sub(r'.fq.gz', '.tmp.fq', cleanFq1)
			tmpFq2 = re.sub(r'.fq.gz', '.tmp.fq', cleanFq2)
			command = '\\\n\t'.join(['%s -q 10 -a %s ' % (self.cutadapt, adapter3),
										'--minimum-length 30 -O 7 ',
										'--quality-base %s ' % self.qualityBase,
										'-o %s ' % tmpFq1,
										'-p %s ' % tmpFq2,
										'%s %s ' % (self.fq1, self.fq2)
									])
			command += '\n' + '\\\n\t'.join(['%s -q 15 -a %s ' % (self.cutadapt, adapter5),
										'--minimum-length 30 -O 7 ',
										'--quality-base %s ' % self.qualityBase,
										'-o %s ' % cleanFq2,
										'-p %s ' % cleanFq1,
										'%s %s ' % (tmpFq2, tmpFq1)
									])
			command += '\nrm %s %s' % (tmpFq1, tmpFq2)
			self.fq1 = cleanFq1
			self.fq2 = cleanFq2
		else:
			# Single-end reads
			command = '\\\n\t'.join(['%s -q 10 --minimum-length 30 -O 7 ' % self.cutadapt,
										'--quality-base %s ' % self.qualityBase,
										'-a %s ' % adapter3,
										'-o %s %s' % (cleanFq1, self.fq1)
									])
			self.fq1 = cleanFq1
		
		writeCommands(command, myOutdir+'/cutadapter_'+self.sampleName+'.sh', run)


	def trim_Galore(self, rrbs, adapter5='AGATCGGAAGAGCGTCGTGTAGGGAAA', adapter3='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', run=True):
		''' trim galore for RRBS sequencing '''

		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)


		if self.qualityBase == '33':
			qualityBase = '--phred33'
		elif self.qualityBase == '64':
			qualityBase = '--phred64'
		else:
			qualityBase = ''
			print 'WARNING: wrong quality scores for you input: ' + self.qualityBase

		if rrbs:
			_RRBS = '--rrbs'
		else:
			_RRBS = ''

		command = '\\\n\t'.join(['%s %s' % (self.trim_galore, qualityBase),
								'--gzip --paired %s' % _RRBS,
								'-a %s' % adapter3,
								'-a2 %s' % adapter5,
								'--output_dir %s' % myOutdir,
								'%s %s' % (self.fq1, self.fq2)
								])

		writeCommands(command, myOutdir+'/trim_galore_'+self.sampleName+'.sh', run)

		self.fq1 = os.path.join( myOutdir, re.sub(r".fq.gz$", '_val_1.fq.gz', os.path.basename(self.fq1)))
		self.fq2 = os.path.join( myOutdir, re.sub(r".fq.gz$", '_val_2.fq.gz', os.path.basename(self.fq2)))


	def fastx_lowQual(self, q=5, p=50, run=True):
		''' remove low quality reads by fastx '''
		
		#!!! this function now is discard !!!
		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)
		
		cleanFq1 = re.sub(r'fq$', 'highQ.fq.gz', self.fq1)
		command = '%s/fastq_quality_filter -q %s -p %s -z -i %s -o %s\nrm %s\n' % (self.fastx, q, p, self.fq1, cleanFq1, self.fq1)

		# NOTICE: the fastq_quality_filter input fq file must be unzipped		

		cleanFq2 = re.sub(r'fq$', 'highQ.fq.gz', self.fq2)
		command += '%s/fastq_quality_filter -q %s -p %s -z -i %s -o %s\nrm %s' % (self.fastx, q, p, self.fq2, cleanFq2, self.fq2)
		
		writeCommands(command, myOutdir+'/rm_lowQ_'+self.sampleName+'.sh', run)
		
		self.fq1 = cleanFq1
		self.fq2 = cleanFq2


	def rm_lowQual(self, q=5, p=50, a=0, run=True):
		''' remove low quality reads by qualityControl.py '''

		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)

		if self.fq2 != '':
			# Pair-end reads

			cleanFq1 = re.sub(r'fq.gz$', 'highQ.fq.gz', self.fq1)
			cleanFq2 = re.sub(r'fq.gz$', 'highQ.fq.gz', self.fq2)

			command = '\\\n\t'.join(['%s %s/qualityControl.py -1 %s ' % (self.python, self.NGSTools, self.fq1),
									'-2 %s ' % self.fq2,
									'-q %s -p %s -a %s ' % (q, p, a),
									'-o1 %s ' % cleanFq1,
									'-o2 %s ' % cleanFq2])
		
			#command += '\n' + 'rm %s %s' % (self.fq1, self.fq2)
			self.fq1 = cleanFq1
			self.fq2 = cleanFq2

		else:
			# Single-end reads:

			cleanFq1 = re.sub(r'fq.gz$', 'highQ.fq.gz', self.fq1)

			command = '\\\n\t'.join(['%s %s/qualityControl.py -1 %s ' % (self.python, self.NGSTools, self.fq1),
									'-q %s -p %s -a %s ' % (q, p, a),
									'-o1 %s ' % cleanFq1])

			command += '\n' + 'rm %s' % self.fq1
			self.fq1 = cleanFq1

		writeCommands(command, myOutdir+'/rm_lowQ_'+self.sampleName+'.sh', run)


	def QC_fastqc(self, run=True):
		'''quality control'''
		
		myOutdir = self.outdir+'/qc/'+self.sampleName
		_mkdir(myOutdir)

		command = '%s -o %s -t 6 -d %s %s %s' % (self.fastqc, myOutdir, myOutdir, self.fq1, self.fq2)
		writeCommands(command, myOutdir+'/fastqc_'+self.sampleName+'.sh', run)


	def tophat2(self, run=True):
		'''mapping with tophat2'''

		myOutdir = self.outdir+'/mapping/'+self.sampleName
		_mkdir(myOutdir)

		_phredQual = ''
		if self.qualityBase == '64':
			_phredQual = '--phred64-quals'

		command = '\\\n\t'.join(['%s -p 6' % self.tophat,
								'-G %s %s' % (self.gtf, _phredQual),
								'--library-type %s ' % self.libType,
								'--rg-id %s' % self.sampleName,
								'--rg-sample %s' % self.sampleName,
								'--rg-library %s' % self.sampleName,
								'-o %s %s %s %s' % (myOutdir, self.bowtie2Index, self.fq1, self.fq2)])

		command += '\nmv %s/accepted_hits.bam %s' % (myOutdir, myOutdir+'/'+self.bam)
		
		writeCommands(command, myOutdir+'/tophat_'+self.sampleName+'.sh', run)
		self.bam = myOutdir+'/'+self.bam
		return self.bam


	def bowtie2(self, mode='--local', run=True):
		'''mapping with bowtie2 '''

		myOutdir = self.outdir+'/mapping/'+self.sampleName
		_mkdir(myOutdir)

		_phredQual = ''
		if self.qualityBase == '64':
			_phredQual = '--phred64'

		command = '\\\n\t'.join(['%s %s' % (self.bowtie, _phredQual),
								'--rg-id %s' % self.sampleName,
								'-p 6 %s -x %s -1 %s' % (mode, self.bowtie2Index, self.fq1)])

		if self.fq2 != '':
			command += ' -2 %s ' % self.fq2

		command += '| %s view -bS -h -t %s - > %s' % (self.samtools, self.genomeFai, myOutdir+'/'+self.bam)

		writeCommands(command, myOutdir+'/bowtie2_'+self.sampleName+'.sh', run)

		self.bam = myOutdir+'/'+self.bam
		return self.bam
		

	def bwa(self, run=True):
		'''mapping with bwa mem'''

		myOutdir = self.outdir+'/mapping/'+self.sampleName
		_mkdir(myOutdir)

		if self.qualityBase == '64':
			print "The BWA-MEM does not support the phred64 scaling. Need to convert."
			if run:
				print "\tConverting to phred33."
				phred64to33(self.fq1, self.fq1.replace('.1.fq', '.phred33.1.fq'))
				self.fq1 = self.fq1.replace('.1.fq', '.phred33.1.fq')
				phred64to33(self.fq2, self.fq2.replace('.2.fq', '.phred33.2.fq'))
				self.fq2 = self.fq2.replace('.2.fq', '.phred33.2.fq')
				print "\tConvert done."

		command = 'bwa mem -t 6 -M -R "@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA" %s %s ' % (self.sampleName, self.sampleName, self.sampleName, self.genome, self.fq1)

		if self.fq2 != '':
			command += self.fq2

		command += ' | %s view -bS -t %s - > %s' % (self.samtools, self.genomeFai, myOutdir+'/'+self.bam)

		writeCommands(command, myOutdir+'/bwa_mem_'+self.sampleName+'.sh', run)
	
		self.bam = myOutdir+'/'+self.bam
		return self.bam

	def Bismark(self, run=True):
		'''mapping WGBS or RRBS reads to referent genome using bismark '''

		myOutdir = self.outdir+'/mapping/'
		_mkdir(myOutdir)

		if self.qualityBase == '33':
			qualityBase = '--phred33-quals'
		else:
			qualityBase = '--phred64-quals'


		command = '\\\n\t'.join(['%s --bowtie2' % self.bismark,
								'%s -p 5' % qualityBase,
								'--non_directional',
								'-o %s' % myOutdir,
								'--temp_dir %s/tmp' % myOutdir,
								'%s' % self.bismark_genome,
								'-1 %s' % self.fq1,
								'-2 %s' % self.fq2
								])
		self.bam = myOutdir + os.path.basename(self.fq1) + '_bismark_bt2_pe.bam'

		writeCommands(command, myOutdir+'/bismark_'+self.sampleName+'.sh', run)

		return self.bam


	def Bs_seeker2(self, rrbs=False, run=True):
		'''mapping WGBS or RRBS reads to referent genome using bs_seeker2 '''

		myOutdir = self.outdir+'/mapping/'
		_mkdir(myOutdir)

		self.bam = myOutdir + self.sampleName +'.bam'

		if rrbs:
			rrbsFlag = '-r'
		else:
			rrbsFlag = ''

		command = '\\\n\t'.join(['%s %s/bs_seeker2-align.py' % (self.python, self.bs_seeker),
								'--aligner=bowtie2 -g hg19.fa %s' % rrbsFlag,
								'--temp_dir ./ --bt-p 15',
								'-1 %s' % self.fq1,
								'-2 %s' % self.fq2,
								'-o %s' % self.bam
								])	

		writeCommands(command, myOutdir+'/bs_seeker2_'+self.sampleName+'.sh', run)

		return self.bam


	def samtools_sort(self, run=True):
		'''sort bam with samtools'''

		myOutdir = self.outdir+'/mapping/'+self.sampleName
		_mkdir(myOutdir)

		sortedBam = re.sub(r'.bam$', '.sort.bam', self.bam)
		command = '%s sort -@3 -m 2G %s %s' % (self.samtools, self.bam, sortedBam.rsplit('.', 1)[0])
		
		self.bam = sortedBam
		writeCommands(command, myOutdir+'/samtools_sort_'+self.sampleName+'.sh', run)

		return self.bam


	def picard_mergebam(self, bamsList, mergedSamplesName, run=True):
		'''merge bam with picard'''

		mergeBam = self.outdir+'/mapping/%s.bam' % mergedSamplesName

		command = picard_merge(bamList, mergeBam, self)

		self.bam = mergeBam
		writeCommands(command, self.outdir+'/mapping/picard_mergebam_'+self.sampleName+'.sh', run)

		return self.bam


	def rmdup(self, remove=True, run=True):
		'''mark or remove PCR duplicates in bam with picard'''

		rmdupBam = re.sub(r'.bam$', '.rmdup.bam', self.bam)

		command = picard_rmdup(self.bam, remove, self)

		self.bam = rmdupBam
		writeCommands(command, self.outdir+'/mapping/'+self.sampleName+'/picard_rmdup_'+self.sampleName+'.sh', run)

		return self.bam


	def picard_reorder(self, run=True):
		'''picard reorder'''

		reorderBam = re.sub(r'.bam$', '.reorder.bam', self.bam)
		
		command = '\\\n\t'.join(['java -Xmx6g -jar %s/ReorderSam.jar' % self.picard,
								'I=%s' % self.bam,
								'O=%s' % reorderBam,
								'R=%s' % self.genome])

		command += '\nrm %s' % self.bam

		self.bam = reorderBam
		writeCommands(command, self.outdir+'/mapping/'+self.sampleName+'/picard_reorder_'+self.sampleName+'.sh', run)

		return self.bam


	def samtools_idxstats(self, run=True):
		'''stats mapping rate with samtools idxstats'''

		command = '%s idxstats %s > %s' % (self.samtools, self.bam, re.sub(r'.bam', '.stats', self.bam))

		writeCommands(command, self.outdir+'/mapping/samtools_idxstats_'+self.sampleName+'.sh', run)


	def HTSeq_count(self, bam='', run=True):
		'''count reads with HTSeq '''

		if bam == '':
			bam = self.bam

		countFile = re.sub(r'.bam$', '.counts', bam)
		command = '%s -r pos -s no -f bam -a 10 %s %s > %s' % (self.htseq, bam, self.gtf, countFile)

		writeCommands(command, self.outdir+'/mapping/'+self.sampleName+'/htseq_count_'+self.sampleName+'.sh', run)

		return countFile

	def gfoldCount(self, run=True):
		''' count reads with gfold'''

		gfoldDir = os.path.join(self.outdir, 'GFold')
		_mkdir(gfoldDir)

		countFile = os.path.join(gfoldDir, self.sampleName+'.gfoldCount')

		command = '\\\n\t'.join(['%s view %s | %s count -ann %s ' % (self.samtools, self.bam, self.gfold, self.gtf),
								'-tag stdin ',
								'-o %s ' %  countFile])

		writeCommands(command, gfoldDir+'/gfold_count_%s.sh' % self.sampleName, run)

	def splitN(self, bam='', run=True):
		'''GATK SplitNCigarReads: Splits reads that contain Ns in their cigar string'''

		if bam == '':
			bam = self.bam

		outbam = re.sub(r'bam$', 'splitN.bam', bam)
		outdir = os.path.join(self.outdir, 'variation')
		_mkdir(outdir)

		command = '%s index %s\n' % (self.samtools, bam)

		command += '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T SplitNCigarReads',
									'-R %s' % self.genome,
									'-I %s' % bam,
									'-U ALLOW_N_CIGAR_READS',
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
									'-o %s/%sIndelRealigner.intervals' % (outdir, self.sampleName) ])
		
		# IndelRealigner
		command += '\n' + '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T IndelRealigner',
									'-R %s' % self.genome,
									'-I %s' % self.bam,
									'-targetIntervals %s' % outdir+'/'+self.sampleName+'IndelRealigner.intervals',
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
									'-o %s/%srecal_data.table' % (outdir, self.sampleName) ])

		# PrintReads
		command += '\n' + '\\\n\t'.join(['java -Xmx5g -jar %s' % self.GATK,
									'-T PrintReads',
									'-R %s' % self.genome,
									'-I %s' % bam,
									'-nct 10',
									'-BQSR %s/%srecal_data.table' % (outdir, self.sampleName),
									'-o %s' % outbam ])

		self.bam = outbam
		writeCommands(command, outdir+'/gatk_recal_'+self.sampleName+'.sh', run)

		return self.bam


	def GATK_call(self,  run=True):
		''' GATK HaplotypeCaller '''

		outdir = os.path.join(self.outdir, 'variation', self.sampleName)
		_mkdir(outdir)

		command = GATK_HC(self.bam, self, outdir, self.sampleName)

		writeCommands(command, outdir+'/gatk_HC_'+self.sampleName+'.sh', run)

		self.rawVcf = '%s/%s.raw.snps.indels.vcf' % (outdir, self.sampleName)

		return self.rawVcf

	def GATK_filter(self, run=True):

		outdir = os.path.join(self.outdir, 'variation', self.sampleName)
		_mkdir(outdir)


		self.fltVcf = '%s/%s.flt.snps.indels.vcf' % (outdir, self.sampleName)
		
		command = GATK_filter(self.bam, self.rawVcf, self.fltVcf, self)

		writeCommands(command, outdir+'/gatk_filter_'+self.sampleName+'.sh', run)

		return self.fltVcf



	def samtools_call(self, run=True):
		'''call snp/indel using samtools '''

		outdir = os.path.join(self.outdir, 'variation', self.sampleName)
		_mkdir(outdir)

		command, rawVcf = bcftools_call(self.bam, self, outdir, self.sampleName)

		writeCommands(command, outdir+'/samtools_call_'+self.sampleName+'.sh', run)

		self.rawVcf = rawVcf

		return self.rawVcf


	def samtools_filter(self, run=True):
		'''use samtools_filter(bam, inputVcf, outputVcf, cfg):'''

		outdir = os.path.join(self.outdir, 'variation', self.sampleName)

		self.filterVcf = re.sub(r'vcf.gz$', 'filter.vcf.gz', self.rawVcf)
		command = bcftools_filter(self.rawVcf, self.filterVcf, self)

		writeCommands(command, outdir+'/samtools_filter_'+self.sampleName+'.sh', run)

		return self.filterVcf

def deseq2(sampleCountPath_Condition, outdir):
	'''
	CAREFUL!!!, did not test! use before test.

	sampleCountPath_Condition is a hash:
	{sample1_count_path : condition|sample1_name,
	sample2_count_path : condition|sample2_name}
	return: the Rscript doing DEG calling using DESeq2
	'''

	Rscript = '''
library("DESeq2")

sampleFiles <- c($sampleF)

sampleCondition <- c($sampleC)
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition)

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("$C1", "$C2"))

dds <- DESeq(ddsHTSeq)
res <- results(dds)
res <- res[order(res$padj), ]

resSig <- subset(res, padj < 0.1)
write.table(as.data.frame(resSig), sep="\\t", file=paste("$outdir", "DESeq.out.xls", sep="/"))

pdf(paste("$outdir", "deseq2_MAplot.png", sep="/"))
plotMA(dds, ylim=c(-2,2), main="deseq2")
dev.off()

'''

	sampleF = '", "'.join(sorted(sampleCountPath_Condition.keys()))
	sampleF = '"' + sampleF + '"'

	sampleC = '"'
	for i in sorted(sampleCountPath_Condition.keys()):
		sampleC += sampleCountPath_Condition[i].split('|')[0] + '","'
	sampleC = sampleC[:-2] 

	try:
		(C1, C2) = set(sampleCountPath_Condition.values())
	except:
		print 'ERROR: condition must be 2.'
		(C1, C2) = list(set(sampleCountPath_Condition.values()))[0:2]

	C1 = C1.split('|')[0]
	C2 = C2.split('|')[0]
	Rscript = Template(Rscript)
	Rscript = Rscript.safe_substitute(sampleF=sampleF, sampleC=sampleC, C1=C1, C2=C2, outdir=outdir)

	return Rscript

def htseq2table(sampleCountPath_Condition, outfile):
	'''merge all sample's gene counts from HTseq-count to one file '''
	
	res = {}
	headerLine = ""
	for path in sampleCountPath_Condition.keys():
		headerLine += "\t" + sampleCountPath_Condition[path].split('|')[1]
		for line in open(path):
			col = line.strip().split()
			if res.has_key(col[0]):
				res[col[0]] += '\t' + col[1]
			else:
				res[col[0]] = '\t' + col[1]

	outputHandle = open(outfile, "w")
	outputHandle.write(headerLine + '\n')
	for i in res:
		outputHandle.write( i + res[i] + '\n')


def deseq(sampleCountPath_Condition, outdir):
	'''deseq to call DEGs '''

	htseq2table(sampleCountPath_Condition, os.path.join(outdir, 'data.table'))

	Rscript = '''
library(DESeq)
setwd("$outdir")

countTable = read.table( "data.table", header=TRUE, row.names=1 )
condition = colnames(countTable)

cds = newCountDataSet( countTable, condition )

##normalization
cds = estimateSizeFactors( cds )
sizeFactors( cds )
normCounts = counts( cds, normalized=TRUE )
cds = estimateDispersions( cds )
#cds = estimateDispersions( cds, method="blind", sharingMode = "fit-only" )

#boxplot
df = as.data.frame(normCounts)
mat = log10(df + 1)
library(gplots)
library(reshape)
mat.melt = melt(mat)
qplot(factor(variable), value, data=mat.melt, geom="boxplot", fill=factor(variable), alpha=I(.5))

# hclust
d = dist(t(df))
h = hclust(d)
plot(h)

df = df[rowSums(df)>0,]
mat=as.matrix(df)
mat.z = t(scale(t(mat)))
hm = heatmap.2(mat.z,main='gene expression heatmap',col=rev(redblue(125)),
          keysize=1,density.info="none", srtCol=45, margins = c(10,5),
          dendrogram="col", #key.title='expression level',key.xlab = "log10(CPM+1)",
          trace="none",labRow=NA, Rowv = TRUE, Colv=FALSE,
          scale="none")

#Calling differential expression
call_DEG = function(condition1, condition2){
  res = nbinomTest( cds, condition1, condition2 )
  #DESeq::plotMA(res)
  resSig = res[ res$pval < 0.05, ]
  resSig = na.omit(resSig)
  resSig <- resSig[ order(resSig$pval), ]
  outFile = paste(condition1, condition2, "DEG.xls", sep="_")
  outFile.all = paste(condition1, condition2, "ALL.xls", sep="_")
  write.table(resSig, file=outFile, quote=FALSE, sep="\\t", row.names=FALSE)
  write.table(res, file=outFile.all, quote=FALSE, sep="\\t", row.names=FALSE)
}

call_DEG("$C1", "$C2")
'''

	C1 = C2 = ''

	for samplePath in sampleCountPath_Condition.keys():
		C1 = sampleCountPath_Condition[samplePath].split('|')[0]
	for samplePath in sampleCountPath_Condition.keys():
		if sampleCountPath_Condition[samplePath].split('|')[0] != C1:
			C2 = sampleCountPath_Condition[samplePath].split('|')[0]
			break

	Rscript = Template(Rscript)

	try:
		Rscript = Rscript.safe_substitute(outdir = outdir, C1 = C1, C2 = C2)
	except:
		Rscript = Rscript.safe_substitute(outdir = outdir, C1 = 'Condition1', C2 = 'Condition2')
	finally:
		pass

	return Rscript

def dexseq():
	'''under construction '''
	pass








