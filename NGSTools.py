#-*- coding: UTF-8 -*-
import os
import re
from ConfigParser import ConfigParser

#############################################
#Name:    Next Generation Sequencing Tools
#Author:  Czh3 <zhangchao3@hotmail.com>
#Data:    2015-2-9
#Version: v0.1
#############################################


def isFile(file):
    if os.file.exists(file):
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
        print 'SUBMIT: '+shell
        os.system('sh %s > %s' % (shell, shellLog))
        print 'DONE: '+shell
    else:pass

def _mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

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


class NGSTools:
    
    def __init__(self, sampleName, outdir, fq1, fq2='', quanlityBase='32', config='~/.NGSTools.cfg'):
        self.sampleName = sampleName
        self.outdir = os.path.abs(outdir)
        _mkdir(outdir+'/data')
        
        # link the raw fastq file to the work dir, and rename it
        assert isFile(fq1)
        if fq1.endswith('.gz'):
            self.fq1 = outdir+'/data/'+sampleName+'.1.fq.gz'
        else:
            self.fq1 = outdir+'/data/'+sampleName+'.1.fq'
        assert not os.system('ln -sf %s %s' % (fq1, self.fq1))
        if fq2 != '':
            assert isfile(fq2)
            if fq2.endswith('.gz'):
                self.fq2 = outdir+'/data/'+sampleName+'.2.fq.gz'
            else:
                self.fq2 = outdir+'/data/'+sampleName+'.2.fq'
            assert not os.system('ln -sf %s %s' % (fq2, outdir+'/data/'+sampleName+'.2.fq'))
        else:
            self.fq2 = ''

        self.quanlityBase = quanlityBase
        
        self.bam = sampleName+'.bam'
        
        
        #read config file
        config = ConfigParser()
        config.read(cfgfile)
        self.genome = config.get('genome', 'fasta')
        self.bowtie2Index = config.get('genome', 'bowtie2Index')
        self.gtf = config.get('genome', 'gtf')
        self.picard = config.get('tools', 'picard')
        self.python = config.get('tools', 'python')
        self.htseq = config.get('tools', 'htseq-count')
    

    def cutadapter(self, adapter5='AGATCGGAAGAGCGTCGTGTAGGGAAA', adapter3='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', run=True):
        '''cut illumina sequencing adapter'''
        myOutdir = self.outdir+'/qc'
        _mkdir(myOutdir)
        
        if self.fq1.endswith('.gz'):
            cleanFq1 = re.sub(r'fq.gz$', 'clean.fq.gz', os.path.abs(self.fq1))
        else:
            cleanFq1 = re.sub(r'fq$', 'clean.fq', self.fq1)

        if self.fq2 != '':
            if self.fq2.endswith('.gz'):
                cleanFq2 = re.sub(r'fq.gz$', 'clean.fq.gz', os.path.abs(self.fq2))
            else:
                cleanFq2 = re.sub(r'fq$', 'clean.fq', self.fq2)

        if self.fq2 != '':
            command1 = 'cutadapt -a %s -e 0.01 -m 30 -O 5 -q 15 -o %s %s' % (adapter5, cleanFq2, self.fq2)
            self.fq2 = cleanFq2
        else:
            command1 = ''
        command2 = 'cutadapt -a %s -e 0.01 -m 30 -O 5 -q 10 -o %s %s' % (adapter3, cleanFq1, self.fq1)
        self.fq1 = cleanFq1
        
        writeCommands(command1+'\n'+command2, myOutdir+'/cutadapter.sh', run)


    def fastqc(self, run=True):
        '''quality control'''
        
        myOutdir = self.outdir+'/qc'
        _mkdir(myOutdir)
        command = 'fastqc -o %s -t 6 -d %s %s %s' % (outdir, outdir, self.fq1, self.fq2)
        writeCommands(command, myOutdir+'/fastqc.sh', run)

    def tophat2(self, run=True):
        '''mapping with tophat2'''
        myOutdir = self.outdir+'/mapping'
        _mkdir(myOutdir)
        
        _phredQual = ''
        if self.quanlityBase == '64':
            _phredQual = '--phred64-quals'
        command = 'tophat -p 6 -G %s %s -o %s %s %s %s' % (self.gtf, _phredQual, myOutdir, self.bowtie2Index, self.fq1, self.fq2)
        command += '\nmv %s/accepted_hits.bam %s' % (myOutdir, myOutdir+'/'+self.bam)
        
        writeCommands(command, myOutdir+'/tophat.sh', run)
        return myOutdir+'/'+self.bam

    def bowtie2(self, mode='--local', run=True):
        '''mapping with bowtie2 '''
        myOutdir = self.outdir+'/mapping'
        _mkdir(myOutdir)
        
        _phredQual = ''
        if self.quanlityBase == '64':
            _phredQual = '--phred64'
        command = 'bowtie2 %s -x %s -p 6 %s -1 %s ' % (_phredQual, self.bowtie2Index, mode, self.fq1)
        if self.fq2 != '':
            command += '-2 %s ' % self.fq2
        command += '| samtools view -bS -h - > %s' % myOutdir+'/'+self.bam

        writeCommands(command, myOutdir+'/bowtie2.sh', run)
        return myOutdir+'/'+self.bam
            
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
        command += ' | samtools view -bS - > %s' % self.bam

        writeCommands(command, myOutdir+'/bwa_mem.sh', run)
        return myOutdir+'/'+self.bam

    def samtools_sort(self, run=True):
        '''sort bam with samtools'''
        sortedBam = re.sub(r'.bam$', '.sort.bam', self.bam)
        command = 'samtools sort -@3 -m 2G %s %s' % (self.bam, sortedBam)
        
        self.bam = sortedBam
        writeCommands(command, self.outdir+'/mapping/samtools_sort.sh', run)
        return self.bam

    def picard_merge(self, bamsList, mergedSamplesName, run=True):
        '''merge bam with picard'''
        mergeBam = self.outdir+'/mapping/%s.bam' % mergedSamplesName
        command = '\\\n\t'.join(['java -Xmx5g -jar %s/MergeSamFiles.jar' % self.picard,
                                 'INPUT='+' INPUT='.join(bamsList),
                                 'OUTPUT='+mergeBam,
                                 'TMP_DIR='+self.outdir,
                                 'USE_THREADING=true',
                                 'SORT_ORDER=coordinate'
                                 'VALIDATION_STRINGENCY=SILENT',
                                 'MAX_RECORDS_IN_RAM=5000000'])
        command += '\nsamtools index '+mergeBam

        self.bam = mergeBam
        writeCommands(command, self.outdir+'/maping/picard_mergebam.sh', run)
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
        writeCommands(command, self.outdir+'/mapping/picard_rmdup.sh', run)
        return self.bam

    def samtools_idxstats(self, run=True):
        '''stats mapping rate with samtools idxstats'''
        command = 'samtools idxstats %s > %s' % (self.bam, re.sub(r'.bam', '.stats', self.bam))
        writeCommands(command, self.outdir+'/mapping/samtools_idxstats.sh', run)

    def HTSeq_count(self, run=True):
        '''count reads with HTSeq '''
        countFile = re.sub(r'.bam$', '.counts', self.bam)
        command = '%s %s -s no -f bam -a 10 %s > %s' % (self.python, self.htseq, self.bam, countFile)
        writeCommands(command, self.outdir+'/mapping/htseq_count.sh', run)
        return countFile






