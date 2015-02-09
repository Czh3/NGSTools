import os
import re
import ConfigParser

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
        return Flase

def safeOpen(file, mode='r'):
    if file.endswith(.gz):
        import gzip
        return gzip.open(file, mode)
    else:
        return open(file, mode)

def _isRun(shell, logic):
    if logic:
        shellLog = shell.replace('.sh', '.log')
        os.system('sh %s > %s' % (shell, shellLog))
    else:pass

def _mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def writeCommands(command, file, run):
    '''write commands to the file'''
    with safeOpen(file, 'w') as shell:
        shell.write(command)
    _isRun(file, run)


class NGSTools：
    
    #init
    cleanFq1 = ''
    cleanFq2 = ''
    config = configParser.ConfigParser()
    config.readfg(cfgfile)
    bam = ''
    
    def __init__(self, outdir, fq1, fq2='', quanlityBase='32'，config='~/.NGSTools.cfg'):
        self.outdir = os.path.abs(outdir)
        assert isFile(fq1)
        if fq2 != ''
            assert isfile(fq2)
        self.fq1 = fq1
        self.fq2 = fq2
        self.quanlityBase = quanlityBase

        NGSTools.cleanFq1 = fq1
        NGSTools.cleanFq1 = fq2
        if fq1.endswith('.gz'):
            NGSTools.bam = re.sub(r'f(ast)?q.gz$', 'bam', os.path.basename(fq1))
        else:
            NGSTools.bam =re.sub(r'f(ast)?q$', 'bam', os.path.basename(fq1))
        
        #read config file
        config = configParser.ConfigParser()
        config.readfg(cfgfile)
        self.genome = config.get('genome', 'fasta')
        self.bowtie2Index = config.get('genome', 'bowtie2Index')
        sefl.gtf = config.get('genome', 'gtf')
    

    def cutadapter(self, adapter5='AGATCGGAAGAGCGTCGTGTAGGGAAA', adapter3='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'，run=True):
        #cut illumina sequencing adapter
        myOutdir = self.outdir+'/qc'
        _mkdir(myOutdir)
        
        if self.fq1.endswith('.gz'):
            NGSTools.cleanFq1 = re.sub(r'f(ast)?q.gz$', 'clean.fq.gz', os.path.abs(self.fq1))
            if self.fq2 != ''
                NGSTools.cleanFq2 = re.sub(r'f(ast)?q.gz$', 'clean.fq.gz', os.path.abs(self.fq2))
        else:
            NGSTools.cleanFq1 = re.sub(r'f(ast)?q$', 'clean.fq', self.fq1)
            if self.fq2 != ''
                NGSTools.cleanFq2 = re.sub(r'f(ast)?q$', 'clean.fq', self.fq2)
        if self.fq2 != '':
            command1 = 'cutadapt -a %s -e 0.01 -m 30 -O 5 -q 15 -o %s %s' % (adapter5, NGSTools.cleanFq1, self.fq2)
        else:
            command1 = ''
        command2 = 'cutadapt -a %s -e 0.01 -m 30 -O 5 -q 10 -o %s %s' % (adapter3, NGSTools.cleanFq2, self.fq1)
    
        writeCommands(command1+'\n'+command2, myOutdir+'/cutadapter.sh', run)


    def fastqc(self, fq1=cleanFq1, fq2=cleanFq2, run=True):
        #quality control
        
        myOutdir = self.outdir+'/qc'
        _mkdir(myOutdir)
        command = 'fastqc -o %s -t 6 -d %s %s %s' % (outdir, outdir, fq1, fq2)
        writeCommands(command, myOutdir+'/fastqc.sh', run)

    def tophat2(self, fq1=cleanFq1, fq2=cleanFq2, run=True):
        '''mapping with tophat2'''
        myOutdir = self.outdir+'/mapping'
        _mkdir(myOutdir)
        
        _phredQual = ''
        if self.quanlityBase == '64'
            _phredQual = '--phred64-quals'
        command = 'tophat -p 6 -G %s %s -o %s %s %s %s' % (self.gtf, _phredQual, myOutdir, self.bowtie2Index, fq1, fq2)
        NGSTools.bam = os.path.join(myOutdir, NGSTools.bam)
        command += '\nmv %s/accepted_hits.bam %s' % (myOutdir, NGSTools.bam)
        writeCommands(command, myOutdir+'/tophat.sh', run)

    def bowtie2(self, fq1=cleanFq1, fq2=cleanFq2, mode='--local', run=True):
        '''bowtie2 mapping'''
        myOutdir = self.outdir+'/mapping'
        _mkdir(myOutdir)
        
        _phredQual = ''
        if self.quanlityBase == '64'
            _phredQual = '--phred64'
        command = 'bowtie2 %s -x %s -p 6 %s -1 %s ' % (_phredQual, self.bowtie2Index, mode, fq1)
        if fq2 != '':
            command += '-2 %s ' % fq2
        command += '| samtools view -bS -h - > %s' % NGSTools.bam
        writeCommands(command, myOutdir+'/bowtie2.sh', run)










