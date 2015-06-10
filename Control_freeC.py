#!/usr/bin/env python

import argparse
from string import Template
import os.path
import os
import sys

ref='/PUBLIC/database/HUMAN/genome/human/human_b37/human_g1k_v37_decoy.fasta'

def bam2pileup(bam,outdir):
	bam_file=os.path.basename(bam)
	assert bam_file.split('.')[-1].upper() == 'BAM'
	pileup_file=bam_file.rsplit('.',1)[0]+'.pileup.gz'
	pileup_path=os.path.join(outdir,pileup_file)
	shell_file=bam_file.rsplit('.',1)[0]+'.bam2pileup.sh'
	shell_path=os.path.join(outdir,shell_file)
	code='samtools mpileup -f %s %s | gzip - > %s' % (ref, bam, pileup_path)
	# code='samtools mpileup %s | gzip - > %s' % (bam, pileup_path)
	open(shell_path,'w').write(code)
	os.system('sh %s' % shell_path)
	return pileup_path
	

parser = argparse.ArgumentParser(description='pipeline for VirusFinder2')
parser.add_argument('--loh',help="If the LOH included", type=int, choices=[0,1], default=0)
parser.add_argument('--o',help="The directory to store software output, default: cwd",default=None)
parser.add_argument('--ploidy', help="ploidy", type=int, default=2)
parser.add_argument('--sex',help="gender of the sample", choices=['XX','XY'], default=None)
parser.add_argument('--n',help="thread number default: 8",type=int,default=8)
parser.add_argument('--sample',help="sample bam file", required=True)
parser.add_argument('--control',help="control bam file", default=None)
parser.add_argument('--purity',help="purity from eg. Absolute", type=float, default=None)
parser.add_argument('--type',help="sequencing type", required=True, choices=['WGS','WES'])
parser.add_argument('--target',help="target region for capture sequencing", default=None)
parser.add_argument('--contamination',help="going to evaluate contamination", type=int, choices=[0,1], default=1)
parser.add_argument('--format',help="alignment file format", choices=['BAM','pileup'], default='BAM')

argv=parser.parse_args()
if argv.o:
	assert os.path.isdir(argv.o)
	o=argv.o
else:
	o=os.getcwd()


sample=argv.sample
control=argv.control
assert os.path.isfile(sample)
if control:
	assert os.path.isfile(control)
	flag='freec.somatic.cnv'
else:
	flag='freec.cnv'

	
loh=argv.loh
format=argv.format
if loh:
	assert os.path.isfile(control)
	if format == 'BAM':
		sample=bam2pileup(argv.sample,o)
		control=bam2pileup(argv.control,o)
		format='pileup'
	
if argv.type == 'WES':
	assert os.path.isfile(argv.target)


ploidy=argv.ploidy
sex=argv.sex
n=argv.n
type=argv.type
target=argv.target
purity=argv.purity
contamination=argv.contamination
	
config='''
[general]

BedGraphOutput=TRUE

breakPointType = 4

chrFiles = /PUBLIC/database/HUMAN/genome/human/human_b37/byChr/

chrLenFile = /PUBLIC/database/HUMAN/genome/human/human_b37/freec/chr.24.length

maxThreads= %s

outputDir = %s

ploidy = %s

# gemMappabilityFile = /PUBLIC/database/HUMAN/genome/human/human_b37/freec/human_g1k_v37_decoy.mappability.100bp.out.mappability

# uniqueMatch=TRUE

''' % (n,o,ploidy)

if sex:
	config+='''
sex= %s
''' % (sex)

if contamination:
	config+='''
contaminationAdjustment = TRUE
'''
	if purity:
		config+='''
contamination = %s
''' % (purity)



if type == 'WGS':
	config+='''
coefficientOfVariation = 0.05

#with targeted sequencing, I would not recommend to use forceGCcontentNormalization=1 or 2 since capture bias can be much stronger than GC-content bias
forceGCcontentNormalization = 1
'''

if type == 'WES':
	config+='''
#breakPointThreshold=1.5

printNA=FALSE

window = 500
'''
	if loh:
		config+='''
#will not have effect since FREEC won't use BAF information to correct predicted copy numbers
noisyData=TRUE
'''


config+='''
[sample]

mateFile = %s

inputFormat = %s

mateOrientation = FR

''' % (sample,format)

if control:
	config+='''
[control]

mateFile = %s

inputFormat = %s

mateOrientation = FR

''' % (control,format)

if loh:
	config+='''
[BAF]
#to calculate BAF values, you need to provide mateFile in SAMtools pileup format

SNPfile = /PUBLIC/database/HUMAN/control_freec/hg19_snp131.SingleDiNucl.1based.txt

minimalCoveragePerPosition = 5

minimalQualityPerPosition = 5

shiftInQuality = 33
'''

if type == 'WES':
	config+='''
[target]

captureRegions = %s

''' % (target)

config_path=os.path.join(o,'controlfreec.conf')
open(config_path,'w').write(config)
# sys.exit()

# if type == 'WGS':
	# if control:
		# if loh:
			# open(config_path,'w').write(general+general_WGS+sample+control+BAF)
		# else:
			# open(config_path,'w').write(general+general_WGS+sample+control)
	# else:
		# open(config_path,'w').write(general+gem+general_WGS+sample)

# if type == 'WES':
	# if control:
		# if loh:
			# open(config_path,'w').write(general+general_WES+sample+control+BAF)
		# else:
			# open(config_path,'w').write(general+general_WES+sample+control)
	# else:
		# open(config_path,'w').write(general+gem+general_WES+sample)


base_path=os.path.join(o,os.path.basename(sample))		
		
code='''
control_freec -conf %s
''' % (config_path)

if loh:
	code+='''
cat /PROJ/HUMAN/share/Cancer/var/makeGraph_freec.R | R --slave --args %s %s %s
''' % (ploidy, base_path+'_ratio.txt',base_path+'_BAF.txt')
else:
	code+='''
cat /PROJ/HUMAN/share/Cancer/var/makeGraph_freec.R | R --slave --args %s %s
''' % (ploidy, base_path+'_ratio.txt')

# postprocessing='''
# echo CNVs to bed format
# awk -v OFS="\t" '$1~/^([0-9]|X|Y)/{print $1,$2,$3,$3-$2,$5,"NA","NA","NA",$4,"NA"}' $bp_CNVs > $bp_CNVs.bed &&

# intersectBed -a $bp_CNVs.bed -b /PUBLIC/database/HUMAN/genome/human/human_b37/freec/all.NBlock.larger1000bp.bed -f 0.5 -v > $bp_CNVs.final.bed &&

# awk -F"\t" -v OFS="\t" 'BEGIN{id=0}{id++;print $1,"FREEC",$5,$2,$3,".",".",".","CopyNumber="$9";Size="$4";CNVID="id";CNVType="$5; print $1,"FREEC",$5,$2,$2,".",".",".","CopyNumber="$9";Size="$4";\
# CNVID="id";CNVType=breakpoint"; print $1,"FREEC",$5,$3,$3,".",".",".","CopyNumber="$9";Size="$4";CNVID="id";CNVType=breakpoint";}' %bp_CNVs.final.bed > %bp.CNVs.gff &&

# freec.cnv.summary.pl -s $s $bp_CNVs > $bp.summary.txt &&

# Var_annotation.sh -t CNVType %s/%s.%s.gff %s &&

# sv_cnv.stat.pl -s %s %s/%s.%s.hg19_multianno.xls > %s/%s.%s.stat.xls
# ''' 

# code+=Template(postprocessing).substitute()

shell_path=os.path.join(o,flag+'.'+os.path.basename(sample)+'.sh')
open(shell_path,'w').write(code)
os.system('sh %s' % (shell_path))
