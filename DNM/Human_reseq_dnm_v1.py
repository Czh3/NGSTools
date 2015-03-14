#/usr/bin/env python
import argparse
import sys
import os
from string import Template

###init
REF = '/PUBLIC/database/HUMAN/genome/human/human_b37/human_g1k_v37_decoy.fasta'
ANNOVAR = "/PUBLIC/software/HUMAN/annovar"


def getChild(ped):
	'''get the chiled ID from the ped file.'''
	ids  = []
	for line in open(ped):
		words = line.split('\t')
		ids.append(words[1])
	for line in open(ped):
		words = line.split('\t')
		if words[1] in ids and words[2] in ids:
			return words[1]

def getFamilyID(ped):
	'''get family ID'''
	for i in open(ped):
		if i.startswith('#'):
			continue
		else:
			return i.split()[0]


parser = argparse.ArgumentParser(description="de novo mutaion calling module.<zhangchao@novogenen.com>", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--ped', type=str, help="The plink-liked ped file of the trio family.", required=True)
parser.add_argument('--bams', type=str, help="The trio's bam file paths divided by comma(,).", required=True)
parser.add_argument('--soft', type=str, help="The sofeware chosed to call DNM.", choices=['samtools', 'denovogear'], default='samtools')
parser.add_argument('--DNG_pp_cutoff', type=int, help="DNG Posterior probability cutoff", default=0.5)
parser.add_argument('--DNG_rd_cutoff', type=int, help="Read depth filter, sites where by DNG either one of the sample have read depth less than this threshold are filtered out.", default=10)
parser.add_argument('--region', type=str, help="Region of the BCF file to perform denovo calling.string of the form \"chr:start-end\"")
parser.add_argument('--familyID', type=str, help="Not use the familyID in ped, but use the given family ID")
parser.add_argument('--outdir', type=str, help="output directory.", default='./')

args = vars(parser.parse_args())
ped = args['ped']
bams = args['bams']
soft = args['soft']
pp_cutoff = args['DNG_pp_cutoff']
rd_cutoff = args['DNG_rd_cutoff']
region = args['region']
outdir = args['outdir']

if not os.path.exists(outdir):
	os.system('mkdir -p ' + outdir)

bams = bams.split(',')
if len(bams) != 3:
	exit('there must be 3 bams')

childID = getChild(ped)
if args['familyID']:
	familyID = args['familyID']
else:
	familyID = getFamilyID(ped)

j = 0
for bam in bams:
	for i in os.popen('samtools view -H %s' % bam):
		if i.find('SM:%s' % childID) != -1:
			del bams[j]
			bams.append(bam)
			break
	else:
		j += 1
		continue
	break
bams = ' '.join(bams)

Commands = '#!/usr/bin/bash\necho "begin\\t"`date`\n'

if soft == 'samtools':
	################	samtools 	#################
	if region:
		Commands += 'samtools mpileup -q 1 -C 50 -t DP,DV -m 2 -r %s -F 0.002 -ugf %s %s | bcftools call  -vmO v -C trio -S %s -o %s/%s.trio.vcf\n' % (region, REF, bams, ped, outdir, familyID)
	else:
		Commands += 'samtools mpileup -q 1 -C 50 -t DP,DV -m 2 -F 0.002 -ugf %s %s | bcftools call  -vmO v -C trio -S %s -o %s/%s.trio.vcf\n' % (REF, bams, ped, outdir, familyID)

	#samtools filter
	Commands += '''
	bcftools filter -s FLTER -i "%%QUAL>30 && DP>30 && MQ>40" $outdir/$familyID.trio.vcf > $outdir/$familyID.trio.tmp.vcf
	awk -F '\\t' '{if($1~/^#/)print;else if($1!~/hs37d5/ $$ $7~/PASS/ && $10~/0\/0/ && $11~/0\/0/ && $12!~/0\/0/)print}' $outdir/$familyID.trio.tmp.vcf > $outdir/$familyID.trio.vcf

	'''
else:
	##########	DNG	############
	if region:
		Commands += 'samtoolsv0.1.19 mpileup -l %s -gDSf %s %s | denovogear dnm auto --ped %s --pp_cutoff %s --rd_cutoff %s --region %s --output_vcf %s/$familyID.trio.vcf --bcf -\n' % (region, REF, bams, ped, pp_cutoff, rd_cutoff, region, outdir)
	else:
		Commands += 'samtoolsv0.1.19 mpileup -gDSf %s %s | denovogear dnm auto --ped %s --pp_cutoff %s --rd_cutoff %s --output_vcf %s/$familyID.trio.tmp.vcf --bcf -\n' % (REF, bams, ped, pp_cutoff, rd_cutoff, outdir)
	#reform DNG's vcf
	Commands += 'python /PROJ/HUMAN/share/Disease/Varition/DNM/reform_denovogear_vcf.py -i $outdir/$familyID.trio.tmp.vcf -o $outdir/$familyID.trio.vcf\n'

#annotation
Commands += '$Annovar/Var_annotation_disease.sh %s/$familyID.trio.vcf $familyID\n' % outdir

#filter
Commands += 'python /PROJ/HUMAN/share/Disease/Varition/Filter/filter_vcf.py --dbsnpNonFlagged --thausandGenomes 0.005 --function --synonymous --sift --polyphen -i %s/$familyID.trio.filter.reformated.vcf.gz -o %s/$familyID.trio\n' % (outdir, outdir)
Commands += 'vcf2tab.py -i %s/$familyID.trio.dbsnp.1000g.func.syn.sift.polyphen.vcf -o %s/$familyID.trio.dbsnp.1000g.func.syn.sift.polyphen.xls\n' % (outdir, outdir)

#annotation
Commands += '''
$Annovar/addAnn_byName.pl -annName OMIM $outdir/$familyID.trio.dbsnp.1000g.func.syn.sift.polyphen.xls \\
	| $Annovar/addAnn_byName.pl -annName GO_BP \\
	| $Annovar/addAnn_byName.pl -annName GO_CC \\
	| $Annovar/addAnn_byName.pl -annName GO_MF \\
	| $Annovar/addAnn_byName.pl -annName KEGG_PATHWAY \\
	| $Annovar/addAnn_byName.pl -annName PID_PATHWAY \\
	| $Annovar/addAnn_byName.pl -annName BIOCARTA_PATHWAY \\
	| $Annovar/addAnn_byName.pl -annName REACTOME_PATHWAY \\
	> $outdir/$familyID.trio.dbsnp.1000g.func.syn.sift.polyphen.anno.xls

echo -e 'done\\t'`date`
'''

tmp = Template(Commands)
Commands = tmp.safe_substitute(outdir=outdir, Annovar=ANNOVAR, familyID=familyID)

with open(outdir+'/%s_DNM.sh' % soft, 'w') as shell:
	shell.write(Commands)








