#!/usr/bin/env python
import argparse
import os
import sys
import re

def safe_open(file, mode='r'):
	if not os.path.exists(file):
		sys.exit('File not exists.\n')
	if file.endswith('.gz'):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)
 

parser = argparse.ArgumentParser(description="convert vcf file to tab.<zhangchao@novogene.com>")
parser.add_argument('-i','--input',help='the input VCF file(required).',required=True)
parser.add_argument('-o','--output',help='the output tab file(required).',required=True)
parser.add_argument('--info',help='output the INFO in vcf', action="store_true", default=False)
args = vars(parser.parse_args())



vcf = args['input']
tab = args['output']

title = 'CHROM   POS     ID      REF     ALT     QUAL    FILTER  GeneName        Func    Gene     GeneDetail      ExonicFunc      AAChange        Gencode cpgIslandExt    cytoBand wgRna   targetScanS     phastConsElements46way  tfbsConsSites   genomicSuperDups dgvMerged       gwasCatalog     Repeat  encodeGm12878   encodeH1hesc    encodeHelas3     encodeHepg2     encodeHuvec     encodeK562      snp138  snp138NonFlagged 1000g2012apr_eur        1000g2012apr_asn        1000g2012apr_afr        1000g2012apr_amr 1000g2012apr_all        hapmapCHB_allele        hapmapCHB_genotype       esp6500si_all   ljb23_sift      ljb23_pp2hvar   ljb23_pp2hdiv   ljb23_mt ljb23_lrt        ljb23_metalr    FORMAT'
#title = 'CHROM   POS     ID      REF     ALT     QUAL    FILTER  1000g2012apr_all  AAChange ExonicFunc  Func Gencode Gene GeneDetail GeneName Repeat cpgIslandExt  cytoBand dgvMerged encodeGm12878 encodeH1hesc encodeHelas3 encodeHepg2 encodeHuvec encodeK562 esp6500si_all genomicSuperDups gwasCatalog hapmapCHB_allele hapmapCHB_genotype ljb23_lrt ljb23_metalr ljb23_mt ljb23_pp2hdiv ljb23_pp2hvar ljb23_sift phastConsElements46way snp138 snp138NonFlagged targetScanS tfbsConsSites wgRna  FORMAT'
title = title.split()
if args['info']:
	title.insert(-1, 'INFO')
out = open(tab, 'w')
out.write('\t'.join(title))

for line in safe_open(vcf):
	line = line.strip()
	if line.startswith('##'):
		continue
	if line.startswith('#CHROM'):
		for i in line.split('\t')[9:]:
			out.write('\t' + i)
		out.write('\tshared\n')
		continue
	line = line.replace('=Name', '')
	line = line.replace('exonic;splicing', 'exonic,splicing')
	line = line.replace('UTR5;UTR3', 'UTR5,UTR3')
	line = line.replace('upstream;downstream', 'upstream,downstream')
	line = line.replace(';Name=', ',Name=')
	line = line.replace(';NM_', ',NM_')
	line = line.replace(';NR_', ',NR_')
	line = line.replace(';ALLELE_END', '')
	line = line.replace('INDEL;', '')
	line = line.replace('gff3', 'Repeat')
	#print line

	tabs = line.split('\t')
	for i in tabs[0:7]:
		out.write(i + '\t')
	for i in title:
		for j in tabs[7].split(';'):
			if not re.search(r'=', j):
				print "Warning:" + j
			j = j.split('=', 1)
			if i == j[0]:
				out.write(j[1] + '\t')
				break
			
			
	'''
	for i in tabs[7].split(';'):
		if not re.search(r'=', i):
			print "Warning:" + i
		j = i.split('=', 1)
		if j[0] in title:
			out.write(j[1] + '\t')
		else:
			pass
	'''
	if args['info']:
		info = tabs[7].split(';ANNOVAR_DATE',1)[0]
		out.write(info + '\t')
	out.write(tabs[8] + '\t')

	shared = 0
	for i in tabs[9:]:
		out.write(i + '\t')
		if not i.startswith('.'):
			shared += 1
	out.write(str(shared) + '\n')	

	
