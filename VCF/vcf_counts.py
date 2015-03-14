#!/usr/bin/env python 
import os
from vcf import VCF
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="vcf counts", formatter_class=RawTextHelpFormatter)
parser.add_argument('--list', help="countlist file, which contain several vcf file path.")
parser.add_argument('--odir', help='out dir')
argv = vars(parser.parse_args())

files = argv['list']

if argv['odir']:
	odir = argv['odir']
else:
	odir = 'counts'

if not os.path.exists(odir):
	assert not os.mkdir(odir)

#chromosome = open(odir+'/chr.xls', 'w')
snp_func = open(odir+'/snp_func.xls', 'w')
snp_other = open(odir+'/snp_other.xls', 'w')
snp_exonicfunc = open(odir+'/snp_exonicfunc.xls', 'w')
indel_func = open(odir+'/indel_func.xls', 'w')
indel_exonicfunc = open(odir+'/indel_exonicfunc.xls', 'w')
indel_other = open(odir+'/indel_other.xls', 'w')
indel_length = open(odir+'/indel_length.xls', 'w')
ts_tv = open(odir+'/ts_tv.xls', 'w')


def mydivide(x,y):
	if y == 0:
		return float(0)
	else:
		return float(x)/float(y)


##init
flag = 1
indel_len = {}


############title
##chr.xls
#title = [str(i) for i in range(1,23)]
#title_chr = ['Sample'] + title + ['X','Y']
#chromosome.write('\t'.join(title_chr)+'\n')


##snp_func.xls
title_snp_func = ['Sample','exonic','intronic','UTR3','UTR5','intergenic','ncRNA_exonic','ncRNA_intronic','upstream','downstream','splicing','ncRNA_UTR3','ncRNA_UTR5','ncRNA_splicing']
snp_func.write('\t'.join(title_snp_func)+'\n')


##snp_exonicfunc.xls
title_snp_exonicfunc = ['Sample','synonymous_SNV','missense_SNV','stopgain','stoploss','unknown']
snp_exonicfunc.write('\t'.join(title_snp_exonicfunc)+'\n')

##indel_func.xls
title_indel_func = ['Sample','exonic','intronic','UTR3','UTR5','intergenic','ncRNA_exonic','ncRNA_intronic','upstream','downstream','splicing','ncRNA_UTR3','ncRNA_UTR5','ncRNA_splicing']
indel_func.write('\t'.join(title_indel_func)+'\n')

##indel_exonicfunc.xls
title_indel_exonicfunc = ['Sample','frameshift_deletion','frameshift_insertion','nonframeshift_deletion','nonframeshift_insertion','stoploss','stopgain','unknown']
indel_exonicfunc.write('\t'.join(title_indel_exonicfunc)+'\n')






for file in open(files, 'r'):
	if file.startswith('#'):continue
	file = file.strip()
	sample_name = os.path.basename(file)
	sample_name = sample_name.split('.')[0]

	myVCF = VCF(file)
	snp = myVCF.filter()
	indel_file = file.replace('snp','indel')
	myVCF = VCF(indel_file)
	indel = myVCF.filter()
	
	"""
	##chr.xls
	chr = myVCF.chr_stat(vcf)
	chromosome.write(sample_name)
	for i in [str(i) for i in range(1,23)]+['X','Y']:
		try:
			chromosome.write('\t'+str(chr[i]))
		except:
			chromosome.write('\t0')
	chromosome.write('\n')
	"""

	#snp_func.xls
	snp_f = myVCF.stat_func(snp)
	snp_func.write(sample_name)
	for i in title_snp_func[1:]:
		if i in snp_f:
			snp_func.write('\t'+str(snp_f[i]))
		else:
			snp_func.write('\t0')
	snp_func.write('\n')

	#snp_exonicfunc.xls
	snp_f = myVCF.stat_exonicfunc(snp)
	snp_exonicfunc.write(sample_name)
	for i in title_snp_exonicfunc[1:]:
		if i in snp_f:
			snp_exonicfunc.write('\t'+str(snp_f[i]))
		else:
			snp_exonicfunc.write('\t0')
	snp_exonicfunc.write('\n')

	#snp_other.xls
	Hom, Het = myVCF.get_Hom_Het(snp)
	novo = myVCF.novo(snp)
	if flag:
		snp_other.write('Sample\tall\tgenotype.Het\tgenotype.Hom\tnovel\tnovel_proportion\n')
	snp_other.writelines('\t'.join([sample_name,str(len(snp)),str(Het),str(Hom),str(len(novo)),'%.9f' % (len(novo)/float(len(snp)))])+'\n')

	#indel_func.xls
	indel_f = myVCF.stat_func(indel)
	indel_func.write(sample_name)
	for i in title_indel_func[1:]:
		if i in indel_f:
			indel_func.write('\t'+str(indel_f[i]))
		else:
			indel_func.write('\t0')
	indel_func.write('\n')

	#indel_exonicfunc.xls
	indel_f = myVCF.stat_exonicfunc(indel)
	indel_exonicfunc.write(sample_name)
	for i in title_indel_exonicfunc[1:]:
		if i in indel_f:
			indel_exonicfunc.write('\t'+str(indel_f[i]))
		else:
			indel_exonicfunc.write('\t0')
	indel_exonicfunc.write('\n')

	#indel_other.xls
	Hom, Het = myVCF.get_Hom_Het(indel)
	novo1 = myVCF.novo(indel)
	if flag:
		indel_other.write('Sample\tall\tgenotype.Het\tgenotype.Hom\tnovel\tnovel_proportion\n')
	indel_other.writelines('\t'.join([sample_name,str(len(indel)),str(Het),str(Hom),str(len(novo1)),'%.9f' % (mydivide(len(novo1), float(len(indel))))])+'\n')

	#indel_length.xls
	indel_len[sample_name] = myVCF.stat_indel_length(indel)

	#ts_tv.xls
	ts, tv, tstv = myVCF.ts_tv(snp)
	ts1, tv1, tstv1 = myVCF.ts_tv(novo)
	if flag:
		ts_tv.writelines('Sample\tnovel_ts\tnovel_ts/tv\tnovel_tv\tts\tts/tv\ttv\n')
	ts_tv.writelines('\t'.join([sample_name,str(ts1),str(tstv1),str(tv1),str(ts),str(tstv),str(tv)])+'\n')

	flag = 0

# deal indel_length.xls
merge_indel = []
for i in indel_len:
	for j in indel_len[i]:
		if j not in merge_indel:
			merge_indel.append(j)

indel_length.write('Sample')
merge_indel.sort()
for i in merge_indel:
	indel_length.write('\t'+str(i))
indel_length.write('\n')


for i in indel_len:
	indel_length.write(i)
	for site in merge_indel:
		if site in indel_len[i]:
			indel_length.write('\t'+str(indel_len[i][site]))
		else:
			indel_length.write('\t0')
	indel_length.write('\n')			
