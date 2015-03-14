#!/usr/bin/env python
import argparse
import re
import os

def safe_open(file, mode='r'):
	'''safe open '''
	try:
		if not file.endswith('.gz'):
			return open(file, mode)
		else:
			import gzip
			return gzip.open(file, mode)
	except IOError:
		print file + " not found."

def dbsnp(vcf, outvcf):
	'''filter the variarions that in dbSNP'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if not re.search(r';snp\d+=rs\d+', line):
			number += 1
			out.write(line)
	return str(number)

def dbsnpNonFlagged(vcf, outvcf):
	'''filter the variarions that in dbsnpNonFlagged'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if not re.search(r';snp\d+NonFlagged=rs\d+', line):
			number += 1
			out.write(line)
	return str(number)

def thausandGenomes(vcf, outvcf, frequence):
	'''filter the frequence from 1000 genomes bigger than given.'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if re.search(r'1000g2012apr_all=[0-9|\.]+', line):
			fre = re.search(r'1000g2012apr_all=([0-9|\.]+)', line).group(1)
			if fre <= frequence:
				number += 1
				out.write(line)
		else:
			number += 1
			out.write(line)
	return str(number)

def function(vcf, outvcf):
	'''filter functions are not exonic or splicing.'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if re.search(r';Func=', line):
			fun = re.search(r';Func=([\w|\d]+);', line).group(1)
			if ('exonic' == fun) or ('splicing' == fun):
				number += 1
				out.write(line)
		#else:
		#	number += 1
		#	out.write(line)
	return str(number)

def synonymous(vcf, outvcf):
	'''filter the synonymous variarions.'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if not re.search(r'ExonicFunc=synonymous_SNV', line):
			number += 1
			out.write(line)
	return str(number)

def sift(vcf, outvcf):
	'''filter the tolerated variations that sift predicts.'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if not re.search(r'_sift=[0-9|\.|,]+T', line):
			number += 1
			out.write(line)
	return str(number)

def polyphen(vcf, outvcf):
	'''filter the benign variations that polyphen2 predicts'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if not ( re.search(r'_pp2hvar=[0-9|\.|,]+B', line) and re.search(r'_pp2hdiv=[0-9|\.|,]+B', line) ):
			number += 1
			out.write(line)
	return str(number)

def mutationTaster(vcf, outvcf):
	'''filter the polymorphism variations that polyphen2 predicts.'''
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf, 'r'):
		if line.startswith('#'):
			out.write(line)
			continue
		if not re.search(r'_mt=[0-9|\.|,]+[P|N]', line):
			number += 1
			out.write(line)
	return str(number)

def caseOrControl(vcf, outvcf, sample, caseNumber, mode):
	"""choose the case or control sample"""
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	_sample = []
	_position = []

	if os.path.exists(sample):
		for s in safe_open(sample):
			_sample.append(s.strip())
	else:
		for s in sample.split(','):
			_sample.append(s)

	for line in safe_open(vcf, 'r'):
		line = line.strip()
		flag = 0
		if line.startswith('##'):
			out.write(line + '\n')
			continue

		if line.startswith('#CHROM'):
			try:
				_samples = line.split('\t')[9:]
				for i in range(len(_samples)):
					if _samples[i] in _sample:
						_position.append(i)
				continue
			except:
				print "Error 1: your vcf file must contains this '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT' line."
				sys.exit()
			finally:
				out.write(line + '\n')

		_samples = line.split('\t')[9:]
		
		for i in range(len(_samples)):
			if i in _position:
				if mode == 'case':
					if _samples[i] != '.':
						flag += 1
				elif mode == 'control':
					if _samples[i] == '.':
						flag += 1
				else:
					sys.exit()
			
		if flag >= int(caseNumber):
			out.write(line + '\n')
			number += 1
		
	return str(number)

def getGeneList(vcfLine):
	"""get gene from a singe vcf line"""
	_res = []
	mat = re.search(r"GeneName=Name=(\S+);Func=", vcfLine)
	if not mat:
		return _res
	for i in mat.group(1).split(','):
		_res.append(i)
	return _res

def outputGeneVCF(vcf, outvcf, geneList):
	"""output a vcf contains the given gene list"""
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	for line in safe_open(vcf):
		flag = 0
		if line.startswith('#'):
			out.write(line)
			continue
		for i in getGeneList(line):
			if i in geneList:
				flag = 1
		if flag:
			out.write(line)
			number += 1
	return str(number)

def dominant(vcf, outvcf, sample, sampleNumber):
	"""though the dominant gentic mode of this disease to get the more likely muations."""
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()
	
	_sample = []
	_position = []
	
	if os.path.exists(sample):
		for s in safe_open(sample):
			_sample.append(s.strip())
	else:
		for s in sample.split(','):
			_sample.append(s)

	for line in safe_open(vcf):
		line = line.strip()
		if line.startswith('##'):
			out.write(line + '\n')
			continue
		if line.startswith('#CHROM'):
			try:
				_samples = line.split('\t')[9:]
				for i in range(len(_samples)):
					if _samples[i] in _sample:
						_position.append(i)
				continue
			except:
				print "Error 1: your vcf file must contains this '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT' line."
				sys.exit()
			finally:
				out.write(line + '\n')

		_samples = line.split('\t')[9:]
		_isHet = 0

		###choose the heterozygous mutations
		for i in range(len(_samples)):
			if i in _position:
				if _samples[i].split(':')[0] != '1/1':
					_isHet += 1
		if _isHet >= int(sampleNumber):
			out.write(line + '\n')
			number += 1
	
	return str(number)


def recessive(vcf, outvcf, sample, sampleNumber):
	"""though the recessive gentic mode of this disease to get the more likely muations."""
	number = 0

	if vcf.find('.vcf'):
		out = open(outvcf, 'w')
	else:
		sys.exit()

	_geneList =[]
	_doubleGene = []
	_sample = []
	_position = []

	if os.path.exists(sample):
		for s in safe_open(sample):
			_sample.append(s.strip())
	else:
		for s in sample.split(','):
			_sample.append(s)

	for line in safe_open(vcf, 'r'):
		line = line.strip()
		if line.startswith('##'):
			out.write(line + '\n')
			continue

		if line.startswith('#CHROM'):
			try:
				_samples = line.split('\t')[9:]
				for i in range(len(_samples)):
					if _samples[i] in _sample:
						_position.append(i)
				continue
			except:
				print "Error 1: your vcf file must contains this '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT' line."
				sys.exit()
			finally:
				out.write(line + '\n')

		_samples = line.split('\t')[9:]
		_isHom = 0

		### choose the homozygous mutation
		for i in range(len(_samples)):
			if i in _position:
				if _samples[i].split(':')[0] == '1/1':
					_isHom += 1
		if _isHom >= int(sampleNumber):
			out.write(line + '\n')
			number += 1

		for i in getGeneList(line):
			if i in _geneList:
				_doubleGene.append(i)
			else:
				_geneList.append(i)

	return str(number), set(_geneList), set(_doubleGene)


		


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Filter VCF file to get rare muations which may associate with disease.\n<zhangchao@novogene.com>")
	parser.add_argument('-i','--input',help="the input vcf file.(.gz file is acceptable)", required=True)
	parser.add_argument('-o','--output',help="output file prefix.This program will output several vcf file", default='out')
	parser.add_argument('--dbsnp',help="filter the variations that are in dbSNP database.",action="store_true", default=False)
	parser.add_argument('--dbsnpNonFlagged',help="filter the variations that are in dbsnpNonFlagged database.",action="store_true", default=False)
	parser.add_argument('--thausandGenomes',help="filter the variations whose frequence (based on 1000genomes project) are high. \nDefault frequence to filter :0.005\nSet 'NULL' to not filter the 1000gemomes.", default='0.005')
	parser.add_argument('--function',help="filter the variations whose functions are not exonic or splicing",action="store_true", default=False)
	parser.add_argument('--synonymous',help="filter the synonymous variarions",action="store_true", default=False)
	parser.add_argument('--sift',help="filter the tolerated variarions that sift predicts",action="store_true", default=False)
	parser.add_argument('--polyphen',help="filter the benign variarions that polyphen predicts",action="store_true", default=False)
	parser.add_argument('--mutationTaster',help="filter the polymorphism and polymorphism_automatic variarions that Mutation Taster predicts",action="store_true", default=False)
	parser.add_argument('--geneList',help="get the variarions in this gene list.")
	parser.add_argument('--case',help="the case sample name divided by comma(,). same as the name in merged VCF.\n"
									"(can be a file contains sample name in each line)")
	parser.add_argument('--caseNumber',help="when the case sample number >= this argument will be output.")
	parser.add_argument('--control',help="the control sample name divided by comma(,). same as the name in merged VCF.\n"
									"(can be a file contains sample name in each line)")
	parser.add_argument('--controlNumber',help="when the control sample number >= this argument will be output.")
	parser.add_argument('--genticModel',help="though the gentic mode of this disease to get the more likely muations\n"
									"can be choosed:	dominant or recessive\n"
									"\tThe dominant choice will get the heterozygous mutations.\n"
									"\tThe recessive choice will get the homozygous mutations or more than 2 heterogeneous mutations in one gene.",choices=['dominant','recessive'])

	args = vars(parser.parse_args())
	vcf = args['input'].strip()
	outvcf = args['output']
	stat = open(outvcf+'.filter.stat.xls', 'w')
	outvcf += '.vcf'
	title = 'total\t'
	number = os.popen('less %s | grep -v "#" | wc -l' % vcf).read().split()[0] + '\t'
	
	if args['dbsnp']:
		outvcf = outvcf.replace('.vcf', '.dbsnp.vcf')
		number += dbsnp(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'dbsnp\t'
	if args['dbsnpNonFlagged']:
		outvcf = outvcf.replace('.vcf', '.dbsnp.vcf')
		number += dbsnpNonFlagged(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'dbsnpNonFlagged\t'
	if args['thausandGenomes'] != 'NULL':
		outvcf = outvcf.replace('.vcf', '.1000g.vcf')
		number += thausandGenomes(vcf, outvcf, args['thausandGenomes']) + '\t'
		vcf = outvcf
		title += '1000G\t'
	if args['function']:
		outvcf = outvcf.replace('.vcf', '.func.vcf')
		number += function(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'function\t'
	if args['synonymous']:
		outvcf = outvcf.replace('.vcf', '.syn.vcf')
		number += synonymous(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'synonymous\t'
	if args['sift']:
		outvcf = outvcf.replace('.vcf', '.sift.vcf')
		number += sift(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'sift\t'
	if args['polyphen']:
		outvcf = outvcf.replace('.vcf', '.polyphen.vcf')
		number += polyphen(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'polyphen\t'
	if args['mutationTaster']:
		outvcf = outvcf.replace('.vcf', '.mutationTaster.vcf')
		number += mutationTaster(vcf, outvcf) + '\t'
		vcf = outvcf
		title += 'mutationTaster\t'
	if args['geneList']:
		outvcf = outvcf.replace('.vcf', '.geneList.vcf')
		myGeneList = []
		if os.path.exists(args['geneList']):
			for i in safe_open(args['geneList']):
				myGeneList.append(i.strip())
		else:
			for i in args['geneList'].split(','):
				myGeneList.append(i)
		number += outputGeneVCF(vcf, outvcf, myGeneList) + '\t'
		vcf = outvcf
		title += 'geneList\t'
	if args['case'] and args['caseNumber']:
		outvcf = outvcf.replace('.vcf', '.case.vcf')
		number += caseOrControl(vcf, outvcf, args['case'], args['caseNumber'], 'case') + '\t'
		vcf = outvcf
		title += 'case\t'
	if args['control'] and args['controlNumber']:
		outvcf = outvcf.replace('.vcf', '.control.vcf')
		number += caseOrControl(vcf, outvcf, args['control'], args['controlNumber'], 'control') + '\t'
		vcf = outvcf
		title += 'control\t'
	if args['genticModel'] == 'dominant':
		if not args['case']:
			sys.exit('genticModel argument must with the case argument.\n')
		outvcf = outvcf.replace('.vcf', '.heterozygous.vcf')
		number += dominant(vcf, outvcf, args['case'], args['caseNumber']) + '\t'
		vcf = outvcf
		title += 'heterozygous\t'
	if args['genticModel'] == 'recessive':
		if not args['case']:
			sys.exit('genticModel argument must with the case argument.\n')
		outvcf = outvcf.replace('.vcf', '.homozygous.vcf')
		a, b, c= recessive(vcf, outvcf, args['case'], args['caseNumber'])
		number += a + '\t'
		title += 'homozygous\t'

		outGeneList = outvcf.replace('.homozygous.vcf', '.homozygous.geneList')
		gene = safe_open(outGeneList, 'w')
		for i in b:
			gene.write(i + '\n')
		gene.close()

		outvcf = outvcf.replace('homozygous', 'compoundHeterozygous')
		number += outputGeneVCF(vcf, outvcf, c) + '\t'
		title += 'compoundHeterozygous\t'

		outGeneList = outvcf.replace('.compoundHeterozygous.vcf', '.compoundHeterozygous.geneList')
		gene = safe_open(outGeneList, 'w')
		for i in c:
			gene.write(i + '\n')
		gene.close()




	stat.write(title.strip('\t') + '\n')
	stat.write(number.strip('\t') + '\n')
	stat.close()
	

