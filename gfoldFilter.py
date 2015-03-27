#!/usr/bin/env python

#
#	A script for GFOLD output file filter to get Different Expression Genes(DEGs).
#	gfold version: V1.1.1
#


import sys
import os

__EMAIL__ = 'zhangchao3@hotmail.com'

def usage():
	''' print usage information '''
	print '''
USAGE:
python gfoldFilter.py gfold_output_file  DEGs_output_prefix
'''

def GfoldFilter(file, outfilePrefix='out', rpkm=5, log2fdc=2):
	''' filter the GFOLD output file to get the most-likely different expression genes
		file: the GFOLD output file
		rpmk: when bother condition are lower than this rpkm will be discard
		log2fdc: log2(fold change) filter.
	'''
	#init
	upRegulateGenesOutfile = open(outfilePrefix+'_up.DEGs', 'w')
	downRegulateGenesOutfile = open(outfilePrefix+'_down.DEGs', 'w')
	upRegulateGeneList = open(outfilePrefix+'_up.DEGs.genelist.txt', 'w')
	downRegulateGeneList = open(outfilePrefix+'_down.DEGs.genelist.txt', 'w')

	for gene in open(file):
		if gene.startswith('#') or gene.startswith('\t'):
			continue

		column = gene.strip().split('\t')

		if column[2] == '0':
			continue

		# filter rpkm
		if float(column[5]) < rpkm and float(column[6]) < rpkm:
			continue

		# filter log2fdc
		if float(column[4]) > log2fdc:
			upRegulateGenesOutfile.write(gene)
			upRegulateGeneList.write(column[0]+'\n')
		elif float(column[4]) < -log2fdc:
			downRegulateGenesOutfile.write(gene)
			downRegulateGeneList.write(column[0]+'\n')

	upRegulateGenesOutfile.close()
	downRegulateGenesOutfile.close()


if __name__ == '__main__':
	if len(sys.argv) == 3 and os.path.exists(sys.argv[1]):
		GfoldFilter(sys.argv[1], sys.argv[2], 1, 2)
	else:
		sys.exit(usage())



