#!/usr/bin/env python

import gzip
import sys


def load_gtf(file):
	'''load genome gtf file '''
	refflat = {}
	for line in gzip.open(file):
		cols = line.split('\t')
		refflat[cols[1]] = "\t".join([cols[2], cols[4], cols[5]])
	return refflat

def annGene(annFile, chr, pos):
	for gene in annFile:
		value = annFile[gene]
		info = value.split('\t')
		if chr == info[0] and int(info[1]) < int(pos) and int(info[2]) > int(pos):
			return gene

	return 'intergenetic'


if __name__ == '__main__':

	annFile = load_gtf("/home/zhangc/reference/Arabidopsis_thaliana/NCBI/TAIR10/Annotation/Genes/refFlat.txt.gz")

	for fusion in open(sys.argv[1]):
		fusion = fusion.strip()
		cols = fusion.split('\t')
		(chr1, chr2) = cols[0].split('-')
		pos1 = cols[1]
		pos2 = cols[2]
		gene1 = annGene(annFile, chr1, pos1)
		gene2 = annGene(annFile, chr2, pos2)
		print gene1+'\t'+gene2+'\t'+fusion

