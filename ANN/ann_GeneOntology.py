#!/usr/bin/env python
import sys
import os

# add gene ontology anntation tab divided files.
# zhangchao3@hotmail.com

def usage():
	print  '\nUSAGE <for cuffdiff output>:\npython %s infile outfile' % sys.argv[0]
	exit()


def safeOpen(file, mode='r'):
	if file.endswith('.gz'):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)


def loadGeneAliases():
	'''load gene aliases '''
	geneSymbol = {}
	for id in safeOpen('/home/zhangc/database/gene_aliases_20131231.gz'):
		cols = id.split('\t')
		geneSymbol[cols[1]] = cols[0]
	return geneSymbol


def loadGO():
	'''load gene ontology database'''
	GOdatabase = {}
	#key: gene
	#value: F(GOtermID):GO term|C(GOtermID):GO term


	for gene in safeOpen('/home/zhangc/database/ATH_GO_GOSLIM.txt.gz'):
		cols = gene.split('\t')

		if GOdatabase.has_key(cols[0]):
			GOdatabase[cols[0]] += '|' + cols[7] + '(' + cols[5] + '):' + cols[4]
		else:
			GOdatabase[cols[0]] = cols[7]+ '(' + cols[5] + '):' + cols[4]

	# reform database
	for key in GOdatabase:
		value = []
		for term in GOdatabase[key].split('|'):
			value.append(term)
		
		GOdatabase[key] = '|'.join(set(value))

	return GOdatabase	



def GOannotate(GOdatabase, gene):
	''' GO annotation'''

	try:
		transGene = geneSymbol[gene]
	except:
		transGene = gene

	try:
		GOannotation = GOdatabase[transGene]
	except:
		GOannotation = 'NULL'
		print 'WARNING: "%s" has no GO annotation term' % gene

	return GOannotation

if __name__ == '__main__':
	
	##Init
	# the Number of gene col  ## 0 base
	GENEIDCOL = 0

	#check input file
	if len(sys.argv) != 3:
		usage()

	infile = sys.argv[1]
	if not os.path.exists(infile):
		usage()

	# load databse for gene ID transformate
	geneSymbol = loadGeneAliases()

	# load GO database 
	GOdatabase = loadGO()

	# output file
	out = open(sys.argv[2], 'w')

	for line in safeOpen(infile):
		cols = line.strip().split('\t')
		gene = cols[GENEIDCOL]
		
		if gene == '-':
			out.write('\t'.join(cols) + '\n')
			continue


		if gene.find(',') == -1:
			# not find ','
			GOannotation = GOannotate(GOdatabase, gene)
			cols.append(GOannotation)
		else:
			for g in gene.split(','):
				GOannotation = GOannotate(GOdatabase, g)
				cols.append(GOannotation)
	
		out.write('\t'.join(cols) + '\n')
		
	out.close()
