import sys
from scipy import stats


global GO
GO = {}
'''
{
	GO_ID1 : [gene1, gene2],
	GO_ID2 : ...
}
'''

global totalGeneNum
totalGene = []

for i in open("NCIB8687.gene_go.txt"):
	j = i.strip().split()
	
	if GO.has_key(j[1]):
		GO[j[1]].append(j[0])
	else:
		GO[j[1]] = [j[0]]

	totalGene.append(j[1])

geneList = []
for i in open(sys.argv[1]):
	geneList.append(i.strip())

genelist = list(set(geneList))


totalGeneNum = len(set(totalGene))

def testEnrichment(geneList, GOterm):
	geneNum = len(geneList)
	inTermGene = []
	otherTermGene = []

	for i in GO[GOterm]:
		if i in geneList:
			inTermGene.append(i)

	inTermNum = len(set(inTermGene))
	outTermNum = geneNum - inTermNum

	termGeneNum = len(GO[GOterm])
	termGeneNotInList = termGeneNum - inTermNum

	otherTermGeneNotInList = totalGeneNum - termGeneNum

	oddsratio, pvalue = stats.fisher_exact([[inTermNum, outTermNum], [termGeneNotInList, otherTermGeneNotInList]])

	if pvalue < 0.05 and inTermNum >= 2:
		print "\t".join([GOterm, str(inTermNum), "%.2f" % (100.0*inTermNum/termGeneNum), ",".join(inTermGene), str(pvalue)])


# title
print "GO_Term\tConut\t%\tGenes\tPvalue"

for i in GO:
	testEnrichment(geneList, i)




