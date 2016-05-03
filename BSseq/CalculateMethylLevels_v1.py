import sys


## calculte the methylation levels from a bed file


def readBed(bedfile):
	"""read bed file , return a dict"""

	'''
	{
		'chrapos1' : lineNumber(bed file terms id)
		'chrapos2' : lineNumber
	}
	'''

	pos = {}
	lineNum = 0
	for i in bedfile:
		j = i.strip().split("\t")
		j[1] = int(j[1])
		j[2] = int(j[2])

		for k in xrange(j[1], j[2]+1):
			ke = j[0] + 'a' + str(k)
			if pos.has_key(ke):
				pos[ke].append(lineNum)
			else:
				pos[ke] = [lineNum]

		lineNum += 1
	
	return pos


def calMethyLevels(bismark, pos, depth=3, cSites=3):
	"""calculate the methylation levels in the input bed file terms """

	##   init
	__CHR = 0
	__POS = 1
	__METHY = 3
	__UNMETHY = 4
	__TYPE = 5


	methylation = {}
	C_sites = {}
	result = {}
	for i in ["CG", "CHG", "CHH"]:
		methylation[i] = {}
		C_sites[i] = {}
		result[i] = {}
	for i in open(bismark):
		cols = i.strip().split('\t')
	
		cols[__METHY] = float(cols[__METHY])
		cols[__UNMETHY] = float(cols[__UNMETHY])

		if cols[__METHY] + cols[__UNMETHY] < int(depth):
			continue

		_key = cols[__CHR] + 'a' + cols[__POS]
		if _key in pos:
			for l in pos[_key]:
				if methylation[cols[__TYPE]].has_key(l):
					methylation[cols[__TYPE]][l] += cols[__METHY]/float(cols[__METHY] + cols[__UNMETHY])
					C_sites[cols[__TYPE]][l] += 1
				else:
					methylation[cols[__TYPE]][l] = cols[__METHY]/float(cols[__METHY] + cols[__UNMETHY])
					C_sites[cols[__TYPE]][l] = 1

	for j in ["CG", "CHG", "CHH"]:
		for i in methylation[j]:
			if C_sites[j][i] >= int(cSites):
				result[j][i] = methylation[j][i]/float(C_sites[j][i])
			else:
				result[j][i] = 'NULL'
	return result



if __name__ == '__main__':
	usage = 'python CalculateMethylLevels.py bismarkFile inputBed  depth minCsites outfile'
	if len(sys.argv) != 6:
		sys.exit(usage)	


	bedfile = open(sys.argv[2]).readlines()
	data = readBed(bedfile)
	result = calMethyLevels(sys.argv[1], data, sys.argv[3], sys.argv[4])


	out = open(sys.argv[5], 'w')

	'''
	for j in ["CG", "CHG", "CHH"]:
		# check
		for i in sorted(result[j]):
			for k in ["CG", "CHG", "CHH"]:
				try:
					assert result[k][i]
				except:
					result[k][i] = 0
	'''

	for j in ["CG", "CHG", "CHH"]:
		#for i in sorted(result[j]):
		for i in xrange(len(bedfile)):
			try:
				bedfile[i] = bedfile[i].strip() + '\t' + str(result[j][i]) + '\n'
			except:
				bedfile[i] = bedfile[i].strip() + '\t0\n'
		#out.write(bedfile[i].strip() + res + '\n')
	out.writelines(bedfile)
	out.close()



