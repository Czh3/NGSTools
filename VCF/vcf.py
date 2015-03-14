import sys
import os
import re

class VCF:
	"""some function to deal with VCF file"""

 	def __init__(self, file):
	#	super(VCF, self).__init__()
		self.file = file

	def _readVCF(self):

		"""read vcf file"""
		if not os.path.exists(self.file):
			sys.exit('%s file not exist\n' % self.file)

		if self.file.endswith('.gz'):
			import gzip
			return gzip.open(self.file, 'r').readlines()
		else:
			return open(self.file, 'r').readlines()

	def filter(self,stat=''):
		"""filter the VCF"""
		vcf = self._readVCF()

		filter_vcf = []
		for i in vcf:
			if i.startswith('#'):
				continue
			words = i.split('\t')
			if stat == 'PASS':
				if words[6] != 'PASS' or words[6] != '.':
					continue

			if words[4].find(','):
				words[4] = words[4].split(',')[0]
			filter_vcf.append(i)
				

		vcf = []	
		return filter_vcf

 	def novo(self, vcf):
		"""get the novo muation"""
		novo = []

		for i in range(len(vcf)):
			words = vcf[i].split('\t')
			if not re.search('snp\d+=rs',words[7]):
				novo.append(vcf[i])
				continue

		return novo

	def get_SNP_INDEL(self, vcf):
		"""get the SNP and InDel muation\n """

		snp = []
		indel = []
		for i in range(len(vcf)):
			words = vcf[i].split('\t')
			if len(words[3]) == 1 and len(words[4].split(',')[0]) == 1:
				snp.append(vcf[i])
			else:
				indel.append(vcf[i])
		return snp, indel

	def chr_stat(self, vcf):
		"""stat the muation in the chromosomo\n
		return a dict:
		chromosomo = {
			'1' : numOfMuationInTheChr1
			'2'	: numOfMuationInTheChr2
			...
		}
		"""

		chromosomo = {}
		for line in vcf:
			words = line.split('\t')
			if words[0] in chromosomo:
				chromosomo[str(words[0])] += 1
			else:
				chromosomo[str(words[0])] = 1

		return chromosomo

	def stat_indel_length(self, vcf):
		"""stat the number per indel length
		return a dict;
		indel_length = {
			1 : number
			2	: number
			-1 : number
			...
		}
		"""
		indel_length = {}

		for line in vcf:
			words = line.split('\t')
			if len(words[3]) == 1 and len(words[4]) == 1:
				continue
			else:
				length = len(words[4]) - len(words[3])
				if length in indel_length:
					indel_length[length] += 1
				else:
					indel_length[length] = 1
		return indel_length


	def stat_func(self, vcf):
		"""stat the num of snp or indel muation that functional
		return a dict:
		func = {
			'intergenic' : number
			'exonic' : number
			...	
		}
		"""
		func = {}
		function = ''
		for line in vcf:
			words = line.split('\t')
			for i in words[7].split(';'):
				if i.startswith('Func=') :#Func=ncRNA_exonic
					function = i.split('=')[1]
					break

			for f in function.split(','):
				if f in func:
					func[f] += 1
				else:
					func[f] = 1
		return func

	def stat_exonicfunc(self, vcf):
		"""stat the num of snp/indel muation that  exonic functional
		return a dict:
		exonicfunc = {
			'missense_SNV' : number
			...	
		}
		"""
		exonicfunc = {}
		function = ''
		flag = 0

		for line in vcf:
			words = line.split('\t')
			for i in words[7].split(';'):
				if i.startswith('Func='):
					function = i.split('=')[1]
					for j in function.split(','):
						if j == 'exonic':
							flag = 1
							break
			if flag == 0:
				continue
			for i in words[7].split(';'):
				if i.startswith('ExonicFunc='):
					function = i.split('=')[1]
					break

			for f in function.split(','):
				if f in exonicfunc:
					exonicfunc[f] += 1
				else:
					exonicfunc[f] = 1
			flag = 0
		return exonicfunc

	def get_Hom_Het(self, vcf):
		"""get the number of Het and Hom muation
		return a array[2]
		[Hom_num,Het_num]
		"""

		Het_num = 0
		Hom_num = 0	
		for line in vcf:
			hehe = line.split('\t')[-2].split(':')
			xixi = line.split('\t')[-1].split(':')
			format = dict(zip(hehe,xixi))
			if not format.has_key('GT'):
				continue
			if format['GT'] == '1/1' or format['GT'] == '0/0':
				Hom_num += 1
			else:
				Het_num += 1
		return Hom_num, Het_num

	def ts_tv(self, vcf):
		"""for snp muation
		get the ts/tv ratio
		return:
		ts,tv,ts/tv
		"""
		ts = 0
		tv = 0
		ts_flag = ['AG', 'GA', 'CT', 'TC']

		for line in vcf:
			words = line.split('\t')
			if words[3].upper()+words[4].split(',')[0].upper() in ts_flag:
				ts += 1
			else:
				tv += 1
		return ts, tv, '%.9f' % (ts/float(tv))
 		
