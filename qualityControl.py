#!/usr/bin/env python
#-*- encoding=utf8 -*-

from Bio import SeqIO
import sys

import argparse


def safe_open(file, mode='r'):
	''' safe open file '''

	if file.endswith('.gz'):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)


def average(List):
	'''calculate the mean of the input list'''

	return sum(List)/float(len(List))


def qualityLessThan(List, quality):
	'''count the bases which quality less than the given number.'''

	lowQualityNumber = 0
	for baseQual in List:
		if baseQual <= quality:
			lowQualityNumber += 1

	return lowQualityNumber


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Illumina Single-end or Pair-end sequencing quality control tools.<zhangchao3@hotmail.com>',)

	parser.add_argument('-1', '--fastq1',
						help='The input fastq1 file.',
						required=True)

	parser.add_argument('-2', '--fastq2',
						help='The input fastq2 file.')

	parser.add_argument('-q', '--minQuality',
						help='Minimum quality score to keep.',
						type=int,
						default=5)

	parser.add_argument('-p', '--minPercent',
						help='Minimum percent of bases that must have [-q] quality.',
						type=int,
						default=50)

	parser.add_argument('-a', '--averageQuality',
						help='The average quality of SE/PE reads less than the averageQuality will be descard.',
						type=float,
						default=10)

	parser.add_argument('-f', '--phredQualtiy',
						help='phred quality score',
						choices=[33, 64],
						default=33)

	parser.add_argument('-o1', '--outFile1',
						help='fastq1 output file.',
						required=True)

	parser.add_argument('-o2', '--outFile2',
						help='fastq2 output file.')

	args = parser.parse_args()

	if (args.outFile2 and (not args.fastq2)) or ((not args.outFile2) and args.fastq2):
		sys.exit('Check your input arguments')
		

	fastqFormat = 'fastq'
	if args.phredQualtiy == 64:
		fastqFormat = 'fastq-illumina'


	

	if not args.fastq2:
		# SE reads
		handle = safe_open(args.fastq1, 'rU')
		fastq1out = safe_open(args.outFile1, 'wb')

		for record in SeqIO.parse(handle, fastqFormat):
			phredQualityList = record.letter_annotations['phred_quality']

			if average(phredQualityList) < args.averageQuality:
				# filter average quality score is to low
				continue

			if qualityLessThan(phredQualityList, args.minQuality) > args.minPercent * len(phredQualityList) / 100.0:
				# filter the low quality base is to much
				continue

			fastq1out.write(record.format(fastqFormat))

		fastq1out.close()

	else:
		# PE reads
		fastq1handle = safe_open(args.fastq1, 'rU')
		fastq2handle = safe_open(args.fastq2, 'rU')
		fastq1out = safe_open(args.outFile1, 'wb')
		fastq2out = safe_open(args.outFile2, 'wb')

		fastq1SeqIO = SeqIO.parse(fastq1handle, fastqFormat)
		fastq2SeqIO = SeqIO.parse(fastq2handle, fastqFormat)
		while True:
			try:
				fastq1record = fastq1SeqIO.next()
				fastq2record = fastq2SeqIO.next()
			except:
				break
			

			phredQualityList1 = fastq1record.letter_annotations['phred_quality']
			phredQualityList2 = fastq2record.letter_annotations['phred_quality']

			if (average(phredQualityList1) < args.averageQuality) or (average(phredQualityList2) < args.averageQuality):
				continue


			if qualityLessThan(phredQualityList1, args.minQuality) > args.minPercent * len(phredQualityList1) / 100.0:
				if qualityLessThan(phredQualityList2, args.minQuality) > args.minPercent * len(phredQualityList2) / 100.0:
					continue

			fastq1out.write(fastq1record.format(fastqFormat))
			fastq2out.write(fastq2record.format(fastqFormat))

		fastq1out.close()
		fastq2out.close()












