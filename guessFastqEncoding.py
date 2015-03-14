#!/usr/bin/env python
#-*- encoding=utf8 -*-

import sys
from re import search

def usage():
	usage = '\tUSAGE:\n\tpython %s YourFastqFile' % sys.argv[0]
	sys.exit(usage)


def safe_open(file, mode='r'):
	if file.endswith('.gz'):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)


def main():
	if len(sys.argv) < 2:
		usage()

	defaultLineNumber = 0
	for line in open(sys.argv[1]):
		defaultLineNumber += 1
		if defaultLineNumber % 4 == 0:
			if search(r'\d', line):
				print '\n*** encoding: phred+33 ***'
				break

		if defaultLineNumber > 10000:
			print '\n*** encoding: phred+64 ***'


if __name__ == '__main__':
	main()
