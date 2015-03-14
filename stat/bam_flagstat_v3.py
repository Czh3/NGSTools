#!/PUBLIC/software/public/System/Python-2.7.6/bin/python
import HTSeq
import argparse
import re
import os

parser = argparse.ArgumentParser(description="samtools flagstat statistics")
parser.add_argument('--bam',help="the input bam file",type=str)
#parser.add_argument('--logdir',help="the path for .e & .o file")
#parser.add_argument('--sample',help="sample name")
args = parser.parse_args()


#logdir = args.logdir
bam = args.bam
#sample = args.sample


flags = {0x1 : "template having multiple segments in sequencing",
         0x2 : "each segment properly aligned according to the aligner",
         0x4 : "segment unmapped",
         0x8 : "next segment in the template unmapped",
         0x10 : "SEQ being reverse complemented",
         0x20 : "SEQ of the next segment in the template being reversed",
         0x40 : "the first segment in the template",
         0x80 : "the last segment in the template",
         0x100 : "secondary alignment",
         0x200 : "not passing quality controls",
         0x400 : "PCR or optical duplicate",
         0x800 : "supplementary alignment"}


mate_mapped_same_chr, mate_mapped_dif_chr, mate_mapped_dif_chr_a5 = 0,0,0
unmapped, paired, read1, read2, properly, duplicate, total = 0, 0, 0, 0, 0, 0, 0

bamfile = HTSeq.BAM_Reader(bam)
for almnt in bamfile:
	if almnt.aligned:
		if almnt.flag & 0x900 == 0:
			total += 1
			if almnt.flag & 0x400 != 0:
				duplicate += 1
			if almnt.mate_aligned:
				paired +=1
				if almnt.proper_pair:
					properly += 1
				if almnt.iv.chrom == almnt.mate_start.chrom:
					mate_mapped_same_chr += 1
				else:
					mate_mapped_dif_chr += 1
					if almnt.aQual >= 5:
						mate_mapped_dif_chr_a5 += 1
			#if almnt.pe_which == 'first':
			#	read1 += 1
			#if almnt.pe_which == 'second':
			#	read2 += 1
	else:
		unmapped += 1
		total += 1
					
print 'Total:\t%d (100%%)' % total
print 'Duplicate (rate of mapped):\t%d (%.2f%%)' % (duplicate, (duplicate*100.0)/(total-unmapped))
print 'Mapped:\t%d (%.2f%%)' % (total-unmapped, (total-unmapped)*100.0/total)
print 'Properly mapped:\t%d (%.2f%%)' % (properly, (properly*100.0)/total)
print 'PE mapped:\t%d (%.2f%%)' % (paired, (paired*100.0)/total)
print 'SE mapped:\t%d (%.2f%%)' % (total-paired, ((total-paired)*100.0)/total)
print 'With mate mapped to a different chr:\t%d (%.2f%%)' % (mate_mapped_dif_chr, (mate_mapped_dif_chr*100.0)/total)
print 'With mate mapped to a different chr ((mapQ>=5)):\t%d (%.2f%%)' % (mate_mapped_dif_chr_a5,(mate_mapped_dif_chr_a5*100.0)/total)

"""
job = 'picard_rmdupBam_%s_e' % sample
flist = os.listdir(logdir)
f = {}
if len(flist) == 1:
	tfile = os.path.join(logdir,flist[0])
else:
	for i in flist:
		if job in i:
			f[int(i[len(job):-4])] = i
	tfile = os.path.join(logdir,f[sorted(f.keys())[-1]])

q = re.compile(r'\s+(\d+)\s+')
f = open(tfile,'r')
for eachLine in f:
	if 'Marking' in eachLine:
		m = q.findall(eachLine)
		if m:
			x = m[0]
		else:
			x = 0
duplicate = int(x)
f.close()

print '%d + 0 in total (QC-passed reads + QC-failed reads)' % total
print '%d + 0 duplicates' % duplicate
print '%d + 0 mapped (%.2f%%:-nan%%)' % (total-unmapped, float(total-unmapped)/total)
print '%d + 0 paired in sequencing' % paired
print '%d + 0 read1' % read1
print '%d + 0 read2' % read2
print '%d + 0 properly paired (%.2f%%:-nan%%)' % (properly,float(properly)/total)
print '%d + 0 with itself and mate mapped' % (mate_mapped_dif_chr + mate_mapped_same_chr)
print '%d + 0 singletons (%.2f%%:-nan%%)' % (total-paired,float(total-paired)/total)
print '%d + 0 with mate mapped to a different chr' % (mate_mapped_dif_chr)
print '%d + 0 with mate mapped to a different chr (mapQ>=5)' % mate_mapped_dif_chr_a5

"""
