#!/PUBLIC/software/public/System/Python-2.7.6/bin/python
import pysam
import argparse
import re
import os

parser = argparse.ArgumentParser(description="samtools flagstat statistics")
parser.add_argument('--bam',help="the input bam file",type=str)
#parser.add_argument('--logdir',help="the path for .e & .o file")
#parser.add_argument('--sample',help="sample name")
#parser.add_argument('--txt',help="out file",type=str)

argv = vars(parser.parse_args())
#logdir = argv['logdir']
bam = argv['bam']
#sample = argv['sample']
#txt = argv['txt']

flags = {1 : "template having multiple segments in sequencing",
         2 : "each segment properly aligned according to the aligner",
         3 : "segment unmapped",
         4 : "next segment in the template unmapped",
         5 : "SEQ being reverse complemented",
         6 : "SEQ of the next segment in the template being reversed",
         7 : "the first segment in the template",
         8 : "the last segment in the template",
         9 : "secondary alignment",
         10 : "not passing quality controls",
         11 : "PCR or optical duplicate",
         12 : "supplementary alignment"}


infile = pysam.Samfile(bam,'rb')
#reads = {}

mate_mapped_same_chr, mate_mapped_dif_chr, mate_mapped_dif_chr_a5 = 0,0,0
unmapped, read1, read2, properly, duplicate, total, PE_unmapped = 0, 0, 0, 0, 0, 0, 0
SE_mapped = 0
for n in infile:
    if n.flag & 0x900 != 0:continue
    total += 1
    if n.rname != '*' and n.rnext != '*':
        if n.rname != n.rnext:
            mate_mapped_dif_chr += 1
            if n.mapq >= 5:
                mate_mapped_dif_chr_a5 += 1
        else:
            mate_mapped_same_chr += 1
    if n.is_duplicate:
        duplicate += 1
    if n.is_proper_pair:
        properly += 1
#    if n.is_read1:
#        read1 += 1
#    if n.is_read2:
#        read2 += 1
    if n.is_unmapped:
        unmapped += 1
        if n.mate_is_unmapped:
            PE_unmapped += 1
    elif n.mate_is_unmapped:
        SE_mapped += 2


print 'Total:\t%d (100%%)' % total
print 'Duplicate:\t%d (%.2f%%)' % (duplicate, (duplicate*100.0)/(total-unmapped))
print 'Mapped:\t%d (%.2f%%)' % (total-unmapped, ((total-unmapped)*100.0)/total)
print 'Properly mapped:\t%d (%.2f%%)' % (properly, (properly*100.0)/total)
print 'PE mapped:\t%d (%.2f%%)' % (total-SE_mapped-PE_unmapped, ((total-SE_mapped-PE_unmapped)*100.0)/total)
print 'SE mapped:\t%d (%.2f%%)' % (SE_mapped, (SE_mapped*100.0)/total)
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
A = open(tfile,'r')
for n in A:
    if 'Marking' in n:
        m = q.findall(n)
        if m:
            x = q.findall(n)[0]
        else:x = 0
duplicate = int(x)

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
