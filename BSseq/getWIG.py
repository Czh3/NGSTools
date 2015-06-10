import sys

d3 = open("/lustre/user/liclab/zhangc/proj/lingte/RRBS/bismark/reorder/C3.CX.txt")
d4 = open("/lustre/user/liclab/zhangc/proj/lingte/RRBS/bismark/reorder/D4.CX.txt")
o3 = open("D3.wig", 'w')
o4 = open("D4.wig", 'w')

# header
o3.write('track type=wiggle_0 name="D3" description="methylation level mC/(mC+Un_mC)" visibility=full autoScale=off viewLimits=0.0:1.0 color=65,105,225 yLineMark=11.76 yLineOnOff=on priority=10\n')
o4.write('track type=wiggle_0 name="D4" description="methylation level mC/(mC+Un_mC)" visibility=full autoScale=off viewLimits=0.0:1.0 color=221,160,221 yLineMark=11.76 yLineOnOff=on priority=10\n')


chrom = "chr"

while True:
	try:
		j3 = d3.next().split()
		j4 = d4.next().split()	
	except StopIteration:
		break

	if j3[5] != "CG":
		continue


	j3[3] = int(j3[3])
	j3[4] = int(j3[4])
	j4[3] = int(j4[3])
	j4[4] = int(j4[4])
	c3 = j3[3]+j3[4] 
	c4 = j4[3]+j4[4]

	if c3 < 10 or c4 < 10:
		continue

	if j3[0] != chrom:
		o3.write('variableStep chrom=' + j3[0] + '\n')
		o4.write('variableStep chrom=' + j3[0] + '\n')
		chrom = j3[0]

	me3 = j3[3]/float(c3)
	me4 = j4[3]/float(c4)
	o3.write(j3[1] + '\t' + str(me3) + '\n')
	o4.write(j4[1] + '\t' + str(me4) + '\n')


