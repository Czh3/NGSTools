import sys

d3 = open("/lustre/user/liclab/zhangc/proj/lingte/RRBS/bismark/reorder/C3.CX.txt")
d4 = open("/lustre/user/liclab/zhangc/proj/lingte/RRBS/bismark/reorder/D4.CX.txt")


while True:
	try:
		j3 = d3.next().split()
		j4 = d4.next().split()	
	except StopIteration:
		break

	j3[3] = int(j3[3])
	j3[4] = int(j3[4])
	j4[3] = int(j4[3])
	j4[4] = int(j4[4])
	c3 = j3[3]+j3[4] 
	c4 = j4[3]+j4[4]
	if c3 < 10 or c4 < 10:
		continue

	me3 = j3[3]/float(c3)
	me4 = j4[3]/float(c4)
	if me3 > 0.7 and me4 < 0.2:
		print '\t'.join([j3[0], j3[1], j3[1], 'C', 'T', str(me3), str(me4), j3[5]])
	elif me4 > 0.7 and me3 < 0.2:
		print '\t'.join([j3[0], j3[1], j3[1], 'C', 'T', str(me3), str(me4), j3[5]])


