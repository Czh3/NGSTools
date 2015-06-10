import sys

d3 = "/lustre/user/liclab/zhangc/proj/lingte/RRBS/bismark/reorder/C3.CX.txt"
d4 = "/lustre/user/liclab/zhangc/proj/lingte/RRBS/bismark/reorder/D4.CX.txt"

sites = []
for i in open("../DMR_0.3"):
	j = i.split()
	for k in xrange(int(j[1]), int(j[2])+1):
		sites.append(j[0]+'a'+str(k))

def get_output(sample, sites):
	#o = open(out, 'w')
	for i in open(sample):
		j = i.split()

		if j[5] != 'CG':
			continue 

		d = int(j[3])+int(j[4])
		if d < 5:
			continue

		if j[0]+'a'+j[1] in sites:
			me = float(j[3])/d
			me = '%.8f' % (me)
			me = me[:4]
			print j[0]+'\t'+j[1]+'\t'+me
get_output(sys.argv[1], sites)


