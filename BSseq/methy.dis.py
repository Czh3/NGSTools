import sys

## calculate every C sites' methylation level to plot a bar_plot
# example: CpGs_methy_dis_global.png

#USAGE
# argv1: input bismark CX report
# argv2: output prefix

Cmethy = {}
CGmethy = {}
CHGmethy = {}
CHHmethy = {}


def addOne(h, k):
	if h.has_key(k):
		h[k] += 1
	else:
		h[k] = 1



for i in open(sys.argv[1]):
	j = i.split("\t")
	d = int(j[3]) + int(j[4])
	if d < 5:
		continue
	me = float(j[3])/d
	#me = '%.2f' % (me-0.005)
	me = '%.10f' % (me)
	me = me[:4]
	addOne(Cmethy, me)

	if j[5] == 'CG':
		addOne(CGmethy, me)
	elif j[5] == 'CHG':
		addOne(CHGmethy, me)
	elif j[5] == 'CHH':
		addOne(CHHmethy, me)

prefix = sys.argv[2].strip()


def output(ofile, hash):
	o = open(ofile, 'w')
	for i in sorted(hash.keys()):
		o.write(i +'\t'+ str(hash[i])+'\n')

output(prefix+".C", Cmethy)
output(prefix+".CG", CGmethy)
output(prefix+".CHG", CHGmethy)
output(prefix+".CHH", CHHmethy)



