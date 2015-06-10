#!/usr/bin/env python


# A windows based DMRs(sifferent metylaiton regions) calling program
# based on bismark output





import sys
from scipy import stats

def safeOpen(file, mode="r"):
	if file.endswith(".gz"):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)



def readWindows(fileHand1, fileHand2, ctype='CG', window=50, minDepth=10, site=2):

	methyList1 = []
	methyList2 = []
	avgDepth1 = 0
	avgDepth2 = 0
	
	line1 = fileHand1.next().split()
	line2 = fileHand2.next().split()

	windowsID = 0
	windowsInfo = {}

	while True:


		line1[3] = int(line1[3])
		line1[4] = int(line1[4])
		line2[3] = int(line2[3])
		line2[4] = int(line2[4])



		depth1 = line1[3] + line1[4]
		depth2 = line2[3] + line2[4]
		

		if depth1 >= minDepth and depth2 >= minDepth and line1[5] == ctype:

			methyList1.append(line1[3]/float(depth1))
			methyList2.append(line2[3]/float(depth2))
			avgDepth1 += depth1
			avgDepth2 += depth2

		newWindowsID = int(line1[1]) / window 
		if windowsID != newWindowsID:

			if len(methyList1) >= site:
				number = len(methyList1)
				windowsInfo[line1[0] +'a'+ str(windowsID)] = (number, avgDepth1, avgDepth2, methyList1, methyList2)

			methyList1 = []
			methyList2 = []
			avgDepth1 = 0
			avgDepth2 = 0
			windowsID =  newWindowsID

			
		try:
			line1 = fileHand1.next().split()
			line2 = fileHand2.next().split()
		except StopIteration:
			break

	return windowsInfo


def testDMR(data, site=5, window=50, diff=0.1, method="T_test", p_val=0.01, outfile="DMR.out"):
	'''read file from bismark output'''


	out = safeOpen(outfile, 'w')

	for k in data:
		chrom, windowID = k.split("a")
		windowID = int(windowID)

		keyWindowPlus1 = chrom +'a'+ str(windowID + 1)

		if not data.has_key(keyWindowPlus1):
			continue

		windowStart = window * windowID + 1
		windowEnd = window * (windowID + 2)

		# merge two step-windows
		number = data[k][0] + data[keyWindowPlus1][0]
		if number < site:
			continue
		depth1 = (data[k][1] + data[keyWindowPlus1][1])/number
		depth2 = (data[k][2] + data[keyWindowPlus1][2])/number

		methyListWindow1 = data[k][3] + data[keyWindowPlus1][3]
		methyListWindow2 = data[k][4] + data[keyWindowPlus1][4]
		avgMethy1 = sum(methyListWindow1) / number
		avgMethy2 = sum(methyListWindow2) / number

		if abs(avgMethy1 - avgMethy2) < 0.1:
			continue

		# DMR test
		if method == "T_test":
			t, p = stats.ttest_ind(methyListWindow1, methyListWindow2)
		elif method == "wilcoxon":
			t, p = stats.wilcoxon(methyListWindow1, methyListWindow2,correction=True)
		elif method == "fisher":
			all1 = depth1 * number
			me1 = all1 * avgMethy1
			unme1 = all1 = me1
			all2 = depth2 * number
			me2 = all2 * avgMethy2
			unme2 = all2 - me2
			odd, p = stats.fisher_exact([[me1, nume1], [me2, unme2]])

		if p <= p_val:
			outList = [chrom, windowStart, windowEnd, number, avgMethy1, avgMethy2, depth1, depth2, p]
			out.write("\t".join([str(i) for i in outList]) + '\n')


# main
if __name__ == '__main__':
	fileHand1 = safeOpen(sys.argv[1])
	fileHand2 = safeOpen(sys.argv[2])

	data = readWindows(fileHand1, fileHand2, window=50, minDepth=10, site=2)

	## background region
	# all sequenced region
	if True:
		back = open("background.bed", 'w')
		for k in data:
			chrom, windowID = k.split("a")
			windowID = int(windowID)
			windowStart = 50 * windowID + 1
			windowEnd = 50 * (windowID + 2)
			back.write(chrom+"\t"+str(windowStart)+'\t'+str(windowEnd)+'\n')



	testDMR(data, method="wilcoxon")



