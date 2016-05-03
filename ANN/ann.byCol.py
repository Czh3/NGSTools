import sys

def ann_col(input_file, db_file, col_input, col_db, col_db_ann):
	out_put = []
	db = {}
	for i in open(db_file):
		j = i.split()
		ann = ''
		for k in col_db_ann.split(","):
			ann += "\t"+j[int(k)]
		db[j[int(col_db)]] = ann.lstrip("\t")

	for i in open(input_file):
		j = i.split()
		if j[int(col_input)] in db:
			ann = db[j[int(col_input)]]
		else:
			ann = "NULL"
		out_put.append(i.strip() + '\t' + ann + '\n')
	return out_put

if len(sys.argv) != 6:
	sys.exit("argv: input_file, db_file, col_input, col_db, col_db_ann")
else:
	out = ann_col(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
	for i in out:
		print i.strip()
