
import os

cwd = os.getcwd()
samples_path = cwd + '/samples_filter_all/'
queryStr = '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
#queryStr = 'foobar'

for f in os.listdir(samples_path):
	fPATH = samples_path + f
	with open(fPATH, 'r') as currFile:
		if queryStr not in currFile.read():
			print(f)