import itertools
import math
import os
import subprocess
import sys

def readfile(filename):
	perf = None
	details = {}
	with open(filename) as fp:
		perf = fp.readline().rstrip()
		for line in fp:
			beg, end, p1, p2 = line.split()
			beg = int(beg)
			end = int(end)
			p1 = float(p1)
			p2 = float(p2)
			details[(beg,end)] = (p1, p2)
	return perf, details

def compare_outputs(file1, file2):
	p1, d1 = readfile(file1)
	p2, d2 = readfile(file2)

	report = []
	if float(p1) == float(p2): return report
#	if math.isclose(float(p1), float(p2), abs_tol=1e-6): return report

	report.append(f'performance mismatch {p1} {p2}')
	for loc in d1:
		if loc not in d2: report.append(f'{loc} missing in {file2}')
	for loc in d2:
		if loc not in d1: report.append(f'{loc} missing in {file1}')
	for loc in d1:
		if loc not in d2: continue
		a1, a2 = d1[loc]
		b1, b2 = d2[loc]
		if a1 != b1 or a2 != b2:
			report.append(f'calculation mismatch {d1[loc]} vs. {d2[loc]}')
	return report

# genes in the test set
genes = [
	'ch.13301',
	'ch.19857',
	'ch.2110',
	'ch.9395',
]

# programs to check (just pairs)
p1 = 'geniso'    # original python implementation
p2 = 'isoformer' # genomikon C implementation

# base parameters not changed (for now)
base = '--min_exon 25 --min_intron 35 --max_splice 3 --flank 99 --icost 20'

# combinatoric parameter swapping -> parameter_sets
opts = [
	"--apwm models/acc.pwm",
	"--dpwm models/don.pwm",
	"--emm models/exon.mm",
	#"--imm models/intron.mm", # imm requires don, acc so can't run solo
	"--elen models/exon.len",
	"--ilen models/intron.len",
]
parameter_sets = []
for k in range(len(opts) +1):
	for t in itertools.combinations(opts, k):
		parameter_sets.append(' '.join(t))

# testing loop
for pset in parameter_sets:
	for gene in genes:
		ff = f'apc/{gene}.fa'
		gff = f'apc/{gene}.gff3'
		pid = os.getpid()
		out1 = f'/tmp/{pid}.1.gff'
		out2 = f'/tmp/{pid}.2.gff'
		cmp1 = f'/tmp/{pid}.1.cmp'
		cmp2 = f'/tmp/{pid}.2.cmp'

		# run both programs
		subprocess.run(f'{p1} {ff} {pset} > {out1}', shell=True)
		subprocess.run(f'{p2} {ff} {pset} > {out2}', shell=True)

		# compare program output to standard set
		subprocess.run(f'cmpiso {gff} {out1} > {cmp1}', shell=True)
		subprocess.run(f'cmpiso {gff} {out2} > {cmp2}', shell=True)

		# compare to each other
		report = compare_outputs(cmp1, cmp2)

		# clean up
		os.remove(out1)
		os.remove(out2)
		os.remove(cmp1)
		os.remove(cmp2)

		if len(report) == 0:
			print(f'{gene} no errors found with parameters: "{pset}"')
		else:
			print('\nDifferences Found')
			print(f'Parameters: "{pset}"')
			for line in report: print(line)
			sys.exit('aborting')
