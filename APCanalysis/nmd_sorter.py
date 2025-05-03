# quick script to sort NMD teste results

import sys
import csv

file = sys.argv[1]

genes = {}
with open(file, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split(',')
		if line[0] not in genes:
			genes[line[0]] = [line[1:4]] 
		else:
			genes[line[0]].append(line[1:4])
# order of distances is c, k, x, y
# start with sorting manhattan distance c
cdist = {}
kdist = {}
xdist = {}
ydist = {}
for g in genes:
	cdist[g] = genes[g][0]
	kdist[g] = genes[g][1]
	xdist[g] = genes[g][2]
	ydist[g] = genes[g][3]
	
# if deltas are 0, no NMD? look at individual genes to verify

with open('mdist_sorted.csv', 'w') as csvfile:
	writer = csv.writer(csvfile)
	for item in sorted(cdist.items(), key=lambda item: float(item[1][0])):
		info = item[1]
		info.insert(0, item[0])
		writer.writerow(info)

# manhattan distance is probably all we need

# quick stats
# distribution of distances, does intron number affect predictability?
# average distance per intron count normalized, many 1s but less of others
# how do predictions cluster? what is driving this trend? if there is one
# median/mean distance
# Does taking NMD into account significantly affect predictions





