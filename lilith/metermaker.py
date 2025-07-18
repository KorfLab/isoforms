import files
import json
import sys



exonfile = sys.argv[1]
ksize = int(sys.argv[2])

kmers = {}
totalks = 0

with files.getfp(exonfile) as fp:
	for line in fp:
		l = line.strip()

		for i in range(len(l) - ksize +1):
			kmer = l[i:i+ksize]

			if kmer not in kmers:
				kmers[kmer] = 0

			kmers[kmer] += 1
			totalks += 1

for km in kmers.keys():
	kmers[km] = kmers[km] / totalks



print(json.dumps(kmers, sort_keys=True, indent=4))