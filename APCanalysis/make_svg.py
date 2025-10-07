import argparse

parser = argparse.ArgumentParser(
	description='create svg file from APC isoforms')
parser.add_argument('apc_gff')
parser.add_argument('wb_gff')

arg = parser.parse_args()


c = 0 
with open(arg.apc_gff, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		line = line.split('\t')
		if len(line) != 9: continue
		if line[2] == 'mRNA': continue
		iso_parent = line[8].split('=')[1]
		print(iso_parent)
		

		c += 1
		if c == 20: break
