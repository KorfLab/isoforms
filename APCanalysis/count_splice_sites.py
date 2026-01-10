import argparse
import glob
import isoform

parser = argparse.ArgumentParser(description=
	'count number of splice sites in APC genes')
parser.add_argument('smallgenes')

arg = parser.parse_args()

flank = 99
minex = 25

site_counts = {}
for fpath in glob.glob(f'{arg.smallgenes}/*.fa'):
	for name, seq in isoform.read_fasta(fpath):
		dons, accs = isoform.gtag_sites(seq, flank, minex)
		site_counts[name] = [len(dons), len(accs), len(seq)]

print(f'dons,accs,seq_len')
for item in site_counts.items():
	print(f'{item[1][0]},{item[1][1]},{item[1][2]}')
	

