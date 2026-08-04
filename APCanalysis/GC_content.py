import argparse
import glob

parser = argparse.ArgumentParser(description='get GC content data from '
	'smallgenes and APC')
parser.add_argument('smallgenes')
parser.add_argument('APC')
	
args = parser.parse_args()

if not args.smallgenes.endswith('/'):
	args.smallgenes = f'{args.smallgenes}/'
	
if not args.APC.endswith('/'):
	args.APC = f'{args.APC}/'

sg_paths = {}
for fasta in glob.glob(f'{args.smallgenes}*.fa'):
		gid = '.'.join(fasta.split('/')[-1].split('.')[:-1])
		gff = f'{args.smallgenes}{gid}.gff3'
		gpaths[gid] = [fasta, gff]
		

