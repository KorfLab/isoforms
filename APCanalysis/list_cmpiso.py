import argparse
import glob
import isoform


parser = argparse.ArgumentParser(description=
	'run cmpiso on all APC genes and organize results')
parser.add_argument('APCisos', nargs='+', type=str, metavar='<file>', 
	help='list of directories with APC iso results')

args = parser.parse_args()

for dpath in args.APCisos:

	if dpath.endswith('/'): dpath = dpath
	else: dpath = f'{dpath}/'
	
	print(dpath)
	for fpath in glob.glob(f'{dpath}*'):
		print(fpath)


