import argparse
import glob
import isoform


parser = argparse.ArgumentParser(description=
	'run cmpiso on all APC genes and organize results')
parser.add_argument('APCisos', nargs='+', type=str, metavar='<file>', 
	help='list of directories with APC iso results, ' 
	'input individual directory paths')

args = parser.parse_args()

paths2gffs = {}
for dpath in args.APCisos:

	if dpath.endswith('/'): dpath = dpath
	else: dpath = f'{dpath}/'
	
	dir_name = dpath.split('/')[-2]
	paths2gffs[dir_name] = []
	for fpath in glob.glob(f'{dpath}*'):
		paths2gffs[dir_name].append(fpath)
		
for item in paths2gffs.items():
	print(item)

