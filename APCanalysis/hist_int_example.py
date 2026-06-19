# create a csv file to import into jupyter notebooks
# example histogram intersection between APC and RNA

import glob
import argparse
import isoform
import copy

parser = argparse.ArgumentParser(description='create histogram intersection '
								'csv file for figures')
parser.add_argument('APC', help='APC generated gff files')
parser.add_argument('RNA', help='RNA-seq database gff files')

args = parser.parse_args()

if not args.APC.endswith('/'):
	args.APC = f'{args.APC}/'
	
if not args.RNA.endswith('/'):
	args.RNA = f'{args.RNA}/'
	
# create identical length histograms
def add_zeroes(introns1, introns2):
	
	i1 = copy.deepcopy(introns1)
	i2 = copy.deepcopy(introns2)
	
	for i in i1:
		if i not in i2: i2[i] = 0
	
	for i in i2:
		if i not in i1: i1[i] = 0
		
	i2_sort = {}
	for i in i1:
		i2_sort[i] = i2[i]
		
	return i1, i2_sort
	
hists = {}
for apc_path in glob.glob(f'{args.APC}*'):
	gID = apc_path.split('.')[-3]
	rna_path = f'{args.RNA}ch.{gID}.gff3'
	apc_ints = isoform.get_introns(apc_path)
	rna_ints = isoform.get_introns(rna_path)
	i1, i2 = add_zeroes(apc_ints, rna_ints)
	hists[gID] = [[i1], [i2]]
	for i in i1.items():
		print(i)
	break

for item in hists.items():
	print(item)
	break
	

