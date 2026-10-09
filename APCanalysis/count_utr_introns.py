import argparse
import glob

parser = argparse.ArgumentParser(description='count number of APC '
	'predictions that have introns in 5\' or 3\' utr')
parser.add_argument('smallgenes')
parser.add_argument('--apc_results', type=str, required=False, nargs=4, 
	help='directories with APC results')
parser.add_argument('--top_isos', type=int, required=False, default=1,
	help='search only the top %(default)s isos')
	
args = parser.parse_args()

longest_cdss = {}
for gpath in glob.glob(f'{args.smallgenes}*gff3'):
	gid = '.'.join(gpath.split('/')[-1].split('.')[:-1])
	rna_cds = {}
	coor_list = []
	with open(gpath, 'r') as fp:
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[1] == 'WormBase' and line[2] == 'CDS':
				cds = [int(line[3]), int(line[4])]
				for c in cds:
					coor_list.append(c)
	longest_cdss[gid] = [min(coor_list), max(coor_list)]

def get_isos(apc_gff):
	
	isoforms = {}
	count = 0
	with open(apc_gff, 'r') as agfp:
		for line in agfp:
			line = line.rstrip()
			line = line.split('\t')
			if len(line) != 9: continue
			if line[2] == 'mRNA' or line[2] == 'gene': continue
			feature = [line[2], int(line[3]), int(line[4]), line[5]]
			iso_parent = line[8].split('=')[1]
			if iso_parent not in isoforms:
				isoforms[iso_parent] = [feature]
			else:
				isoforms[iso_parent].append(feature)
				
	return isoforms

by_apc_type = {}
for res in args.apc_results:
	apc_type = res.split('/')[-2]
	int_in_utrs = {}
	for file in glob.glob(f'{res}*'):
		isos = get_isos(file)
		icount = 0
		gid = '.'.join(file.split('/')[-1].split('.')[:3])
		int_in_utrs[gid] = {}
		for item in isos.items():
			if icount < args.top_isos:
				gid = item[0].split('-')[1]
				int_coors = []
				for feature in item[1]:
					if feature[0] == 'intron':
						int_coors.append(feature[1])
						int_coors.append(feature[2])
			else: break
			icount += 1
			cl, cr = longest_cdss[gid][0], longest_cdss[gid][1]
			il, ir = min(int_coors), max(int_coors)
			if il < cl or ir > cr:
				# over utr
				int_in_utrs[gid][icount] = True
			else:
				# not over utr
				int_in_utrs[gid][icount] = False
				
	by_apc_type[apc_type] = int_in_utrs
	
# summarize results
utr_counts = {}
for apc in by_apc_type.items():
	for gene in apc[1].items():
		# just do top isoform, no need for args.top_isos
		t_or_f = gene[1][1]
		if apc[0] not in utr_counts:
			utr_counts[apc[0]] = {True: 0, False: 0}
			utr_counts[apc[0]][t_or_f] += 1
		else:
			utr_counts[apc[0]][t_or_f] += 1
			
# print this for summary
#print(utr_counts)

# csv format with gene ids for pandas
print(f'gene_id,top_int_utr,apc_type')
for apc in by_apc_type.items():
	if apc[0] == 'APCisos_base':
		atype = 'base'
	if apc[0] == 'APCisos_nmd':
		atype = 'nmd'
	if apc[0] == 'APCisos_optiso':
		atype = 'optiso'
	if apc[0] == 'APCisos_optiso_nmd':
		atype = 'optiso_nmd'
	for gene in apc[1].items():
		print(f'{gene[0]},{gene[1][1]},{atype}')




