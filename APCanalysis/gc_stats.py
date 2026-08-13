import argparse
import glob
import json

parser = argparse.ArgumentParser(description='get gc content summary '
					'statistics into csv format for jupyter notebook')
parser.add_argument('gcc_res', help='directory with gc content '
					'results')
					
args = parser.parse_args()

if not args.gcc_res.endswith('/'):
	args.gcc_res = f'{args.gcc_res}/'

# ce.1.2, first field in isodif ranks

avg_gc = {}
for file in glob.glob(f'{args.gcc_res}*'):
	gid = '.'.join(file.split('/')[-1].split('.')[:-2])
	with open(file, 'r') as fp:
		data = json.load(fp)
		for item1 in data.items():
			cds_gc = 0
			cds_n = 0
			int_gc = 0
			int_n = 0
			for item2 in item1[1].items():
				if not item2[1]: continue
				for item3 in item2[1].items():
					if item2[0] == 'cds':
						cds_gc += item3[1]
						cds_n += 1
					if item2[0] == 'intron':
						int_gc += item3[1]
						int_n += 1
			if int_n == 0:
				cds_avg = round(cds_gc/cds_n, 2)
				int_avg = None
			else:
				cds_avg = round(cds_gc/cds_n, 2)
				int_avg = round(int_gc/int_n, 2)
			print(gid, cds_avg, int_avg)
			if gid not in avg_gc:
				avg_gc[gid] = {'cds': [cds_avg], 'int': [int_avg]}
			else:
				avg_gc[gid]['cds'].append(cds_avg)
				avg_gc[gid]['int'].append(int_avg)
				
avg_gc2 = {}
for item1 in avg_gc.items():
	for item2 in item1[1].items():
		vals = []
		item_len = 0
		for x in item2[1]:
			if x:
				vals.append(x)
				item_len += 1
		if item1[0] not in avg_gc2:
			avg_gc2[item1[0]] = {}
			avg_gc2[item1[0]][item2[0]] = round(sum(vals)/item_len, 2)
		else:
			avg_gc2[item1[0]][item2[0]] = round(sum(vals)/item_len, 2)
		
for item in avg_gc2.items():
	print(item)


